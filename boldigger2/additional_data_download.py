import more_itertools, requests_html, datetime, time, json
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from boldigger2.exceptions import APIOverload
from boldigger2.exceptions import ProxyNotWorking
from requests.exceptions import ReadTimeout
from requests.exceptions import ChunkedEncodingError
from requests.exceptions import ProxyError
from requests.exceptions import ConnectTimeout
from fp.fp import FreeProxy
from fp.errors import FreeProxyException
from fp.errors import FreeProxyException
from json.decoder import JSONDecodeError
from urllib3.exceptions import ReadTimeoutError


# function to sort the hdf dataframe according to the order in the fasta file
# also removes duplicate entries from malformed requests in the previous step
def read_and_order(fasta_path, hdf_name_top_100_hits, read_fasta):
    # read in the fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # the keys are in perfect order, use those to order the hdf file
    sorter = {name: idx for idx, name in enumerate(fasta_dict.keys())}

    # read the hdf data that needs to be sorted
    top_100_hits = pd.read_hdf(hdf_name_top_100_hits, key="top_100_hits_unsorted")

    # remove duplicate entries from malformed responses here
    # generate a unique id
    top_100_hits["unique_id"] = (
        top_100_hits["ID"] + top_100_hits["database"] + top_100_hits["request_date"]
    )

    # create a filter
    filter = top_100_hits.drop_duplicates(subset=["ID", "database", "request_date"])[
        ["ID", "database", "request_date"]
    ].reset_index(drop=True)

    # drop the duplicates
    filter_no_duplicates = filter.drop_duplicates(subset=["ID", "database"]).index

    # select only the entries that are no duplicates or in case of duplication select the 1st answer
    filter = filter.iloc[filter_no_duplicates]
    filter["unique_id"] = filter["ID"] + filter["database"] + filter["request_date"]

    # only select those hits from the top 100 hits, drop the unique id after
    top_100_hits = top_100_hits.loc[top_100_hits["unique_id"].isin(filter["unique_id"])]
    top_100_hits = top_100_hits.drop(labels=["unique_id"], axis=1)

    # add a sorter helper column
    top_100_hits["sorter"] = top_100_hits["ID"].map(sorter)

    # name the index for sorting
    top_100_hits.index.name = "index"

    # sort the results, remove the sorter, reset the index
    top_100_hits = (
        top_100_hits.sort_values(
            by=["sorter", "database", "index"],
            ascending=[True, False, True],
        )
        .drop(labels=["sorter"], axis=1)
        .reset_index(drop=True)
    )

    # add the sorted dataframe to the original hdf storage
    # add the results to the hdf storage
    # set size limits for the columns
    item_sizes = {
        "ID": 100,
        "Phylum": 80,
        "Class": 80,
        "Order": 80,
        "Family": 80,
        "Genus": 80,
        "Species": 80,
        "Subspecies": 80,
        "Status": 15,
        "Process_ID": 25,
        "database": 20,
        "request_date": 30,
    }

    # only have to write the results once
    try:
        pd.read_hdf(hdf_name_top_100_hits, key="top_100_hits_sorted")
        # give user output
        print(
            "{}: Hits are already ordered from a previous run.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
    except KeyError:
        # append results to hdf
        with pd.HDFStore(
            hdf_name_top_100_hits, mode="a", complib="blosc:blosclz", complevel=9
        ) as hdf_output:
            hdf_output.append(
                "top_100_hits_sorted",
                top_100_hits,
                format="t",
                data_columns=True,
                min_itemsize=item_sizes,
                complib="blosc:blosclz",
                complevel=9,
            )

        print(
            "{}: Hits ordered successfully.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

    # drop process IDs that are empty
    with pd.option_context("future.no_silent_downcasting", True):
        process_ids = top_100_hits["Process_ID"].replace("", np.nan).dropna()

    # return process ids
    return top_100_hits, process_ids


# function to check if for any of the process IDs the additional data has already been downloaded
# also removes duplicate entries from the process ids to prepare the download
# split the data up into managable batches of 100 processids
def data_already_downloaded(process_ids, hdf_name_top_100_hits):
    # check if the hdf storage already contains additional data
    try:
        already_downloaded = pd.read_hdf(hdf_name_top_100_hits, key="additional_data")
        already_downloaded = already_downloaded["processid"]
    except:
        already_downloaded = []

    # filter all ids that are already have been downloaded
    filter = set(already_downloaded)
    process_ids_to_download = [id for id in process_ids if id not in filter]
    process_ids_to_download = pd.Series(process_ids_to_download).unique()

    # return in batches of 100 ids
    return more_itertools.chunked(process_ids_to_download, 100)


# function to generate an API download link from a batch of process ids
def generate_download_link(process_id_batch):
    url = "http://www.boldsystems.org/index.php/API_Public/specimen?ids={}&format=json".format(
        "|".join(process_id_batch)
    )

    return url


# function to parse the xml returned by the BOLD api into a dataframe
def json_response_to_dataframe(response, process_id_batch, hdf_name_top_100_hits):
    if "You have exceeded" in response.text:
        raise APIOverload
    if "REMOTE_ADDR" in response.text:
        raise ProxyNotWorking

    # load the json response
    response_data = json.loads(response.text)["bold_records"]["records"]

    # collect all results of one process id batch here
    process_id_batch_results = []

    # collect the data from json
    for process_id in process_id_batch:
        id_values = response_data.get(process_id, {})
        # get all values that are not nested first
        process_id = id_values.get("processid", "")
        record_id = id_values.get("record_id", "")
        bin_uri = id_values.get("bin_uri", "")
        # get data from the specimen identifiers
        specimen_identifiers = id_values.get("specimen_identifiers", {})
        institution_storing = specimen_identifiers.get("institution_storing", "")
        # get data from the specimen description
        specimen_desc = id_values.get("specimen_desc", {})
        sex = specimen_desc.get("sex", "")
        lifestage = specimen_desc.get("lifestage", "")
        # get data from the collection event
        collection_event = id_values.get("collection_event", {})
        country = collection_event.get("country", "")
        # get data from the taxonomy
        taxonomy = id_values.get("taxonomy", {})
        identification_provided_by = taxonomy.get("identification_provided_by", "")
        identification_method = taxonomy.get("identification_method", "")

        # all values together form one new line for the output dataframe
        process_id_batch_results.append(
            [
                process_id,
                record_id,
                bin_uri,
                institution_storing,
                sex,
                lifestage,
                country,
                identification_provided_by,
                identification_method,
            ]
        )

    # generate a dataframe
    process_id_batch_results = pd.DataFrame(
        process_id_batch_results,
        columns=[
            "processid",
            "record_id",
            "bin_uri",
            "institution_storing",
            "sex",
            "lifestage",
            "country",
            "identification_provided_by",
            "identification_method",
        ],
    )

    # append the data to a new key in the hdf storage
    item_sizes = {
        "processid": 30,
        "record_id": 10,
        "bin_uri": 15,
        "institution_storing": 150,
        "sex": 8,
        "lifestage": 80,
        "country": 80,
        "identification_provided_by": 80,
        "identification_method": 150,
    }

    # append results to hdf
    with pd.HDFStore(
        hdf_name_top_100_hits, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            "additional_data",
            process_id_batch_results,
            format="t",
            data_columns=True,
            min_itemsize=item_sizes,
            complib="blosc:blosclz",
            complevel=9,
        )


# function to generate a new proxy
def fresh_proxy():
    while True:
        try:
            proxy = FreeProxy(https=True, rand=True).get()
            tqdm.write(
                "{}: Proxy set to {}.".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), proxy
                )
            )
            return proxy
        except FreeProxyException:
            tqdm.write(
                "{}: No proxy available, waiting.".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )
            time.sleep(60)
            return ""


def download_data(process_ids_to_download, hdf_name_top_100_hits):
    with requests_html.HTMLSession() as session:
        # create a proxy for later
        proxy = ""
        for id_batch in tqdm(
            list(process_ids_to_download), desc="Downloading additional data"
        ):
            # create a url for the id batch
            url = generate_download_link(id_batch)

            # define a timeout counter so set a fresh proxy from time to time
            timeout_counter = 0
            # run until getting a valid response
            while True:
                try:
                    # as long as the original IP is working, use this one
                    if not proxy:
                        response = session.get(url, timeout=60)
                    else:
                        response = session.get(
                            url, timeout=60, proxies={"http": proxy, "https": proxy}
                        )
                    # parse the response
                    json_response_to_dataframe(
                        response, id_batch, hdf_name_top_100_hits
                    )
                    break
                except (
                    ReadTimeout,
                    ConnectTimeout,
                    ChunkedEncodingError,
                    ConnectionError,
                    ReadTimeoutError,
                ):
                    tqdm.write(
                        "{}: Read timed out, retrying.".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )
                    timeout_counter += 1
                    if timeout_counter >= 10:
                        proxy = fresh_proxy()
                        timeout_counter = 0
                    continue
                except APIOverload:
                    tqdm.write(
                        "{}: API overloaded. Switching proxy.".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )
                    # set ip adress via a proxy
                    proxy = fresh_proxy()
                    continue
                except (ProxyError, ProxyNotWorking):
                    # set ip adress via a proxy
                    proxy = fresh_proxy()
                    continue
                except JSONDecodeError:
                    tqdm.write(
                        "{}: Malformed response. Switching proxy.".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )
                    # set ip adress via a proxy
                    proxy = fresh_proxy()
                    continue


# function to add the additional data to the top 100 hits
def add_additional_data(hdf_name_top_100_hits, top_100_hits, process_ids):
    # load the additional data downloaded
    additional_data = pd.read_hdf(hdf_name_top_100_hits, key="additional_data")
    # transform additional data to dict, retain the column names
    additional_data = additional_data.to_dict("tight")
    column_names = additional_data["columns"][1:]
    additional_data = additional_data["data"]

    # parse the additional data into a dict in the form of process_id : [data fields] to rebuild the dataframe
    additional_data = [(record[0], record[1:]) for record in additional_data]
    additional_data = {record: data for record, data in additional_data}

    # recreate to dataframe with duplicate values
    additional_data = pd.DataFrame(
        data=[additional_data[record] for record in process_ids],
        columns=column_names,
        index=process_ids.index,
    )

    # add specimen page links
    additional_data["specimen_page_url"] = [
        "http://www.v4.boldsystems.org/index.php/MAS_DataRetrieval_OpenSpecimen?selectedrecordid={}".format(
            record_id
        )
        for record_id in additional_data["record_id"]
    ]

    # merge the additional data and the top 100 hits on index
    top_100_hits = pd.concat([top_100_hits, additional_data], axis=1)

    # add the top 100 hits with additional data to the hdf storage
    # in this case we can infer the size of the columns since we won't append to this file anymore
    with pd.HDFStore(
        hdf_name_top_100_hits, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            "top_100_hits_additional_data",
            top_100_hits,
            format="t",
            data_columns=True,
            complib="blosc:blosclz",
            complevel=9,
        )


def excel_converter(hdf_name_top_100_hits):
    top_100_hits = pd.read_hdf(
        hdf_name_top_100_hits, key="top_100_hits_additional_data"
    )

    # split the dataframe by 1.000.000 entries
    idx_parts = more_itertools.chunked(top_100_hits.index, 1000000)

    # generate an excel savename
    excel_savename = Path(hdf_name_top_100_hits).with_suffix("").with_suffix("")

    for idx, idx_part in enumerate(idx_parts):
        excel_savename = "{}_part_{}.xlsx".format(excel_savename, idx)
        top_100_hits.iloc[idx_part].to_excel(excel_savename, index=False)


# function to check if the additional data has already been downloaded
# download can be skipped if that is the case --> returns True
def additional_data_present(hdf_name_top_100_hits):
    try:
        pd.read_hdf(hdf_name_top_100_hits, key="top_100_hits_additional_data")
        # give user output
        print(
            "{}: Additional data has already been downloaded.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        return True
    except KeyError:
        # if no additional data can be found return False
        return False


# main function to run the additional data download
def main(fasta_path, hdf_name_top_100_hits, read_fasta):
    # give user output
    print(
        "{}: Trying to order the top 100 hits.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # read and sort the hdf file according to the order in the fasta file
    top_100_hits, process_ids = read_and_order(
        fasta_path, hdf_name_top_100_hits, read_fasta
    )

    # check if some of the ids have already been downloaded
    process_ids_to_download = data_already_downloaded(
        process_ids, hdf_name_top_100_hits
    )

    # skip the download if the data is already present
    if not additional_data_present(hdf_name_top_100_hits):
        # download the data
        download_data(process_ids_to_download, hdf_name_top_100_hits)

        # add the metadata to the top 100 hits, push to a new hdf table
        top_100_hits = add_additional_data(
            hdf_name_top_100_hits, top_100_hits, process_ids
        )

    # give user output
    print(
        "{}: Additional data successfully downloaded and saved.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    print(
        "{}: Saving top 100 hits to excel, this may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # run the excel converter in the end
    excel_converter(hdf_name_top_100_hits)


# run only if called as a toplevel script
if __name__ == "__main__":
    main()
