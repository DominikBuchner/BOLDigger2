import more_itertools, asyncio, requests_html, datetime, time, sys
import pandas as pd
import numpy as np
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm.asyncio import tqdm_asyncio
from pathlib import Path
from tqdm import tqdm
from boldigger2.exceptions import APIOverload
from requests.exceptions import ReadTimeout
from requests.exceptions import ConnectionError
from bs4 import BeautifulSoup as BSoup


# function to sort the hdf dataframe according to the order in the fasta file
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

    # override the hdf file with the new, sorted dataframe
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

    # drop process IDs that are empty
    with pd.option_context("future.no_silent_downcasting", True):
        process_ids = top_100_hits["Process_ID"].replace("", np.nan).dropna()

    # return process ids
    return top_100_hits, process_ids


# function to check if some of the additional data has already been downloaded
def check_already_downloaded(process_ids, hdf_name_top_100_hits):
    # try to see if additional data has already been downloaded partly
    try:
        already_downloaded = pd.read_hdf(hdf_name_top_100_hits, key="additional_data")
        already_downloaded = list(already_downloaded["processid"])
    except KeyError:
        already_downloaded = []

    # filter all ids that are already done
    process_ids = [id for id in process_ids if id not in already_downloaded]

    # return the remaining process ids
    return pd.Series(process_ids)


# function to return download links for the BOLD API from the process IDs
def generate_download_links(process_ids):
    # chunk the process IDs into batches of 100 ids remove duplicates, since we only need to query everything once
    process_ids = more_itertools.chunked(process_ids.unique(), 100)

    # join the download batches
    download_batches = [
        "http://www.boldsystems.org/index.php/API_Public/specimen?ids={}".format(
            "|".join(process_ids_chunk)
        )
        for process_ids_chunk in process_ids
    ]

    return download_batches


# function to parse the xml returned by the BOLD api into a dataframe
# of the required dimensions
def xml_to_dataframe(response, hdf_name_top_100_hits):
    # check if the API is overloaded
    if "You have exceeded" in response.text:
        raise APIOverload

    # parse the response and transform to dataframe if it is valid
    soup = BSoup(response.text, "xml")
    # extract all records
    records = soup.find_all("record")

    # define the tags to look for
    tags = [
        "processid",
        "record_id",
        "bin_uri",
        "sex",
        "lifestage",
        "country",
        "identification_provided_by",
        "identification_method",
        "institution_storing",
    ]

    # extract the relevant data's text fields
    process_id_data = [
        [record.find(tag).text if record.find(tag) != None else "" for tag in tags]
        for record in records
    ]

    # transfer to dataframe to return
    process_id_data = pd.DataFrame(process_id_data, columns=tags)

    # append the data to a new key in the hdf storage
    item_sizes = {
        "processid": 18,
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
            process_id_data,
            format="t",
            data_columns=True,
            min_itemsize=item_sizes,
            complib="blosc:blosclz",
            complevel=9,
        )


# asynchronous request code to send n requests at once
async def as_request(url, as_session, hdf_name_top_100_hits):
    while True:
        try:
            # request the api
            response = await as_session.get(url, timeout=60)
            # wait for 5 seconds to not overload the API
            # time.sleep(5)

            xml_to_dataframe(response, hdf_name_top_100_hits)
            break
        except APIOverload:
            # give user output
            tqdm.write(
                "{}: BOLD API overloaded. Waiting a few minutes.".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )

            time.sleep(300)
            continue
        except (ReadTimeout, ConnectionError):
            tqdm.write(
                "{}: BOLD API did not respond. Retrying.".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )
            continue


# function to limit the maximum concurrent downloads
async def limit_concurrency(url, as_session, semaphore, hdf_name_top_100_hits):
    async with semaphore:
        return await as_request(url, as_session, hdf_name_top_100_hits)


# function to create the asynchronous session
async def as_session(download_links, semaphore, hdf_name_top_100_hits):
    as_session = requests_html.AsyncHTMLSession()
    as_session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.0.4758.82 Safari/537.36"
        }
    )
    retry_strategy = Retry(
        total=15,
        status_forcelist=[400, 401, 403, 404, 413, 429, 502, 503, 504],
        backoff_factor=1,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    as_session.mount("https://", adapter)
    as_session.mount("http://", adapter)

    # generate the tasks to perform
    tasks = (
        limit_concurrency(url, as_session, semaphore, hdf_name_top_100_hits)
        for url in download_links
    )

    # return the result
    return await tqdm_asyncio.gather(*tasks, desc="Downloading additional data")


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
        "http://www.boldsystems.org/index.php/MAS_DataRetrieval_OpenSpecimen?selectedrecordid={}".format(
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
        "{}: Ordering top 100 hits.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # read and sort the hdf file according to the order in the fasta file
    top_100_hits, process_ids = read_and_order(
        fasta_path, hdf_name_top_100_hits, read_fasta
    )

    # filter out all process ids that have already been downloaded
    process_ids_to_download = check_already_downloaded(
        process_ids, hdf_name_top_100_hits
    )

    # give user output
    print(
        "{}: Generating download links for additional data.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # generate download links for the process ids, fetch them with async session
    download_batches = generate_download_links(process_ids_to_download)

    # request the data asynchronously
    # set the concurrent download limit here
    sem = asyncio.Semaphore(1)

    # skip the download if the data is already present
    if not additional_data_present(hdf_name_top_100_hits):
        # start the download
        asyncio.run(
            as_session(
                download_links=download_batches,
                semaphore=sem,
                hdf_name_top_100_hits=hdf_name_top_100_hits,
            )
        )

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
