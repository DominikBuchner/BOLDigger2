import more_itertools, asyncio, requests_html, datetime, time
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm.asyncio import tqdm_asyncio
from pathlib import Path
from tqdm import tqdm


# function to sort the hdf dataframe according to the order in the fasta file
def read_and_order(fasta_path, hdf_name_top_100_hits, read_fasta):
    # read in the fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # the keys are in perfect order, use those to order the hdf file
    sorter = {name: idx for idx, name in enumerate(fasta_dict.keys())}

    # read the hdf data that needs to be sorted
    top_100_hits = pd.read_hdf(hdf_name_top_100_hits, key="top_100_hits_unsorted")

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
def xml_to_dataframe(xml_string):
    # extract the root of the xml tree
    root = ET.fromstring(xml_string)
    records = root.findall("record")

    # store data here, append data on the go
    xml_dataframe = pd.DataFrame()

    # add process id, record id, bin uri
    for column_name in ["processid", "record_id", "bin_uri"]:
        xml_dataframe[column_name] = [
            (
                record.find(column_name)
                if record.find(column_name) == None
                else record.find(column_name).text
            )
            for record in records
        ]

    # gather the remaining data from the API
    # collect subheaders in a dict, and subelements as lists of the respective key
    remaining_elements = {
        "specimen_identifiers": ["institution_storing"],
        "specimen_desc": ["sex", "lifestage"],
        "collection_event": ["country"],
        "taxonomy": ["identification_provided_by", "identification_method"],
    }

    for subheader in remaining_elements.keys():
        # find the respective subheader in the xml
        subheader_data = [record.find(subheader) for record in records]

        # collect the respective subelements and push them into the xml dataframe
        for column_name in remaining_elements[subheader]:
            # extract the xml elements first
            elements = [
                record.find(column_name) if record != None else None
                for record in subheader_data
            ]
            # dual loop to not loose the None values,to be able to append to dataframe
            elements = [
                element.text if element != None else None for element in elements
            ]

            xml_dataframe[column_name] = elements

    # return the xml dataframe
    return xml_dataframe


# asynchronous request code to send n requests at once
async def as_request(url, as_session):
    while True:
        try:
            # request the api
            response = await as_session.get(url, timeout=60)
            # wait for 3 seconds before sending next request, otherwise BOLD API will block additional requests
            time.sleep(3)

            xml_dataframe = xml_to_dataframe(response.text)
            break
        except ET.ParseError:
            # give user output
            tqdm.write(
                "{}: BOLD API overloaded. Waiting a few minutes.".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )

            time.sleep(600)
            continue

    return xml_dataframe


# function to limit the maximum concurrent downloads
async def limit_concurrency(url, as_session, semaphore):
    async with semaphore:
        return await as_request(url, as_session)


# function to create the asynchronous session
async def as_session(download_links, semaphore):
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
    tasks = (limit_concurrency(url, as_session, semaphore) for url in download_links)

    # return the result
    return await tqdm_asyncio.gather(*tasks, desc="Downloading additional data")


# function to add the additional data to the top 100 hits
def add_additional_data(
    hdf_name_top_100_hits, top_100_hits, process_ids, additional_data
):
    # concat the additional data downloaded
    additional_data = pd.concat(additional_data, axis=0).reset_index(drop=True)

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

    # give user output
    print(
        "{}: Generating download links for additional data.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # generate download links for the process ids, fetch them with async session
    download_batches = generate_download_links(process_ids)

    # request the data asynchronously
    # set the concurrent download limit here
    sem = asyncio.Semaphore(1)

    # skip the download if the data is already present
    if not additional_data_present(hdf_name_top_100_hits):
        # start the download
        additional_data = asyncio.run(
            as_session(download_links=download_batches, semaphore=sem)
        )

        # add the metadata to the top 100 hits, push to a new hdf table
        top_100_hits = add_additional_data(
            hdf_name_top_100_hits, top_100_hits, process_ids, additional_data
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
