import datetime, sys, more_itertools, datetime, requests_html, asyncio, time
import pandas as pd
import numpy as np
from boldigger2 import login, additional_data_download, digger_hit
from Bio import SeqIO
from pathlib import Path
from bs4 import BeautifulSoup as BSoup
from tqdm import tqdm
from tqdm.asyncio import tqdm_asyncio
from requests.exceptions import ReadTimeout
from requests.exceptions import ConnectionError
from io import StringIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from string import punctuation, digits
from boldigger2.exceptions import BadResponseError


# function to read the fasta file into a dictionary
def read_fasta(fasta_path):
    # extract the directory to work in from the fasta path
    fasta_path = Path(fasta_path)
    fasta_name = fasta_path.stem
    project_directory = fasta_path.parent

    # use SeqIO to read the data into a dict - not neccessary to check for fasta type

    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    # trim headers to maximum allowed length of 99 characters, names are preserverd in the SeqRecord object
    fasta_dict = {key[:99]: value for key, value in fasta_dict.items()}

    valid_chars = {
        "A",
        "C",
        "G",
        "T",
        "M",
        "R",
        "W",
        "S",
        "Y",
        "K",
        "V",
        "H",
        "D",
        "B",
        "X",
        "N",
    }

    # check for invalid sequences (invalid characters or sequences that are too short)
    raise_false_fasta = False

    for key in fasta_dict.keys():
        if len(fasta_dict[key].seq) < 80:
            print(
                "{}: Sequence {} is too short (< 80 bp).".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), key
                )
            )
            raise_false_fasta = True
        # check if the sequences contain invalid chars
        elif not set(fasta_dict[key].seq.upper()).issubset(valid_chars):
            print(
                "{}: Sequence {} contains invalid characters.".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), key
                )
            )
            raise_false_fasta = True

    if not raise_false_fasta:
        return fasta_dict, fasta_name, project_directory
    else:
        sys.exit()


def check_already_downloaded(fasta_dict, hdf_name_top_100_hits, database):
    try:
        # load the hdf file
        hdf_top_100_hits = pd.read_hdf(
            hdf_name_top_100_hits, key="top_100_hits_unsorted"
        )

        # filter for the correct database
        hdf_top_100_hits = hdf_top_100_hits.loc[
            hdf_top_100_hits["database"] == database
        ]

        # collect all unique IDs and remove them from the fasta dict
        unique_ids = hdf_top_100_hits["ID"].unique()
        fasta_dict = {
            id: seq for (id, seq) in fasta_dict.items() if id not in unique_ids
        }
    except FileNotFoundError:
        # do nothing if nothing has been downloaded yet
        pass

    # return the (updated) fasta dict
    return fasta_dict


# function to gather download links and append them to the hdf storage
def gather_download_links(session, fasta_dict, query_size, database):
    # extract query-size elements from the fasta dict
    bold_query = dict(more_itertools.take(query_size, fasta_dict.items()))

    # generate the query string
    bold_query_string = ""

    for key in bold_query.keys():
        bold_query_string += ">{}\n".format(key)
        bold_query_string += "{}\n".format(bold_query[key].seq)

    # generate the data for the post request
    if database == "species":
        post_request_data = {
            "tabtype": "animalTabPane",
            "historicalDB": "",
            "searchdb": "COX1_SPECIES",
            "sequence": bold_query_string,
        }
    else:
        post_request_data = {
            "tabtype": "animalTabPane",
            "historicalDB": "",
            "searchdb": "COX1",
            "sequence": bold_query_string,
        }

    # post the request, reduce timeout to 5 minutes, decrease query size instead of just retrying
    response = session.post(
        "https://boldsystems.org/index.php/IDS_IdentificationRequest",
        data=post_request_data,
        timeout=300,
    )

    # extract the download links from the response
    soup = BSoup(response.text, "html5lib")
    download_links = soup.find_all("span", style="text-decoration: none")
    download_links = [
        "http://boldsystems.org" + download_links[i].get("result")
        for i in range(len(download_links))
    ]

    # check if the number of download links matches the query size
    if len(download_links) != len(bold_query):
        # raise a BadResponseError if something happened on BOLDs end
        raise BadResponseError

    # gather the results in a dataframe to easily append them to hdf
    download_dataframe = pd.DataFrame(
        data=zip(bold_query.keys(), download_links), columns=["id", "url"]
    )

    return download_dataframe


# function to update the query size
# accepts that query size and an increase argument. negative increase values will lead to a decrease
def update_query_size(query_size, increase):
    # update the query size via increase
    query_size = query_size + increase

    # return the updated value. can only be in the range of 1 to 50
    if query_size < 1:
        query_size = 1
        return query_size
    elif query_size > 50:
        query_size = 50
        return query_size
    else:
        return query_size


# asynchronous request code to send n requests at once
# database is a string specifying where the data comes from
async def as_request(species_id, url, as_session, database, hdf_name_top_100_hits):
    # add all requests to the eventloop
    # request top 100 hits
    # retry in case of connection error
    while True:
        try:
            response = await as_session.get("{}&display=100".format(url), timeout=60)
            break
        except ConnectionError:
            continue

    # parse the response and pass it to pandas
    response = BSoup(response.text, "html5lib")

    # check for broken records already here in the raw html, since a valid and a broken record both return 4 tables
    broken_record = response.find_all("div", id="kohana_error")

    if len(broken_record) == 1:
        result = pd.DataFrame(
            [[species_id] + ["BrokenRecord"] * 7 + [0.0] + [""] * 2],
            columns=[
                "ID",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
                "Subspecies",
                "Similarity",
                "Status",
                "Process_ID",
            ],
        )

    # read the tables from the response if the result table is not broken
    response_table = pd.read_html(
        StringIO(str(response)),
        header=0,
        converters={"Similarity (%)": float},
        flavor="html5lib",
    )

    # code to generate the no match table
    if len(response_table) == 2:
        result = pd.DataFrame(
            [[species_id] + ["NoMatch"] * 7 + [0.0] + [""] * 2],
            columns=[
                "ID",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
                "Subspecies",
                "Similarity",
                "Status",
                "Process_ID",
            ],
        )
    else:
        # further strip down the html to only collect the top 100 hits irrespective of database used
        if database == "species":
            response = response.select('h3:-soup-contains("Top 100 Matches")')[
                -1
            ].find_next("table", class_="table resultsTable noborder")
        else:
            response = response.select('h3:-soup-contains("Top 100 Matches")')[
                -1
            ].find_next("table", class_="resultsTable noborder")

        # finally scrape the correct response table
        result = pd.read_html(
            StringIO(str(response)),
            header=0,
            converters={"Similarity (%)": float},
            flavor="html5lib",
        )[-1]

        ids = [
            tag.get("id") for tag in response.find_all(class_="publicrecord")
        ]  # collect process ids for public records
        result.columns = [
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
            "Subspecies",
            "Similarity",
            "Status",
        ]
        result["Process_ID"] = [
            ids.pop(0) if status else np.nan
            for status in np.where(result["Status"] == "Published", True, False)
        ]

        # add an identifier column to be able to sort the table
        result.insert(0, "ID", species_id)

    # fill na values with empty strings to make frames compatible with hdf format
    result = result.fillna("")

    # add the database and a timestamp to the result table
    result["database"] = database
    result["request_date"] = pd.Timestamp.now().strftime("%Y-%m-%d %X")

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

    # append results to hdf
    with pd.HDFStore(
        hdf_name_top_100_hits, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            "top_100_hits_unsorted",
            result,
            format="t",
            data_columns=True,
            min_itemsize=item_sizes,
            complib="blosc:blosclz",
            complevel=9,
        )

    if database == "species":
        # give user output
        tqdm.write(
            "{}: Downloaded top 100 species level records for {}".format(
                datetime.datetime.now().strftime("%H:%M:%S"), species_id
            )
        )

    else:
        tqdm.write(
            "{}: Downloaded top 100 hits of all records for {}".format(
                datetime.datetime.now().strftime("%H:%M:%S"), species_id
            )
        )


# function to limit the maximum concurrent downloads
async def limit_concurrency(
    species_id, url, as_session, database, hdf_name_top_100_hits, semaphore
):
    async with semaphore:
        return await as_request(
            species_id, url, as_session, database, hdf_name_top_100_hits
        )


# function to create the asynchronous session
async def as_session(
    download_links_species, database, hdf_name_top_100_hits, semaphore
):
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

    # create all requests
    tasks = download_links_species.copy()
    tasks = (
        limit_concurrency(
            id, url, as_session, database, hdf_name_top_100_hits, semaphore
        )
        for id, url in zip(tasks["id"], tasks["url"])
    )

    # return the result
    return await asyncio.gather(*tasks)


def check_valid_species_records(fasta_dict, hdf_name_top_100_hits, thresholds):
    # read the hdf
    top_100_hits_species = pd.read_hdf(
        hdf_name_top_100_hits, key="top_100_hits_unsorted"
    )
    # filter for species level database only
    top_100_hits_species = top_100_hits_species.loc[
        top_100_hits_species["database"] == "species"
    ]
    # fill empty string with NA values for easy filtering
    with pd.option_context("future.no_silent_downcasting", True):
        top_100_hits_species = top_100_hits_species.replace("", np.nan)

    # remove digits and punctuation in species column
    specials = punctuation + digits
    top_100_hits_species["Species"] = np.where(
        top_100_hits_species["Species"].str.contains("[{}]".format(specials)),
        np.nan,
        top_100_hits_species["Species"],
    )
    # drop na values from species column
    top_100_hits_species = top_100_hits_species.dropna(subset="Species")
    # only keep values with similarity >= 97%
    top_100_hits_species = top_100_hits_species.loc[
        top_100_hits_species["Similarity"] >= thresholds[0]
    ]

    # pop those values from the fasta dict
    fasta_dict = {
        key: value
        for key, value in fasta_dict.items()
        if key not in top_100_hits_species["ID"].unique()
    }

    return fasta_dict


def main(fasta_path, username="", password="", thresholds=[]):
    # log in to BOLD to generate the session, initialize the query size
    session, username, password = login.bold_login(username=username, password=password)
    query_size = 1

    # read the input fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # generate a name for the top hits hdf file
    hdf_name_top_100_hits = project_directory.joinpath(
        "{}_top_100_hits.h5.lz".format(fasta_name)
    )

    # check if any of the ids have been downloaded and saved already. If so remove them from the fasta dict
    fasta_dict = check_already_downloaded(fasta_dict, hdf_name_top_100_hits, "species")

    # start the download for the species level database
    # give user output
    print(
        "{}: Starting to download from the species level database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # request the server until all links have been generated
    if fasta_dict:
        with tqdm(total=len(fasta_dict), desc="Generating download links") as pbar:
            # generate download links first
            while fasta_dict:
                try:
                    # gather the returned download links to download them straight away
                    download_dataframe = gather_download_links(
                        session, fasta_dict, query_size, database="species"
                    )

                    # set the semaphore to the query size of the original query
                    sem = asyncio.Semaphore(query_size)
                    pbar_update = query_size
                    # update the query size by 5
                    query_size = update_query_size(query_size, 5)

                    # give user output
                    if query_size != 50:
                        tqdm.write(
                            "{}: Query size updated to {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )
                    else:
                        tqdm.write(
                            "{}: Query size kept at {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )

                    # download the data from the generated links, update the fasta dict after successfull download
                    # update the description of the progress bar
                    pbar.set_description("Downloading data")

                    # catch sometimes malformed urls here, produces duplicates in the top 100 download, will be removed
                    # when downloading additional data / ordering the hits
                    try:
                        # run the control loop
                        asyncio.run(
                            as_session(
                                download_dataframe,
                                database="species",
                                hdf_name_top_100_hits=hdf_name_top_100_hits,
                                semaphore=sem,
                            )
                        )
                    except IndexError:
                        continue

                    # update the progress bar
                    pbar.update(pbar_update)
                    pbar.set_description("Generating download links")

                    # update the fasta dict
                    fasta_dict = {
                        id: seq
                        for (id, seq) in fasta_dict.items()
                        if id not in set(download_dataframe["id"])
                    }

                except (ReadTimeout, ConnectionError):
                    # repeat if there is no response
                    # update the query size
                    query_size = update_query_size(query_size, -5)
                    # give user output
                    if query_size != 1:
                        tqdm.write(
                            "{}: BOLD did not respond. Retrying with reduced query size of {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )
                    else:
                        tqdm.write(
                            "{}: BOLD did not respond. Keeping query size at {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )
                except BadResponseError:
                    tqdm.write(
                        "{}: BOLD did not return a sufficient number of download links. Retrying".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )
                    # wait for 3 minutes to give the BOLD Server a break
                    time.sleep(180)

                    # login again
                    session, username, password = login.bold_login(
                        username=username, password=password
                    )

    # give user output
    print(
        "{}: Starting to download from the all records database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # give user output
    print(
        "{}: Performing second login for requesting links from the all records database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # reread the fasta to generate a fresh fasta dict
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # filter the fasta dict for hits no having a species level hit, perform a second log in
    session, username, password = login.bold_login(username=username, password=password)
    fasta_dict = check_valid_species_records(
        fasta_dict, hdf_name_top_100_hits, thresholds=thresholds
    )

    # gather download links at all barcode records level until all download links are requested
    # give user output
    print(
        "{}: Starting to gather download links from the all records database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # check if any of the ids have been downloaded and saved already. If so remove them from the fasta dict
    fasta_dict = check_already_downloaded(
        fasta_dict, hdf_name_top_100_hits, "all_records"
    )

    # request the server until all links have been generated
    if fasta_dict:
        with tqdm(total=len(fasta_dict), desc="Generating download links") as pbar:
            # generate download links first
            while fasta_dict:
                try:
                    # gather the returned download links to download them straight away
                    download_dataframe = gather_download_links(
                        session, fasta_dict, query_size, database="all_records"
                    )

                    # set the semaphore to the query size of the original query
                    sem = asyncio.Semaphore(query_size)
                    pbar_update = query_size
                    # update the query size by 5
                    query_size = update_query_size(query_size, 5)

                    # give user output
                    if query_size != 50:
                        tqdm.write(
                            "{}: Query size updated to {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )
                    else:
                        tqdm.write(
                            "{}: Query size kept at {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )

                    # download the data from the generated links, update the fasta dict after successfull download
                    # update the description of the progress bar
                    pbar.set_description("Downloading data")

                    try:
                        # run the control loop
                        asyncio.run(
                            as_session(
                                download_dataframe,
                                database="all_records",
                                hdf_name_top_100_hits=hdf_name_top_100_hits,
                                semaphore=sem,
                            )
                        )
                    except IndexError:
                        continue

                    # update the progress bar
                    pbar.update(pbar_update)
                    pbar.set_description("Generating download links")

                    # update the fasta dict
                    fasta_dict = {
                        id: seq
                        for (id, seq) in fasta_dict.items()
                        if id not in set(download_dataframe["id"])
                    }

                except (ReadTimeout, ConnectionError):
                    # repeat if there is no response
                    # update the query size
                    query_size = update_query_size(query_size, -5)
                    # give user output
                    if query_size != 1:
                        tqdm.write(
                            "{}: BOLD did not respond. Retrying with reduced query size of {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )
                    else:
                        tqdm.write(
                            "{}: BOLD did not respond. Keeping query size at {}.".format(
                                datetime.datetime.now().strftime("%H:%M:%S"), query_size
                            )
                        )
                except BadResponseError:
                    tqdm.write(
                        "{}: BOLD did not return a sufficient number of download links. Retrying".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )
                    # wait for 3 minutes to give the BOLD Server a break
                    time.sleep(180)

                    # login again
                    session, username, password = login.bold_login(
                        username=username, password=password
                    )

    # give user output
    print(
        "{}: All records top 100 records successfully downloaded.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # download the additional data if it is not present yet
    additional_data_download.main(fasta_path, hdf_name_top_100_hits, read_fasta)

    # filter for the top hits
    digger_hit.main(
        hdf_name_top_100_hits, project_directory, fasta_name, thresholds=thresholds
    )


# run only if called as a toplevel script
if __name__ == "__main__":
    main()
