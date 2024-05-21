import datetime, sys, more_itertools, datetime, requests_html, asyncio, login
import pandas as pd
import numpy as np
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
    for key in fasta_dict.keys():
        if len(fasta_dict[key].seq) < 80:
            print(
                "{}: Sequence {} is too short (< 80 bp).".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), key
                )
            )
            sys.exit()
        # check if the sequences contain invalid chars
        elif not set(fasta_dict[key].seq.upper()).issubset(valid_chars):
            print(
                "{}: Sequence {} contains invalid characters.".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), key
                )
            )
            sys.exit()

    return fasta_dict, fasta_name, project_directory


# function to check for which OTUs download links have already been generated
def check_already_done(fasta_dict, fasta_name, project_directory, database):
    # generate a filename for the hdf link storage file
    hdf_name = project_directory.joinpath("{}_download_links.h5.lz".format(fasta_name))

    # try to read the file. If it does not exist continue with the full seq dict.
    # If there are already download links remove those from the fasta dict
    try:
        download_links = pd.read_hdf(hdf_name)
        download_links = download_links.loc[download_links["database"] == database]
        # remove already saved id from the fasta dict
        for id in download_links["id"]:
            fasta_dict.pop(id, None)

        # return fasta dict and hdf name
        return fasta_dict, hdf_name
    except FileNotFoundError:
        # simply return the fasta dict, since no modifications are needed
        return fasta_dict, hdf_name


# function to gather download links and append them to the hdf storage
def gather_download_links(session, fasta_dict, hdf_name, query_size, database):
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

    # gather the results in a dataframe to easily append them to hdf
    download_dataframe = pd.DataFrame(
        data=zip(bold_query.keys(), download_links), columns=["id", "url"]
    )

    # add the database column
    download_dataframe["database"] = database

    # set size limits for the columns
    item_sizes = {"id": 100, "url": 130, "database": 15}

    # append results to hdf
    with pd.HDFStore(
        hdf_name, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            "download_links",
            download_dataframe,
            format="t",
            data_columns=True,
            min_itemsize=item_sizes,
            complib="blosc:blosclz",
            complevel=9,
        )

    # remove finished ids from the fasta dict
    for id in bold_query.keys():
        fasta_dict.pop(id, None)

    # return the fasta dict again to continue
    return fasta_dict


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
    response = await as_session.get("{}&display=100".format(url), timeout=60)
    # parse the response and pass it to pandas
    response = BSoup(response.text, "html5lib")

    # check for broken records already here in the raw html, since a valid and a broken record both return 4 tables
    broken_record = response.find_all("div", id="kohana_error")

    if len(broken_record) == 1:
        result = pd.DataFrame(
            [[species_id] + ["Broken record"] * 7 + [0.0] + [""] * 2],
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
            [[species_id] + ["No Match"] * 7 + [0.0] + [""] * 2],
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

    # add the database to the result table
    result["database"] = database

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
    }

    # append results to hdf
    with pd.HDFStore(
        hdf_name_top_100_hits, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            "top_100_hits",
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
    return await tqdm_asyncio.gather(*tasks, desc="Downloading top 100 hits")


# function to check of some of the top 100 hits have already been downloaded
def top_100_downloaded(hdf_name_top_100_hits, hdf_name_download_links, database):
    # load the download links
    download_links = pd.read_hdf(hdf_name_download_links)
    # filter for the correct database
    download_links = download_links.loc[download_links["database"] == database]
    try:
        # read already downloaded top 100 hits
        top_100_hits = pd.read_hdf(hdf_name_top_100_hits)
        # only select the correct database
        top_100_hits = top_100_hits.loc[top_100_hits["database"] == database]
        # find all ids that already have been selected
        ids_done = top_100_hits["ID"].unique()

        # remove all links that are already done
        download_links = download_links[~download_links["id"].isin(ids_done)]
    except FileNotFoundError:
        # if there is no data downloaded yet just return all download links
        pass

    return download_links


# function to check which species level hits can used straight away (similarity > 97%, valid name)
# and where additional data needs to be requested, returns a fresh fasta dict with everything that needs to
# be identified again
def check_valid_species_records(fasta_dict, hdf_name_top_100_hits):
    # read the hdf
    top_100_hits_species = pd.read_hdf(hdf_name_top_100_hits)
    # filter for species level database only
    top_100_hits_species = top_100_hits_species.loc[
        top_100_hits_species["database"] == "species"
    ]
    # fill empty string with NA values for easy filtering
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
        top_100_hits_species["Similarity"] >= 97
    ]

    # pop those values from the fasta dict
    for key in top_100_hits_species["ID"]:
        fasta_dict.pop(key, None)

    return fasta_dict


def main(fasta_path):
    # log in to BOLD to generate the session, initialize the query size
    session = login.bold_login()
    query_size = 1

    # read the input fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # check if download links have already been generated for any OTUs
    fasta_dict, hdf_name_download_links = check_already_done(
        fasta_dict, fasta_name, project_directory, database="species"
    )

    # gather download links at species level until all download links are requested
    # give user output
    print(
        "{}: Starting to gather download links at from species level database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # request the server until all links have been generated
    if fasta_dict:
        with tqdm(total=len(fasta_dict), desc="Generating download links") as pbar:
            while fasta_dict:
                try:
                    fasta_dict = gather_download_links(
                        session,
                        fasta_dict,
                        hdf_name_download_links,
                        query_size,
                        database="species",
                    )
                    # update the progress bar
                    pbar.update(query_size)
                    # update the query size by 5
                    query_size = update_query_size(query_size, 5)
                except (ReadTimeout, ConnectionError):
                    # repeat if there is no response
                    # give user output
                    tqdm.write(
                        "{}: BOLD did not respond. Retrying with reduced query size.".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )

                    # update the query size
                    query_size = update_query_size(query_size, -5)

    # give user output
    print(
        "{}: Species level download links successfully generated.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # generate a name for the top hits hdf file
    hdf_name_top_100_hits = project_directory.joinpath(
        "{}_top_100_hits.h5.lz".format(fasta_name)
    )

    # check if some of the links have already been downloaded
    download_links_species = top_100_downloaded(
        hdf_name_top_100_hits, hdf_name_download_links, database="species"
    )

    # set the concurrent download limit here
    # code has been timed to find the optimal download time
    sem = asyncio.Semaphore(6)

    # continue to download the top 100 hits asynchronosly for species level
    if not len(download_links_species.index) == 0:
        asyncio.run(
            as_session(
                download_links_species, "species", hdf_name_top_100_hits, semaphore=sem
            )
        )

    # give user output
    print(
        "{}: Species level top 100 records successfully downloaded.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # reread the fasta to generate a fresh fasta dict
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # filter the fasta dict for hits no having a species level hit, reset query size to 1
    fasta_dict = check_valid_species_records(fasta_dict, hdf_name_top_100_hits)
    query_size = 1
    # gather download links at all barcode records level until all download links are requested
    # give user output
    print(
        "{}: Starting to gather download links from the all records database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # check if download links have already been generated for any OTUs
    fasta_dict, hdf_name_download_links = check_already_done(
        fasta_dict, fasta_name, project_directory, database="all_records"
    )

    # request the server until all links have been generated
    if fasta_dict:
        with tqdm(total=len(fasta_dict), desc="Generating download links") as pbar:
            while fasta_dict:
                try:
                    fasta_dict = gather_download_links(
                        session,
                        fasta_dict,
                        hdf_name_download_links,
                        query_size,
                        database="all_records",
                    )
                    # update the progress bar
                    pbar.update(query_size)
                    # update the query size by 5
                    query_size = update_query_size(query_size, 5)
                except (ReadTimeout, ConnectionError):
                    # repeat if there is no response
                    # give user output
                    tqdm.write(
                        "{}: BOLD did not respond. Retrying with reduced query size.".format(
                            datetime.datetime.now().strftime("%H:%M:%S")
                        )
                    )

                    # update the query size
                    query_size = update_query_size(query_size, -5)

    # give user output
    print(
        "{}: All records download links successfully generated.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # check if some of the links have already been downloaded
    download_links_all_records = top_100_downloaded(
        hdf_name_top_100_hits, hdf_name_download_links, database="all_records"
    )

    # set the concurrent download limit here
    # code has been timed to find the optimal download time
    sem = asyncio.Semaphore(6)

    # continue to download the top 100 hits asynchronosly for species level
    if not len(download_links_all_records.index) == 0:
        asyncio.run(
            as_session(
                download_links_all_records,
                "all_records",
                hdf_name_top_100_hits,
                semaphore=sem,
            )
        )

    # give user output
    print(
        "{}: All records top 100 records successfully downloaded.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )


main(
    "C:\\Users\\Dominik\\Documents\\GitHub\\BOLDigger2\\test_otus.fasta",
)
