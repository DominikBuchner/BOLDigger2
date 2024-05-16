import datetime, sys, more_itertools, datetime, math
import pandas as pd
import login
from Bio import SeqIO
from pathlib import Path
from bs4 import BeautifulSoup as BSoup
from tqdm import tqdm
from requests.exceptions import ReadTimeout
from requests.exceptions import ConnectionError


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
def check_already_done(fasta_dict, fasta_name, project_directory):
    # generate a filename for the hdf link storage file
    hdf_name = project_directory.joinpath("{}_download_links.h5.lz".format(fasta_name))

    # try to read the file. If it does not exist continue with the full seq dict.
    # If there are already download links remove those from the fasta dict
    try:
        download_links = pd.read_hdf(hdf_name)
        # remove already saved id from the fasta dict
        for id in download_links["id"]:
            fasta_dict.pop(id, None)

        # return fasta dict and hdf name
        return fasta_dict, hdf_name
    except FileNotFoundError:
        # simply return the fasta dict, since no modifications are needed
        return fasta_dict, hdf_name


# function to gather download links and append them to the hdf storage
def gather_download_links_species_level(session, fasta_dict, hdf_name, query_size):
    # extract query-size elements from the fasta dict
    bold_query = dict(more_itertools.take(query_size, fasta_dict.items()))

    # generate the query string
    bold_query_string = ""

    for key in bold_query.keys():
        bold_query_string += ">{}\n".format(key)
        bold_query_string += "{}\n".format(bold_query[key].seq)

    # generate the data for the post request
    post_request_data = {
        "tabtype": "animalTabPane",
        "historicalDB": "",
        "searchdb": "COX1_SPECIES",
        "sequence": bold_query_string,
    }

    # post the request
    response = session.post(
        "https://boldsystems.org/index.php/IDS_IdentificationRequest",
        data=post_request_data,
        timeout=900,
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

    # set size limits for the columns
    item_sizes = {"id": 100, "url": 130}

    # append results to hdf
    with pd.HDFStore(
        hdf_name, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            "urls_species_level",
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


def main(fasta_path, query_size):
    # log in to BOLD to generate the session
    session = login.bold_login()

    # read the input fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # check if download links have already been generated for any OTUs
    fasta_dict, hdf_name = check_already_done(fasta_dict, fasta_name, project_directory)

    # gather download links at species level until all download links are requested
    # give user output
    print(
        "{}: Starting to gather download links at from species level database.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # calculate stats for the progress bar
    total_requests = math.ceil(len(fasta_dict) / query_size)

    # request the server until all links have been generated
    with tqdm(total=total_requests, desc="Generating download links") as pbar:
        while fasta_dict:
            try:
                fasta_dict = gather_download_links_species_level(
                    session, fasta_dict, hdf_name, query_size
                )
                # update the progress bar
                pbar.update(1)
            except (ReadTimeout, ConnectionError):
                # repeat if there is no response
                # give user output
                tqdm.write(
                    "{}: Starting to gather download links at from species level database.".format(
                        datetime.datetime.now().strftime("%H:%M:%S")
                    )
                )

    # give user output
    print(
        "{}: Species level download links successfully generated.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # continue to download the top 100 hits asynchronosly


main(
    "C:\\Users\\Dominik\\Documents\\GitHub\\BOLDigger2\\test_otus.fasta",
    50,
)
