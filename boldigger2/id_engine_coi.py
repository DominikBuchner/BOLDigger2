import datetime, sys, more_itertools, datetime
import pandas as pd
import login
from Bio import SeqIO
from pathlib import Path
from bs4 import BeautifulSoup as BSoup


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
        ## CONTINUE CODE HERE IF SOMETHING EXISTS ALREADY
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

    # give user output
    print(
        "{}: Gathering download links. This will take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

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

    print(download_links)


def main(fasta_path, query_size):
    # log in to BOLD to generate the session
    session = login.bold_login()

    # read the input fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)

    # check if download links have already been generated for any OTUs
    fasta_dict, hdf_name = check_already_done(fasta_dict, fasta_name, project_directory)

    # gather download links at species level
    gather_download_links_species_level(session, fasta_dict, hdf_name, query_size)


main("C:\\Users\\Dominik\\Documents\\GitHub\\BOLDigger2\\test_otus.fasta", 1)
