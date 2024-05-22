import more_itertools
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup as BSoup
from io import StringIO
import requests
from xml.etree import ElementTree as ET


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
    process_ids = top_100_hits["Process_ID"].replace("", np.nan).dropna()

    # return process ids
    return top_100_hits, process_ids


# function to return download links for the BOLD API from the process IDs
def generate_download_links(process_ids):
    # chunk the process IDs into batches of 100 ids
    process_ids = more_itertools.chunked(process_ids, 100)

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

    # return the xml dataframe with None values replaced
    return xml_dataframe.fillna(np.nan)


# main function to run the additional data download
def main(fasta_path, hdf_name_top_100_hits, read_fasta):
    # read and sort the hdf file according to the order in the fasta file
    top_100_hits, process_ids = read_and_order(
        fasta_path, hdf_name_top_100_hits, read_fasta
    )

    # generate download links for the process ids, fetch them with async session
    download_batches = generate_download_links(process_ids)

    # request the data asynchronously


# run only if called as a toplevel script
if __name__ == "__main__":
    main()
