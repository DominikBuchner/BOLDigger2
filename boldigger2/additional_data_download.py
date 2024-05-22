import more_itertools
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup as BSoup
from io import StringIO
import requests


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


# main function to run the additional data download
def main(fasta_path, hdf_name_top_100_hits, read_fasta):
    # read and sort the hdf file according to the order in the fasta file
    top_100_hits, process_ids = read_and_order(
        fasta_path, hdf_name_top_100_hits, read_fasta
    )

    # generate download links for the process ids, fetch them with async session
    download_batches = generate_download_links(process_ids)

    # request the data
    r = requests.get(download_batches[0])

    print(
        pd.read_xml(
            StringIO(r.text),
        )
    )


# run only if called as a toplevel script
if __name__ == "__main__":
    main()
