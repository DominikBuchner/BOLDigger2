import pandas as pd
import numpy as np
import datetime
from tqdm import tqdm
from string import punctuation, digits
from joblib import Parallel, delayed
from tqdm_joblib import tqdm_joblib
from pathlib import Path


# funnction to read the sorted top 100 hits including additional data
# remove punctuation and digits from the hits
# only keep the first name of the species column
def read_clean_data(hdf_name_top_100):
    # read the data
    top_100_hits = pd.read_hdf(hdf_name_top_100, key="top_100_hits_additional_data")

    # remove punctuationa and numbers from the taxonomy
    specials = punctuation + digits
    levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    for level in levels:
        top_100_hits[level] = np.where(
            top_100_hits[level].str.contains("[{}]".format(specials)),
            np.nan,
            top_100_hits[level],
        )

    # if there are more than 2 names in the species column only keep the first
    top_100_hits["Species"] = top_100_hits["Species"].str.split(" ").str[0]

    return top_100_hits


# accepts a dataframe for any individual id
# return the threshold to filter for and a taxonomic level
def get_threshold(hit_for_id):
    # find the highest similarity value for the threshold
    threshold = hit_for_id["Similarity"].max()

    # check for no matches and broken records first
    if threshold == 0:
        if hit_for_id["Species"][0] == "NoMatch":
            return 0, "NoMatch"
        elif hit_for_id["Species"][0] == "BrokenRecord":
            return 0, "BrokenRecord"
    else:
        # move through the taxonomy if it is no nomatch hit or broken record
        if threshold >= 97:
            return 97, "Species"
        elif threshold >= 95:
            return 95, "Genus"
        elif threshold >= 90:
            return 90, "Family"
        elif threshold >= 85:
            return 85, "Order"
        elif threshold >= 50:
            return 50, "Class"


## function to move the treshold one level up if no hit is found, also return the new tax level
def move_threshold_up(threshold):
    thresholds = [97, 95, 90, 85, 50]
    levels = ["Species", "Genus", "Family", "Order", "Class"]
    return (
        thresholds[thresholds.index(threshold) + 1],
        levels[thresholds.index(threshold) + 1],
    )


# function to flag the hits
def flag_hits(top_hits, hits_for_id_above_similarity, top_hit):
    # predefine the flags to return
    flags = [False, False, False, False, False]

    # flag 1: Reverse BIN taxonomy
    id_method = top_hits["identification_method"]

    if (
        id_method.str.startswith("BOLD").any()
        or id_method.str.startswith("ID").any()
        or id_method.str.startswith("Tree").any()
    ):
        flags[0] = "1"

    # flag 2: more than one group above the selected threshold
    if len(hits_for_id_above_similarity.index) > 1:
        flags[1] = "2"

    # flag 3: all of the selected top hits are private or early release
    if top_hits["Status"].isin(["Private", "Early-Release"]).all():
        flags[2] = "3"

    # flag 4: top hit is only represented by one record
    if len(top_hits.index) == 1:
        flags[3] = "4"

    # flag 5: top hit is represented by multiple bins
    if len(top_hit["BIN"].str.split(";").item()) > 1:
        flags[4] = "5"

    flags = [i for i in flags if i]
    flags = "_".join(flags)

    return flags


# function to find the top hit for a given ID
def find_top_hit(top_100_hits, idx):
    # only select the respective id
    hits_for_id = (
        top_100_hits.loc[top_100_hits["ID"] == idx].copy().reset_index(drop=True)
    )

    # get the threshold and taxonomic level
    threshold, level = get_threshold(hits_for_id)

    # if NoMatch return the NoMatch, if broken record return BrokenRecord
    if threshold == 0:
        if level == "NoMatch":
            return_value = hits_for_id.query("Species == 'NoMatch'").head(1)
        elif level == "BrokenRecord":
            return_value = hits_for_id.query("Species == 'BrokenRecord'").head(1)

        return_value = return_value[
            [
                "ID",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
                "Similarity",
                "Status",
            ]
        ]
        for value in ["records", "selected_level", "BIN", "flags", "Status"]:
            return_value[value] = np.nan

        return return_value

    # loop through the thresholds until a hit is found
    while True:
        # copy the hits for the respective ID to perform modifications
        hits_for_id_above_similarity = hits_for_id.copy()

        with pd.option_context("future.no_silent_downcasting", True):
            hits_for_id_above_similarity = hits_for_id_above_similarity.replace(
                "", np.nan
            )

        # only select hits above the selected threshold
        hits_for_id_above_similarity = hits_for_id_above_similarity.loc[
            hits_for_id_above_similarity["Similarity"] >= threshold
        ]

        # if no hit remains move up one level until class
        if len(hits_for_id_above_similarity.index) == 0:
            threshold, level = move_threshold_up(threshold)
            continue

        # define the levels for the groupby. care about the selector string later
        all_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
        levels = all_levels[: all_levels.index(level) + 1]

        # only select interesting levels (all levels above and including the selected level)
        hits_for_id_above_similarity = hits_for_id_above_similarity[levels].copy()

        # group the hits by level and then count the appearence
        hits_for_id_above_similarity = pd.DataFrame(
            {
                "count": hits_for_id_above_similarity.groupby(
                    by=levels,
                    sort=False,
                ).size()
            }
        ).reset_index()

        # if the hits still contained np.nan values, groupby will drop them:
        # if theres nothing left after the gruoupby move up one level and continue
        if len(hits_for_id_above_similarity.index) == 0:
            threshold, level = move_threshold_up(threshold)
            continue

        # sort the hits by count
        hits_for_id_above_similarity = hits_for_id_above_similarity.sort_values(
            "count", ascending=False
        )

        # select the hit with the highest count from the dataframe
        # also return the count to display in the top hit table in the end
        top_hits, top_count = (
            hits_for_id_above_similarity.head(1),
            hits_for_id_above_similarity.head(1)["count"].item(),
        )

        # generate the selector string based on the selected level
        query_string = [
            "{} == '{}'".format(level, top_hits[level].item()) for level in levels
        ]
        query_string = " and ".join(query_string)

        # query for the top hits
        top_hits = hits_for_id.query(query_string)

        # collect the bins from the selected top hit
        if threshold == 97:
            top_hit_bins = top_hits["bin_uri"].dropna().unique()
        else:
            top_hit_bins = []

        # select the first match from the top hits table as the top hit
        top_hit = top_hits.head(1).copy()

        # add the record count to the top hit
        top_hit["records"] = top_count

        # add the selected level
        top_hit["selected_level"] = level

        # add the BINs to the top hit
        top_hit["BIN"] = ";".join(top_hit_bins)

        # define level to remove them from low level hits
        levels = ["Class", "Order", "Family", "Genus", "Species"]

        # return species level information if similarity is high enough
        # else remove higher level information form output depending on level
        if threshold == 97:
            break
        else:
            top_hit = top_hit.assign(
                **{k: np.nan for k in levels[levels.index(level) + 1 :]}
            )
            break

    # add flags to the hits
    top_hit["flags"] = flag_hits(top_hits, hits_for_id_above_similarity, top_hit)

    # remove all data that is not needed
    top_hit = top_hit[
        [
            "ID",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
            "Similarity",
            "Status",
            "records",
            "selected_level",
            "BIN",
            "flags",
        ]
    ]

    # return the top hit
    return top_hit


# function to finally save the results
def save_results(
    project_directory,
    fasta_name,
    all_top_hits,
):
    # generate an excel savename
    savename_excel = Path(project_directory).joinpath(
        "{}_identification_result.xlsx".format(fasta_name)
    )
    savename_parquet = Path(project_directory).joinpath(
        "{}_identification_result.parquet.snappy".format(fasta_name)
    )

    # save to excel
    print(
        "{}: Saving results to Excel. This may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    all_top_hits.to_excel(savename_excel, index=False)

    # save to parquet
    print(
        "{}: Saving results to parquet.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    all_top_hits.to_parquet(savename_parquet)


# main function to run the script
def main(hdf_name_top_100, project_directory, fasta_name):
    # give user output
    print(
        "{}: Loading hits to select top hits.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # collect the top 100 hits with additional data
    top_100_hits = read_clean_data(hdf_name_top_100)

    # collect the top hits
    with tqdm_joblib(
        desc="Calculating top hits", total=len(top_100_hits["ID"].unique())
    ) as progress_bar:
        all_top_hits = Parallel(n_jobs=1)(
            delayed(find_top_hit)(top_100_hits, idx)
            for idx in top_100_hits["ID"].unique()
        )

    # concat all the top hits as a final result
    all_top_hits = pd.concat(all_top_hits, axis=0).reset_index(drop=True)

    # save to excel and parquet
    save_results(project_directory, fasta_name, all_top_hits)

    print("{}: Finished.".format(datetime.datetime.now().strftime("%H:%M:%S")))


if __name__ == "__main__":
    main()
