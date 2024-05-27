import pandas as pd
import numpy as np
import datetime
from tqdm import tqdm
from string import punctuation, digits


# funnction to read the sorted top 100 hits including additional data
# remove punctuation and digits from the hits
# only keep the first name of the species column
def read_clean_data(hdf_name_top_100):
    # read the data
    # top_100_hits = pd.read_hdf(hdf_name_top_100, key="top_100_hits_additional_data")
    top_100_hits = pd.read_excel(
        "C:\\Users\\Dominik\\Documents\\GitHub\\BOLDigger2\\test_otus_top_100_hits_part_0.xlsx"
    )

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
            return 98, "Species"
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
    thresholds = [98, 95, 90, 85, 50]
    levels = ["Species", "Genus", "Family", "Order", "Class"]
    return (
        thresholds[thresholds.index(threshold) + 1],
        levels[thresholds.index(threshold) + 1],
    )


# function to find the top hit for a given ID
def find_top_hit(top_100_hits, idx):
    # only select the respective id
    hits_for_id = (
        top_100_hits.loc[top_100_hits["ID"] == idx].copy().reset_index(drop=True)
    )

    # use a placeholder in all empty cells for pd.query to work properly
    hits_for_id = hits_for_id.replace(np.nan, "placeholder")

    # get the threshold and taxonomic level
    threshold, level = get_threshold(hits_for_id)

    # if NoMatch return the NoMatch, if broken record return BrokenRecord
    if threshold == 0:
        if level == "NoMatch":
            return hits_for_id.query("Species == 'NoMatch'").head(1)
        elif level == "BrokenRecord":
            return hits_for_id.query("Species == 'BrokenRecord'").head(1)

    # loop through the thresholds until a hit is found
    while True:
        # copy the hits for the respective ID to perform modifications
        hits_for_id_above_similarity = hits_for_id.copy()
        # only select hits above the selected threshold
        hits_for_id_above_similarity = hits_for_id_above_similarity.loc[
            hits_for_id_above_similarity["Similarity"] >= threshold
        ]

        # group the hits by level and then count the appearence
        hits_for_id_above_similarity = pd.DataFrame(
            {
                "count": hits_for_id_above_similarity.groupby(
                    ["Phylum", "Class", "Order", "Family", "Genus", "Species"],
                    sort=False,
                ).size()
            }
        ).reset_index()

        # sort the hits by count
        hits_for_id_above_similarity = hits_for_id_above_similarity.sort_values(
            "count", ascending=False
        )

        # drop na values at the respective level, replace the placeholder first, then put the placeholder back in place
        with pd.option_context("future.no_silent_downcasting", True):
            hits_for_id_above_similarity = (
                hits_for_id_above_similarity.replace("placeholder", np.nan)
                .dropna(subset=level)
                .fillna("placeholder")
            )

        # if no hit remains move up one level until class
        if len(hits_for_id_above_similarity.index) == 0:
            threshold, level = move_threshold_up(threshold)
            continue
        else:
            # select the hit with the highest count from the dataframe
            # also return the count to display in the top hit table in the end
            top_hit, top_count = (
                hits_for_id_above_similarity.head(1),
                hits_for_id_above_similarity.head(1)["count"],
            )
            top_hit = hits_for_id.query(
                "Class == '{}' and Order == '{}' and Family == '{}' and Genus == '{}' and Species == '{}'".format(
                    top_hit["Class"].item(),
                    top_hit["Order"].item(),
                    top_hit["Family"].item(),
                    top_hit["Genus"].item(),
                    top_hit["Species"].item(),
                )
            )

            # return the BIN if it is unique for all selected hits
            print(top_hit["bin_uri"].unique(), top_count)


# main function to run the script
def main(hdf_name_top_100):
    # give user output
    print(
        "{}: Loading hits to select top hits.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # collect the top 100 hits with additional data
    top_100_hits = read_clean_data(hdf_name_top_100)

    # PARALLIZE THIS IN THE END
    for idx in top_100_hits["ID"].unique():
        print(find_top_hit(top_100_hits, idx))


if __name__ == "__main__":
    main()
