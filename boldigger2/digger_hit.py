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
    else:
        if hit_for_id["Species"][0] == "NoMatch":
            return 0, "NoMatch"
        elif hit_for_id["Species"][0] == "BrokenRecord":
            return 0, "BrokenRecord"


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
    hits_for_id = top_100_hits.loc[top_100_hits["ID"] == idx].copy()

    # get the threshold and taxonomic level
    threshold, level = get_threshold(hits_for_id)

    # if NoMatch return the NoMatch, if broken record return BrokenRecord
    if threshold == 0:
        if level == "NoMatch":
            return hits_for_id.query("Species == NoMatch").head(1)
        elif level == "BrokenRecord":
            return hits_for_id.query("Species == BrokenRecord").head(1)


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
        # break


if __name__ == "__main__":
    main()
