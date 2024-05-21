import pandas as pd


# function to sort the hdf dataframe according to the order in the fasta file
def read_and_order(fasta_path, hdf_name_top_100_hits, read_fasta):
    # read in the fasta
    fasta_dict, fasta_name, project_directory = read_fasta(fasta_path)
    print(fasta_dict.keys())


# main function to run the additional data download
def main(fasta_path, hdf_name_top_100_hits, read_fasta):
    # read and sort the hdf file according to the order in the fasta file
    read_and_order(fasta_path, hdf_name_top_100_hits, read_fasta)


# run only if called as a toplevel script
if __name__ == "__main__":
    main()
