import datetime, sys
from Bio import SeqIO
from pathlib import Path


# function to read the fasta file into a dictionary
def read_fasta(fasta_path):
    # use SeqIO to read the data into a dict - not neccessary to check for fasta type
    fasta_dict = SeqIO.to_dict(SeqIO.parse(Path(fasta_path), "fasta"))

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

    return fasta_dict


read_fasta("C:\\Users\\Dominik\\Documents\\GitHub\\BOLDigger2\\test_otus.fasta")
