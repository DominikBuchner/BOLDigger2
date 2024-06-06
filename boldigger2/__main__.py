import argparse, sys, datetime
from boldigger2 import id_engine_coi


# main function for the command line interface
def main():
    # initialize the parse and display default behavior if called without arguments
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=35)
    # define the parser
    parser = argparse.ArgumentParser(
        prog="boldigger2",
        description="a Python package to identify and organise sequences with the Barcode of Life Data systems",
    )

    # display help when no argument is called
    parser.set_defaults(func=lambda x: parser.print_help())

    # add the subparsers
    subparsers = parser.add_subparsers(dest="function")

    # add the identify parser
    parser_identify = subparsers.add_parser(
        "identify", help="Run the COI identification engine."
    )

    # add the only argument (fasta path)
    parser_identify.add_argument(
        "fasta_file",
        help="Path to the fasta file or fasta file in current working directory.",
    )

    # add the optional argument username
    parser_identify.add_argument(
        "-username",
        default="",
        help="BOLD username",
    )

    # add the optional argument password
    parser_identify.add_argument(
        "-password",
        default="",
        help="BOLD password",
    )

    # add the optional argument thresholds
    parser_identify.add_argument(
        "-thresholds",
        nargs="+",
        type=int,
        help="BOLD password",
    )

    # add version control NEEDS TO BE UPDATED
    parser.add_argument("--version", action="version", version="1.0.6")

    # parse the arguments
    arguments = parser.parse_args()

    # only use the threshold provided by the user replace the rest with defaults
    default_thresholds = [97, 95, 90, 85, 50]
    thresholds = []

    for i in range(5):
        try:
            thresholds.append(arguments.thresholds[i])
        except (IndexError, TypeError):
            thresholds.append(default_thresholds[i])

    if arguments.thresholds:
        # give user output
        print(
            "{}: Default thresholds changed!\n{}: Species: {}, Genus: {}, Family: {}, Genus: {}, Species: {}".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                datetime.datetime.now().strftime("%H:%M:%S"),
                *thresholds
            )
        )

    # print help if no argument is provided
    if len(sys.argv) == 1:
        arguments.func(arguments)

    # run the identification engine
    if arguments.function == "identify":
        id_engine_coi.main(
            arguments.fasta_file,
            username=arguments.username,
            password=arguments.password,
            thresholds=thresholds,
        )


# run only if called as a top level script
if __name__ == "__main__":
    main()
