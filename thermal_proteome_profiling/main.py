import argparse
from pathlib import Path

import pandas as pd


# custom type converter for argparse using pathlib
# checks if the given path exists and is a file. Returns a Path object.
def check_filepath(filepath_str: str) -> Path:
    path = Path(filepath_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(
            f"Error: The path '{filepath_str}' does not exist."
        )
    if not path.is_file():
        raise argparse.ArgumentTypeError(
            f"Error: The path '{filepath_str}' is not a file."
        )
    return path


# custom type converter for argparse using pathlib
# checks if the given path exists and is a directory. Returns a Path object.
def check_directory(dir_path_str: str) -> Path:
    path = Path(dir_path_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(
            f"Error: The path '{dir_path_str}' does not exist."
        )
    if not path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Error: The path '{dir_path_str}' is not a directory."
        )
    return path


# Definition of the parser argparse object
parser = argparse.ArgumentParser(
    prog="TPP_Plotter",
    description="A script for the plotting & filtering of thermal proteome profiling melt curves.",
    epilog="E.g. tpp.py -cf home/data/conditions.csv -d home/data/data.csv",
)

parser.add_argument(
    "-d",
    "--data",
    required=True,
    type=check_filepath,
    help="Path to the data file (must exist).",
)

parser.add_argument(
    "-o",
    "--output_folder",
    required=True,
    type=check_directory,
    help="Path to the output folder (must exist).",
)


# read in the data file
def read_data(path: Path) -> pd.DataFrame:
    data_df = pd.read_csv(path)

    filtered_df = data_df[
        data_df["PG.NrOfStrippedSequencesIdentified (Experiment-wide)"] >= 2
    ]

    pass


def main():
    # initiate argparse parser
    args = parser.parse_args()

    # data file
    data = args.data

    # output folder
    out_dir = args.output_folder


if __name__ == "__main__":
    main()
