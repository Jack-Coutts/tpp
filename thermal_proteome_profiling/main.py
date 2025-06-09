import argparse
import logging
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Configure logging to stdout with basic settings
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)


# custom type converter for argparse using pathlib
# checks if the given path exists and is a file. Returns a Path object.
def check_filepath(filepath_str: str) -> Path:
    path = Path(filepath_str)
    3
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

# Add this after your existing parser arguments, before the functions
mode_group = parser.add_mutually_exclusive_group()

mode_group.add_argument("-g", "--gene", type=str, help="Plot a single gene by name")

mode_group.add_argument(
    "-gl",
    "--gene-list",
    type=check_filepath,
    help="Path to a text file containing gene names (one per line)",
)

mode_group.add_argument(
    "-f",
    "--filter",
    action="store_true",
    help="Filter proteins and add to a gene_list file.",
)


def read_gene_list(file_path: Path) -> list:
    """Read gene names from a text file, one per line"""
    try:
        with open(file_path, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        logging.info(f"{len(genes)} genes read from gene list")
        return genes
    except Exception as e:
        logging.error(f"Error reading {Path}:")
        logging.error(f"{e}")


# read in the data file
def read_data(path: Path) -> Tuple[pd.DataFrame, float, float]:
    try:
        data_df = pd.read_csv(path, delimiter="\t")
        filtered_df = data_df[data_df["# Unique Total Peptides"] >= 2]
        highest_fld_chng = filtered_df["AVG Log2 Ratio"].max()
        lowest_fld_chng = filtered_df["AVG Log2 Ratio"].min()
        logging.info(f"{Path} read success")
        return filtered_df, highest_fld_chng, lowest_fld_chng
    except Exception as e:
        logging.error(f"Error reading {Path}:")
        logging.error(f"{e}")


def plotter(df, y_min, y_max, cell_line, gene_name, output_dir, variants):
    try:
        # Set seaborn style
        sns.set_style("whitegrid")
        plt.rcParams["font.size"] = 11

        # ensure line colour consistency
        sorted_variants = sorted(variants, key=lambda x: x[0])
        variant_palette = {sorted_variants[0]: "#1f77b4", sorted_variants[1]: "#d62728"}

        sns.lineplot(
            data=df,
            x="concentration",  # categorical, equally spaced
            y="log_fld_change",
            hue="variant",  # separate lines per variant
            palette=variant_palette,
            hue_order=[sorted_variants[0], sorted_variants[1]],
            marker="o",
            linewidth=2,
            markersize=6,
        )

        plt.ylim(y_min, y_max)
        plt.xlabel("Concentration", fontsize=12)
        plt.ylabel("AVG Log2 Ratio", fontsize=12)
        plt.title(f"{cell_line} {gene_name}", fontsize=14, fontweight="bold", pad=20)
        # Add horizontal line at y=0
        plt.axhline(y=0, color="gray", linewidth=0.8, alpha=0.7)
        # Legend styling
        plt.legend(frameon=True, fancybox=True, shadow=True)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{gene_name}_{cell_line}.png")
        plt.close()
        logging.info(f"{gene_name} from {cell_line} plotted")

    except Exception as e:
        logging.error(f"Error plotting {gene_name} from {cell_line}:")
        logging.error(f"{e}")


def filter_or_plot_gene(
    gene_name,
    data_df,
    output_dir,
    highest_y_value,
    lowest_y_value,
    filter=False,
    filter_flat=False,
):
    # Filter data for the specific gene
    gene_data = data_df[data_df["Genes"] == gene_name]

    # check gene data exists
    if gene_data.empty:
        print(f"No data found for gene: {gene_name}")
        return None

    # cell line name
    cell_line = gene_data["Comparison (group1/group2)"].iloc[0].split("_")[0]

    # protein variants
    variants = gene_data["Comparison (group1/group2)"].str.split("_").str[1].unique()

    # concentrations
    concentrations = (
        gene_data["Comparison (group1/group2)"]
        .str.split("_")
        .str[2]  # get the third element of each split
        .str.split(" ")
        .str[0]  # get the part before the first space
        .unique()  # get unique values
    )

    sorted_concentrations = sorted(concentrations, key=lambda x: float(x.rstrip("uM")))

    # create dictionary containing each line to be plotted
    lines = {}

    # create the line for each variant
    for variant in variants:
        lines[variant] = []
        var_mask = (
            gene_data["Comparison (group1/group2)"].str.split("_").str[1] == variant
        )
        variant_df = gene_data[var_mask]

        for concentration in sorted_concentrations:
            log_fld_change = variant_df.loc[
                variant_df["Comparison (group1/group2)"].str.contains(
                    concentration, na=False
                ),
                "AVG Log2 Ratio",
            ].iloc[0]
            lines[variant].append(log_fld_change)

    plot_df = pd.DataFrame(lines, index=sorted_concentrations).reset_index()
    plot_df = plot_df.melt(
        id_vars="index", var_name="variant", value_name="log_fld_change"
    )
    plot_df.rename(columns={"index": "concentration"}, inplace=True)

    # filter
    if filter:
        if filter_proteins(lines, sorted_concentrations, variants, flat=filter_flat):
            logging.info(f"{gene_name} from {cell_line} passed filtering")
            with open("data/gene_list_from_filtering.txt", "a") as f:
                # Only add a newline if the file is not empty
                if f.tell() != 0:
                    f.write("\n")
                f.write(gene_name)
        else:
            logging.info(f"{gene_name} from {cell_line} was removed in filtering")

    else:
        plotter(
            plot_df,
            lowest_y_value,
            highest_y_value,
            cell_line,
            gene_name,
            output_dir,
            variants,
        )


def filter_proteins(line_dictionary, concentrations, variants, threshold=2, flat=False):
    # NOTE: Concentrations are treated as categorical data to equally weight each concentration change
    logging.info("filtering started")
    logging.info(f"Area between curves threshold: {threshold}")

    # Convert to arrays
    y1 = np.array(line_dictionary[variants[0]])
    y2 = np.array(line_dictionary[variants[1]])
    x = np.arange(len(concentrations))  # [0, 1, 2, 3] for equidistant x values

    # Compute absolute differences and area
    abs_diff = np.abs(y1 - y2)
    area = np.trapezoid(abs_diff, x)

    if flat:
        # Check that control curve is flat
        flatness_threshold = 0.35
        logging.info(
            f"Standard deviation/Control curve flatness threshold: {flatness_threshold}"
        )
        control_curve_name = "R"
        logging.info(f"control curve name: {control_curve_name}")
        ctrl_curve = np.array(line_dictionary[control_curve_name])  # change R as needed
        std_dev = np.std(ctrl_curve)
        is_flat = std_dev < flatness_threshold

        if area > threshold and is_flat:
            return True
    else:
        if area > threshold:
            return True
        return False


def main():
    # initiate argparse parser
    args = parser.parse_args()

    # data file
    data = args.data

    # output folder
    out_dir = args.output_folder

    # convert data file to dataframe
    data_df, y_max, y_min = read_data(data)

    # Mode 1: Single gene plotting
    if args.gene:
        logging.info(f"Plotting single gene: {args.gene}")
        filter_or_plot_gene(args.gene, data_df, out_dir, y_max, y_min)

    # Mode 2: Gene list plotting
    elif args.gene_list:
        gene_list = read_gene_list(args.gene_list)
        available_genes = set(data_df["Genes"].unique())
        logging.info(f"plotting {len(available_genes)} genes from list...")
        found_genes = [gene for gene in gene_list if gene in available_genes]
        missing_genes = [gene for gene in gene_list if gene not in available_genes]
        logging.info(f"Found {len(found_genes)} genes, {len(missing_genes)} missing")

        filtered_df = data_df[data_df["Genes"].isin(found_genes)]
        y_max = filtered_df["AVG Log2 Ratio"].max()
        y_min = filtered_df["AVG Log2 Ratio"].min()

        if missing_genes:
            logging.info(f"Missing genes: {missing_genes[:5]}...")  # Show first 5

        for i, gene_name in enumerate(found_genes, 1):
            logging.info(f"Plotting {i}/{len(found_genes)}: {gene_name}")
            filter_or_plot_gene(gene_name, data_df, out_dir, y_max, y_min)

    # Mode 3: Filter proteins and create a gene list of filtered
    elif args.filter:
        genes = data_df["Genes"].unique()
        logging.info("Applying filter and adding genes to gene list...")
        for i, gene_name in enumerate(genes, 1):
            logging.info(f"protein {i}/{len(genes)}")
            filter_or_plot_gene(gene_name, data_df, out_dir, y_max, y_min, True, True)

    # Mode 4: Plot all proteins
    else:
        genes = data_df["Genes"].unique()
        logging.info("Plotting all genes...")
        for i, gene_name in enumerate(genes, 1):
            logging.info(f"protein {i}/{len(genes)}")
            filter_or_plot_gene(gene_name, data_df, out_dir, y_max, y_min)


if __name__ == "__main__":
    main()
