import argparse
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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

# Add this after your existing parser arguments, before the functions
mode_group = parser.add_mutually_exclusive_group()

mode_group.add_argument("-g", "--gene", type=str, help="Plot a single gene by name")

mode_group.add_argument(
    "-gl",
    "--gene-list",
    type=check_filepath,
    help="Path to a text file containing gene names (one per line)",
)


def read_gene_list(file_path: Path) -> list:
    """Read gene names from a text file, one per line"""
    with open(file_path, "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes


# read in the data file
def read_data(path: Path) -> Tuple[pd.DataFrame, float, float]:
    data_df = pd.read_csv(path, delimiter="\t")

    filtered_df = data_df[data_df["# Unique Total Peptides"] >= 2]
    highest_fld_chng = filtered_df["AVG Log2 Ratio"].max()
    lowest_fld_chng = filtered_df["AVG Log2 Ratio"].min()

    return filtered_df, highest_fld_chng, lowest_fld_chng


def plotter(df, y_min, y_max, cell_line, gene_name, output_dir):
    # Set seaborn style
    sns.set_style("whitegrid")
    plt.rcParams["font.size"] = 11

    sns.lineplot(
        data=df,
        x="concentration",  # categorical, equally spaced
        y="log_fld_change",
        hue="variant",  # separate lines per variant
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


def single_gene_tpp_plot(
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
            plotter(
                plot_df,
                lowest_y_value,
                highest_y_value,
                cell_line,
                gene_name,
                output_dir,
            )
            print(f"{gene_name}: passed filter & plotted.")
        else:
            print(f"{gene_name}: removed by filter.")

    else:
        plotter(
            plot_df, lowest_y_value, highest_y_value, cell_line, gene_name, output_dir
        )
        print(f"{gene_name}: plotted.")


def filter_proteins(line_dictionary, concentrations, variants, threshold=3, flat=False):
    # NOTE: Concentrations are treated as categorical data to equally weight each concentration change

    # Convert to arrays
    y1 = np.array(line_dictionary[variants[0]])
    y2 = np.array(line_dictionary[variants[1]])
    x = np.arange(len(concentrations))  # [0, 1, 2, 3] for equidistant x values

    # Compute absolute differences and area
    abs_diff = np.abs(y1 - y2)
    area = np.trapezoid(abs_diff, x)

    if flat:
        # Check that control curve is flat
        flatness_threshold = 0.3
        ctrl_curve = np.array(line_dictionary["R"])  # change R as needed
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
        print(f"Plotting single gene: {args.gene}")
        single_gene_tpp_plot(args.gene, data_df, out_dir, y_max, y_min)

    # Mode 2: Gene list plotting
    elif args.gene_list:
        gene_list = read_gene_list(args.gene_list)
        print(f"Plotting {len(gene_list)} genes from list...")

        available_genes = set(data_df["Genes"].unique())
        found_genes = [gene for gene in gene_list if gene in available_genes]
        missing_genes = [gene for gene in gene_list if gene not in available_genes]
        print(f"Found {len(found_genes)} genes, {len(missing_genes)} missing")

        filtered_df = data_df[data_df["Genes"].isin(found_genes)]
        y_max = filtered_df["AVG Log2 Ratio"].max()
        y_min = filtered_df["AVG Log2 Ratio"].min()

        if missing_genes:
            print(f"Missing genes: {missing_genes[:5]}...")  # Show first 5

        for i, gene_name in enumerate(found_genes, 1):
            print(f"Plotting {i}/{len(found_genes)}: {gene_name}")
            single_gene_tpp_plot(gene_name, data_df, out_dir, y_max, y_min)

    # Mode 3: Filter and plot all (default behavior)
    else:
        genes = data_df["Genes"]
        print("Applying filter and plotting all passing genes...")
        for i, gene_name in enumerate(genes, 1):
            single_gene_tpp_plot(gene_name, data_df, out_dir, y_max, y_min, True, True)


if __name__ == "__main__":
    main()
