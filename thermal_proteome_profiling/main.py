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
    print(f"low: {lowest_fld_chng}, high: {highest_fld_chng}")

    return filtered_df, highest_fld_chng, lowest_fld_chng


def single_gene_tpp_plot(
    gene_name, data_df, output_dir, highest_y_value, lowest_y_value
):
    # Set seaborn style
    sns.set_style("whitegrid")
    plt.rcParams["font.size"] = 11

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

    plot_df = pd.DataFrame(lines, index=concentrations).reset_index()
    plot_df = plot_df.melt(
        id_vars="index", var_name="variant", value_name="log_fld_change"
    )
    plot_df.rename(columns={"index": "concentration"}, inplace=True)

    sns.lineplot(
        data=plot_df,
        x="concentration",  # categorical, equally spaced
        y="log_fld_change",
        hue="variant",  # separate lines per variant
        marker="o",
        linewidth=2,
        markersize=6,
    )

    plt.ylim(lowest_y_value, highest_y_value)
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


def filter_by_area_between_curves(R_df, S_df, threshold=0.3):
    """
    Filter genes based on area between R and S curves and return filtered DataFrames

    Parameters:
    R_df (pd.DataFrame): DataFrame containing R variant data
    S_df (pd.DataFrame): DataFrame containing S variant data
    threshold (float): Minimum area between curves to keep a gene

    Returns:
    tuple: (filtered_R_df, filtered_S_df)
    """
    concentration_columns = [
        col
        for col in R_df.columns
        if col not in ["PG.Genes", "PG.UniProtIds", "variant"]
    ]
    concentrations = [float(col) for col in concentration_columns]

    # Get all unique genes present in both dataframes
    R_genes = set(R_df["PG.Genes"].unique())
    S_genes = set(S_df["PG.Genes"].unique())
    all_genes = R_genes.union(S_genes)

    filtered_genes = []
    removed_genes = []

    for gene in all_genes:
        R_gene = R_df[R_df["PG.Genes"] == gene]
        S_gene = S_df[S_df["PG.Genes"] == gene]

        if not R_gene.empty and not S_gene.empty:
            R_values = R_gene[concentration_columns].iloc[0].values
            S_values = S_gene[concentration_columns].iloc[0].values

            # Calculate area between curves using trapezoidal rule
            area_between = np.trapezoid(np.abs(R_values - S_values), concentrations)

            if area_between >= threshold:
                filtered_genes.append(gene)
            else:
                removed_genes.append(gene)
        else:
            # Gene missing in one of the dataframes
            removed_genes.append(gene)

    # Filter the dataframes
    filtered_R_df = R_df[R_df["PG.Genes"].isin(filtered_genes)].copy()
    filtered_S_df = S_df[S_df["PG.Genes"].isin(filtered_genes)].copy()

    # Print statistics
    total_genes = len(all_genes)
    genes_kept = len(filtered_genes)
    genes_removed = len(removed_genes)

    print("Filtering Results:")
    print(f"Total genes: {total_genes}")
    print(f"Genes kept: {genes_kept}")
    print(f"Genes removed: {genes_removed}")
    print(f"Percentage kept: {genes_kept / total_genes * 100:.1f}%")
    print(f"Threshold used: {threshold}")

    # Optional: print some examples of removed genes
    if len(removed_genes) > 0:
        print(f"\nFirst 5 removed genes: {removed_genes[:5]}")

    return filtered_R_df, filtered_S_df


def plot_all_filtered(filtered_R, filtered_S, output_dir):
    R_genes = set(filtered_R["PG.Genes"].unique())

    print(f"Plotting {len(R_genes)} genes...")
    for i, gene_name in enumerate(R_genes, 1):
        print(f"Plotting {i}/{len(R_genes)}: {gene_name}")
        plot_r_s_lines(gene_name, filtered_R, filtered_S, output_dir)


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

    """
    # Mode 3: Filter and plot all (default behavior)
    else:
        print("Applying filter and plotting all passing genes...")
        filter_r, filter_s = filter_by_area_between_curves(r, s)
        plot_all_filtered(filter_r, filter_s, out_dir)
    """


if __name__ == "__main__":
    main()
