import argparse
from pathlib import Path
from statistics import mean

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
def read_data(path: Path) -> pd.DataFrame:
    data_df = pd.read_csv(path)

    filtered_df = data_df[
        data_df["PG.NrOfStrippedSequencesIdentified (Experiment-wide)"] >= 2
    ]

    return filtered_df


# group concentrations and protein variants
def combine_samples(data_df: pd.DataFrame) -> pd.DataFrame:
    # plot data
    plot_data_R = []
    plot_data_S = []

    # iterate over the rows - manually specifying columns due to inconsistencies
    for ind, row in data_df.iterrows():
        dmso_R = mean(
            [
                *row[data_df.filter(like="E01").columns].values,
                *row[data_df.filter(like="F01").columns].values,
                *row[data_df.filter(like="G01").columns].values,
                *row[data_df.filter(like="H01").columns].values,
            ]
        )
        dmso_S = mean(
            [
                *row[data_df.filter(like="E07").columns].values,
                *row[data_df.filter(like="F07").columns].values,
                *row[data_df.filter(like="G07").columns].values,
                *row[data_df.filter(like="H07").columns].values,
            ]
        )

        plot_data_R.append(
            {
                "PG.Genes": row["PG.Genes"],
                "PG.UniProtIds": row["PG.UniProtIds"],
                "variant": "R",
                "0.1": mean(
                    [
                        *row[data_df.filter(like="E02").columns].values,
                        *row[data_df.filter(like="F02").columns].values,
                        *row[data_df.filter(like="G02").columns].values,
                        *row[data_df.filter(like="H02").columns].values,
                    ]
                )
                - dmso_R,
                "0.3": mean(
                    [
                        *row[data_df.filter(like="E03").columns].values,
                        *row[data_df.filter(like="G03").columns].values,
                        *row[data_df.filter(like="H03").columns].values,
                    ]
                )
                - dmso_R,
                "1": mean(
                    [
                        *row[data_df.filter(like="E04").columns].values,
                        *row[data_df.filter(like="F04").columns].values,
                        *row[data_df.filter(like="G04").columns].values,
                        *row[data_df.filter(like="H04").columns].values,
                    ]
                )
                - dmso_R,
                "3": mean(
                    [
                        *row[data_df.filter(like="F05").columns].values,
                        *row[data_df.filter(like="G05").columns].values,
                        *row[data_df.filter(like="H05").columns].values,
                    ]
                )
                - dmso_R,
                "10": mean(
                    [
                        *row[data_df.filter(like="E06").columns].values,
                        *row[data_df.filter(like="F06").columns].values,
                        *row[data_df.filter(like="G06").columns].values,
                        *row[data_df.filter(like="H06").columns].values,
                    ]
                )
                - dmso_R,
            }
        )

        plot_data_S.append(
            {
                "PG.Genes": row["PG.Genes"],
                "PG.UniProtIds": row["PG.UniProtIds"],
                "variant": "S",
                "0.1": mean(
                    [
                        *row[data_df.filter(like="E08").columns].values,
                        *row[data_df.filter(like="F08").columns].values,
                        *row[data_df.filter(like="G08").columns].values,
                        *row[data_df.filter(like="H08").columns].values,
                    ]
                )
                - dmso_S,
                "0.3": mean(
                    [
                        *row[data_df.filter(like="E09").columns].values,
                        *row[data_df.filter(like="F09").columns].values,
                        *row[data_df.filter(like="G09").columns].values,
                        *row[data_df.filter(like="H09").columns].values,
                    ]
                )
                - dmso_S,
                "1": mean(
                    [
                        *row[data_df.filter(like="E10").columns].values,
                        *row[data_df.filter(like="F10").columns].values,
                        *row[data_df.filter(like="G10").columns].values,
                        *row[data_df.filter(like="H10").columns].values,
                    ]
                )
                - dmso_S,
                "3": mean(
                    [
                        *row[data_df.filter(like="E11").columns].values,
                        *row[data_df.filter(like="F11").columns].values,
                        *row[data_df.filter(like="G11").columns].values,
                        *row[data_df.filter(like="H11").columns].values,
                    ]
                )
                - dmso_S,
                "10": mean(
                    [
                        *row[data_df.filter(like="E12").columns].values,
                        *row[data_df.filter(like="F12").columns].values,
                        *row[data_df.filter(like="G12").columns].values,
                        *row[data_df.filter(like="H12").columns].values,
                    ]
                )
                - dmso_S,
            }
        )

    R_df = pd.DataFrame(plot_data_R)
    S_df = pd.DataFrame(plot_data_S)

    return R_df, S_df


def filter_by_area_between_curves(R_df, S_df, threshold=300):
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


def plot_r_s_lines(gene_name, R_df, S_df, output_dir):
    """
    Create a single seaborn-styled plot for one gene showing R and S variant response curves.

    Parameters:
    gene_name (str): The gene name to plot
    R_df (pd.DataFrame): DataFrame containing R variant data
    S_df (pd.DataFrame): DataFrame containing S variant data
    """

    # Set seaborn style
    sns.set_style("whitegrid")
    plt.rcParams["font.size"] = 11

    # Filter data for the specific gene
    R_gene = R_df[R_df["PG.Genes"] == gene_name]
    S_gene = S_df[S_df["PG.Genes"] == gene_name]

    if R_gene.empty and S_gene.empty:
        print(f"No data found for gene: {gene_name}")
        return None

    # Extract concentration columns
    concentration_columns = [
        col
        for col in R_df.columns
        if col not in ["PG.Genes", "PG.UniProtIds", "variant"]
    ]
    concentrations = sorted(concentration_columns, key=lambda x: float(x))

    # Create equally spaced x positions (0, 1, 2, 3, 4, etc.)
    x_positions = list(range(len(concentrations)))

    # Create concentration labels
    conc_labels = [f"{conc}uM" for conc in concentrations]

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot R variant (blue line)
    if not R_gene.empty:
        R_values = R_gene[concentrations].iloc[0].values
        plt.plot(
            x_positions,
            R_values,
            color="#4472C4",
            marker="o",
            linewidth=2,
            markersize=6,
            label="R",
        )

    # Plot S variant (red/orange line)
    if not S_gene.empty:
        S_values = S_gene[concentrations].iloc[0].values
        plt.plot(
            x_positions,
            S_values,
            color="#E15759",
            marker="o",
            linewidth=2,
            markersize=6,
            label="S",
        )

    # Formatting to match the reference image
    plt.title(gene_name, fontsize=14, fontweight="bold", pad=20)

    # Let y-axis limits be automatic based on data

    # Set x-axis ticks and labels with equal spacing
    plt.xticks(x_positions, conc_labels)

    # Add horizontal line at y=0
    plt.axhline(y=0, color="gray", linewidth=0.8, alpha=0.7)

    # Grid styling
    plt.grid(True, alpha=0.3, linewidth=0.5)

    # Labels
    plt.xlabel("Concentration", fontsize=12)
    # plt.ylabel("Response Value", fontsize=12)

    # Legend
    plt.legend(frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{gene_name}.png")
    plt.close()


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
    data_df = read_data(data)
    s, r = combine_samples(data_df)

    # Mode 1: Single gene plotting
    if args.gene:
        print(f"Plotting single gene: {args.gene}")
        plot_r_s_lines(args.gene, r, s, out_dir)

    # Mode 2: Gene list plotting
    elif args.gene_list:
        gene_list = read_gene_list(args.gene_list)
        print(f"Plotting {len(gene_list)} genes from list...")

        available_genes = set(r["PG.Genes"].unique())
        found_genes = [gene for gene in gene_list if gene in available_genes]
        missing_genes = [gene for gene in gene_list if gene not in available_genes]

        print(f"Found {len(found_genes)} genes, {len(missing_genes)} missing")
        if missing_genes:
            print(f"Missing genes: {missing_genes[:5]}...")  # Show first 5

        for i, gene_name in enumerate(found_genes, 1):
            print(f"Plotting {i}/{len(found_genes)}: {gene_name}")
            plot_r_s_lines(gene_name, r, s, out_dir)

    # Mode 3: Filter and plot all (default behavior)
    else:
        print("Applying filter and plotting all passing genes...")
        filter_r, filter_s = filter_by_area_between_curves(r, s)
        plot_all_filtered(filter_r, filter_s, out_dir)


if __name__ == "__main__":
    main()
