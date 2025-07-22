import logging
from pathlib import Path
from typing import Tuple, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Configure logging to stdout with basic settings
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

# --- Interactive Helper Functions ---

def prompt_for_path(prompt_text: str, check_for_file: bool) -> Path:
    """Repeatedly prompts the user for a path until a valid one is given."""
    while True:
        path_str = input(prompt_text).strip()
        path = Path(path_str)
        if not path.exists():
            logging.error(f"Error: The path '{path_str}' does not exist. Please try again.")
            continue
        if check_for_file and not path.is_file():
            logging.error(f"Error: The path '{path_str}' is not a file. Please try again.")
            continue
        if not check_for_file and not path.is_dir():
            logging.error(f"Error: The path '{path_str}' is not a directory. Please try again.")
            continue
        return path

def ask_yes_no(prompt_text: str) -> bool:
    """Asks a yes/no question and returns a boolean."""
    while True:
        answer = input(f"{prompt_text} (y/n): ").strip().lower()
        if answer in ["y", "yes"]:
            return True
        if answer in ["n", "no"]:
            return False
        logging.warning("Invalid input. Please enter 'y' or 'n'.")

# --- Core Logic (Largely Unchanged) ---

def read_gene_list(file_path: Path) -> list:
    """Read gene names from a text file, one per line."""
    try:
        with open(file_path, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        logging.info(f"{len(genes)} genes read from gene list")
        return genes
    except Exception as e:
        logging.error(f"Error reading {file_path}:\n{e}")
        return []

def read_data(path: Path) -> Optional[Tuple[pd.DataFrame, float, float]]:
    """Reads and performs initial filtering on the data file."""
    try:
        data_df = pd.read_csv(path, delimiter="\t")
        filtered_df = data_df[data_df["# Unique Total Peptides"] >= 2]
        if filtered_df.empty:
            logging.error("No data remains after filtering for '# Unique Total Peptides' >= 2.")
            return None
        highest_fld_chng = filtered_df["AVG Log2 Ratio"].max()
        lowest_fld_chng = filtered_df["AVG Log2 Ratio"].min()
        logging.info(f"{path} read success")
        return filtered_df, highest_fld_chng, lowest_fld_chng
    except Exception as e:
        logging.error(f"Error reading or processing {path}:\n{e}")
        return None

def plotter(df, y_min, y_max, cell_line, gene_name, output_dir, variants, error_df=None):
    """Generates and saves a single plot."""
    try:
        sns.set_style("whitegrid")
        plt.rcParams["font.size"] = 11

        sorted_variants = sorted(variants)
        if len(sorted_variants) < 2:
            logging.error(f"Cannot plot {gene_name}, requires at least 2 variants, found {len(sorted_variants)}")
            return

        variant_palette = {sorted_variants[0]: "#1f77b4", sorted_variants[1]: "#d62728"}

        sns.lineplot(
            data=df,
            x="concentration",
            y="log_fld_change",
            hue="variant",
            palette=variant_palette,
            hue_order=sorted_variants,
            marker="o",
            linewidth=2,
            markersize=6,
        )

        if error_df is not None:
            ax = plt.gca()
            unique_concentrations = df["concentration"].unique()
            concentration_to_x = {conc: i for i, conc in enumerate(unique_concentrations)}
            for variant in sorted_variants:
                variant_data = df[df["variant"] == variant]
                error_data = error_df[error_df["variant"] == variant]
                x_positions = [concentration_to_x[conc] for conc in variant_data["concentration"]]
                ax.errorbar(
                    x=x_positions,
                    y=variant_data["log_fld_change"],
                    yerr=error_data["Standard Error"],
                    fmt="none",
                    color=variant_palette[variant],
                    capsize=4,
                )

        plt.ylim(y_min, y_max)
        plt.xlabel("Concentration", fontsize=12)
        plt.ylabel("AVG Log2 Ratio", fontsize=12)
        plt.title(f"{cell_line} {gene_name}", fontsize=14, fontweight="bold", pad=20)
        plt.axhline(y=0, color="gray", linewidth=0.8, alpha=0.7)
        plt.legend(frameon=True, fancybox=True, shadow=True)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{gene_name}_{cell_line}.png", dpi=300, bbox_inches="tight")
        plt.close()
        logging.info(f"{gene_name} from {cell_line} plotted")
    except Exception as e:
        logging.error(f"Error plotting {gene_name} from {cell_line}:\n{e}")

def get_lines(variants, gene_data, sorted_concentrations, in_column, out_column):
    """Reshapes data for plotting by extracting values for each variant curve."""
    lines = {}
    for variant in variants:
        lines[variant] = []
        var_mask = gene_data["Comparison (group1/group2)"].str.split("_").str[1] == variant
        variant_df = gene_data[var_mask]
        for concentration in sorted_concentrations:
            try:
                row_mask = variant_df["Comparison (group1/group2)"].str.contains(concentration, na=False)
                column_value = variant_df.loc[row_mask, in_column].iloc[0]
                lines[variant].append(column_value)
            except IndexError:
                logging.warning(f"Missing data for variant '{variant}' at concentration '{concentration}'")
                lines[variant].append(np.nan)

    plot_df = pd.DataFrame(lines, index=sorted_concentrations).reset_index()
    plot_df = plot_df.melt(id_vars="index", var_name="variant", value_name=out_column)
    plot_df.rename(columns={"index": "concentration"}, inplace=True)
    return plot_df.dropna(), lines


def filter_proteins(line_dictionary: dict, variants: list, flat: bool = False) -> bool:
    """Filters proteins based on area between curves and control curve flatness."""
    threshold = 2.0
    y1 = np.array(line_dictionary[variants[0]])
    y2 = np.array(line_dictionary[variants[1]])

    # Ignore NaN values for calculation
    valid_indices = ~np.isnan(y1) & ~np.isnan(y2)
    if np.sum(valid_indices) < 2: # Need at least 2 points to calculate area
        return False

    y1, y2 = y1[valid_indices], y2[valid_indices]
    x = np.arange(len(y1))

    area = np.trapz(np.abs(y1 - y2), x)

    if not flat:
        return area > threshold

    # Check that control curve is flat
    flatness_threshold = 0.35
    control_curve_name = "R"  # NOTE: Assuming 'R' is the control variant name
    if control_curve_name not in line_dictionary:
        logging.warning(f"Control variant '{control_curve_name}' not found for flatness check.")
        return False

    ctrl_curve = np.array(line_dictionary[control_curve_name])
    ctrl_curve = ctrl_curve[~np.isnan(ctrl_curve)] # remove NaNs
    if len(ctrl_curve) < 2:
        return False

    is_flat = np.std(ctrl_curve) < flatness_threshold
    return area > threshold and is_flat


def filter_or_plot_gene(gene_name, data_df, output_dir, highest_y, lowest_y, do_filter=False, filter_flat=False, error_bars=False):
    """Main processing function for a single gene to either filter it or plot it."""
    gene_data = data_df[data_df["Genes"] == gene_name]
    if gene_data.empty:
        logging.warning(f"No data found for gene: {gene_name}")
        return

    try:
        cell_line = gene_data["Comparison (group1/group2)"].iloc[0].split("_")[0]
        variants = sorted(gene_data["Comparison (group1/group2)"].str.split("_").str[1].unique())
        concentrations = gene_data["Comparison (group1/group2)"].str.split("_").str[2].str.split(" ").str[0].unique()
        sorted_concentrations = sorted(concentrations, key=lambda x: float(x.rstrip("uM")))
    except IndexError:
        logging.error(f"Could not parse metadata for {gene_name}. Skipping.")
        return

    plot_df, lines = get_lines(variants, gene_data, sorted_concentrations, "AVG Log2 Ratio", "log_fld_change")
    if plot_df.empty:
        logging.warning(f"No valid data to process for gene {gene_name}.")
        return

    if do_filter:
        if filter_proteins(lines, variants, flat=filter_flat):
            logging.info(f"{gene_name} from {cell_line} PASSED filtering")
            with open(output_dir / "gene_list_from_filtering.txt", "a") as f:
                f.write(f"{gene_name}\n")
        else:
            logging.info(f"{gene_name} from {cell_line} was removed by filtering")
    else:
        error_df = None
        if error_bars:
            error_df, _ = get_lines(variants, gene_data, sorted_concentrations, "Standard Error", "Standard Error")
        plotter(plot_df, lowest_y, highest_y, cell_line, gene_name, output_dir, variants, error_df)


# --- Main Application Logic ---

def main():
    """Runs the interactive command-line application."""
    print("--- Thermal Proteome Profiling (TPP) Plotter ---")

    # 1. Get input and output paths
    data_path = prompt_for_path("Enter the path to your data file (e.g., data.tsv): ", check_for_file=True)
    output_dir = prompt_for_path("Enter the path to your output folder: ", check_for_file=False)

    # 2. Read and process data
    data_tuple = read_data(data_path)
    if not data_tuple:
        return
    data_df, y_max_global, y_min_global = data_tuple

    # 3. Present menu and get user choice
    print("\nWhat would you like to do?")
    print("  1. Plot a single gene")
    print("  2. Plot a list of genes from a file")
    print("  3. Filter all proteins and create a gene list file")
    print("  4. Plot all proteins in the dataset (can be slow!)")

    while True:
        choice = input("Enter your choice (1-4): ").strip()
        if choice in ["1", "2", "3", "4"]:
            break
        logging.warning("Invalid choice. Please enter a number between 1 and 4.")

    # 4. Execute the chosen action
    if choice == "1":
        gene_name = input("Enter the gene name to plot: ").strip()
        error_bars = ask_yes_no("Add error bars to the plot?")
        filter_or_plot_gene(gene_name, data_df, output_dir, y_max_global, y_min_global, error_bars=error_bars)

    elif choice == "2":
        gene_list_path = prompt_for_path("Enter the path to your gene list file: ", check_for_file=True)
        gene_list = read_gene_list(gene_list_path)
        if not gene_list: return

        error_bars = ask_yes_no("Add error bars to the plots?")
        # Recalculate y-axis limits for just the genes in the list for better scaling
        list_df = data_df[data_df["Genes"].isin(gene_list)]
        y_max_list = list_df["AVG Log2 Ratio"].max()
        y_min_list = list_df["AVG Log2 Ratio"].min()

        for i, gene_name in enumerate(gene_list, 1):
            logging.info(f"Plotting {i}/{len(gene_list)}: {gene_name}")
            filter_or_plot_gene(gene_name, data_df, output_dir, y_max_list, y_min_list, error_bars=error_bars)

    elif choice == "3":
        logging.info("Applying filter and adding passed genes to 'gene_list_from_filtering.txt' in your output folder.")
        genes = data_df["Genes"].unique()
        for i, gene_name in enumerate(genes, 1):
            logging.info(f"Filtering protein {i}/{len(genes)}: {gene_name}")
            # The last two arguments enable filtering with the flatness check
            filter_or_plot_gene(gene_name, data_df, output_dir, y_max_global, y_min_global, do_filter=True, filter_flat=True)

    elif choice == "4":
        if not ask_yes_no("This will plot all proteins and may take a long time. Are you sure?"):
            print("Operation cancelled.")
            return
        genes = data_df["Genes"].unique()
        for i, gene_name in enumerate(genes, 1):
            logging.info(f"Plotting protein {i}/{len(genes)}: {gene_name}")
            filter_or_plot_gene(gene_name, data_df, output_dir, y_max_global, y_min_global)

    print("\n Operation complete. Check your output folder for results.")

if __name__ == "__main__":
    main()
