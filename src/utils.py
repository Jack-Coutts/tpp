import logging
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas.core.generic import NoReturn


# custom type converter for argparse using pathlib
# checks if the given path exists and is a file. Returns a Path object.
def check_filepath(filepath_str: str) -> Path:
    path = Path(filepath_str)
    if not path.exists():
        raise ValueError(f"Error: The path '{filepath_str}' does not exist.")
    if not path.is_file():
        raise ValueError(f"Error: The path '{filepath_str}' is not a file.")
    return path


# custom type converter for argparse using pathlib
# checks if the given path exists and is a directory. Returns a Path object.
def check_directory(dir_path_str: str) -> Path:
    path = Path(dir_path_str)
    if not path.exists():
        raise ValueError(f"Error: The path '{dir_path_str}' does not exist.")
    if not path.is_dir():
        raise ValueError(f"Error: The path '{dir_path_str}' is not a directory.")
    return path


# read gene names from a text file, one per line
def read_gene_list(file_path: Path) -> List:
    try:
        with open(file_path, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        logging.info(f"{len(genes)} genes read from gene list")
        return genes
    except Exception as e:
        logging.error(f"Error reading {file_path}:")
        logging.error(f"{e}")
        return []


# read in the data file
def read_data(path: Path) -> Tuple[pd.DataFrame, float, float]:
    try:
        data_df = pd.read_csv(path, delimiter="\t")
        filtered_df = data_df[data_df["# Unique Total Peptides"] >= 2]
        highest_fld_chng = filtered_df["AVG Log2 Ratio"].max()
        lowest_fld_chng = filtered_df["AVG Log2 Ratio"].min()
        logging.info(f"{path} read success")
        return filtered_df, highest_fld_chng, lowest_fld_chng
    except Exception as e:
        logging.error(f"Error reading {path}:")
        logging.error(f"{e}")


def get_variants(
    df: pd.DataFrame, single_mode: bool, line_name: str
) -> List[str] | NoReturn:
    variants = list(df["Comparison (group1/group2)"].str.split("_").str[1].unique())

    if single_mode:
        if line_name in variants:
            single_variant = [line_name]
            logging.info(f"one variant found successfuly: {line_name}")
            return single_variant
        else:
            logging.error(f"'{line_name}' not found in {variants}")
            raise AssertionError("Variant not found in data")
    else:
        if len(variants) == 2:
            logging.info(f"two variants found successfully: {variants}")
            return variants
        else:
            var_num = len(variants)
            logging.error(f"{var_num} variants found but only two were expected")
            raise AssertionError("Incorrect number of variants found in data")


def plotter(
    df, y_min, y_max, cell_line, gene_name, output_dir, variants, error_df=None
):
    import matplotlib

    # static backend for rendering pngs
    # since only writing graphs to files
    matplotlib.use("agg")

    try:
        # Set seaborn style
        sns.set_style("whitegrid")
        plt.rcParams["font.size"] = 11

        # ensure line colour consistency
        if len(variants) > 1:
            try:
                sorted_variants = sorted(variants, key=lambda x: x[0])
                variant_palette = {
                    sorted_variants[0]: "#1f77b4",
                    sorted_variants[1]: "#d62728",
                }
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
            except Exception as e:
                logging.error(
                    f"Less than 2 protein variants exist: Variants: {variants}\n Error: {e}"
                )
                return
        else:
            sorted_variants = variants
            variant_palette = {sorted_variants[0]: "#d62728"}

            sns.lineplot(
                data=df,
                x="concentration",  # categorical, equally spaced
                y="log_fld_change",
                hue="variant",  # separate lines per variant
                palette=variant_palette,
                marker="o",
                linewidth=2,
                markersize=6,
            )

        # Add error bars if error_df is provided
        if error_df is not None:
            ax = plt.gca()
            unique_concentrations = df["concentration"].unique()
            concentration_to_x = {
                conc: i for i, conc in enumerate(unique_concentrations)
            }

            for variant in sorted_variants:
                variant_data = df[df["variant"] == variant]
                error_data = error_df[error_df["variant"] == variant]
                x_positions = [
                    concentration_to_x[conc] for conc in variant_data["concentration"]
                ]

                ax.errorbar(
                    x=x_positions,
                    y=variant_data["log_fld_change"],
                    yerr=error_data["Standard Error"],
                    fmt="none",
                    color=variant_palette[variant],
                    capsize=4,
                    capthick=1.5,
                    alpha=0.8,
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
        plt.savefig(
            f"{output_dir}/{gene_name}_{cell_line}.png", dpi=300, bbox_inches="tight"
        )
        plt.close()

    except Exception as e:
        logging.error(f"Error plotting {gene_name} from {cell_line}:")
        logging.error(f"{e}")


def get_lines(variants, gene_data, sorted_concentrations, in_column, out_column):
    # check column exists
    if in_column not in gene_data.columns:
        raise ValueError(f"Column {in_column} not found in dataset")

    lines = {}
    for variant in variants:
        lines[variant] = []
        var_mask = (
            gene_data["Comparison (group1/group2)"].str.split("_").str[1] == variant
        )
        variant_df = gene_data[var_mask]

        for concentration in sorted_concentrations:
            try:
                column_value = variant_df.loc[
                    variant_df["Comparison (group1/group2)"].str.contains(
                        concentration, na=False
                    ),
                    in_column,
                ].iloc[0]
                lines[variant].append(column_value)
            except Exception as e:
                logging.error(f"Error getting line for {variant} at {concentration}:")
                logging.error(f"{e}")
                lines[variant].append(None)

    plot_df = pd.DataFrame(lines, index=sorted_concentrations).reset_index()
    plot_df = plot_df.melt(id_vars="index", var_name="variant", value_name=out_column)
    plot_df.rename(columns={"index": "concentration"}, inplace=True)

    return plot_df, lines


def plot_gene(
    gene_name,
    data_df,
    output_dir,
    variants,
    highest_y_value,
    lowest_y_value,
    error_bars=False,
):
    # Filter data for the specific gene
    gene_data = data_df[data_df["Genes"] == gene_name]

    # check gene data exists
    if gene_data.empty:
        print(f"No data found for gene: {gene_name}")
        return None

    # cell line name
    cell_line = gene_data["Comparison (group1/group2)"].iloc[0].split("_")[0]

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

    # Create line df
    plot_df, lines = get_lines(
        variants,
        gene_data,
        sorted_concentrations,
        "AVG Log2 Ratio",
        "log_fld_change",
    )

    error_df = None
    if error_bars:
        error_df, error_lines = get_lines(
            variants,
            gene_data,
            sorted_concentrations,
            "Standard Error",
            "Standard Error",
        )

    plotter(
        plot_df,
        lowest_y_value,
        highest_y_value,
        cell_line,
        gene_name,
        output_dir,
        variants,
        error_df,
    )


def filter_area_between_curves(
    line_dictionary, concentrations, variants, threshold=2, flat=False
):
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


def filter_for_flatness():
    pass


def filter_against_flatness():
    pass
