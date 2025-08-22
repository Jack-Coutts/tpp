import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")  # Set backend early for better performance
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Configuration constants
MIN_UNIQUE_PEPTIDES = 2
FLATNESS_THRESHOLD = 0.35
CONTROL_CURVE_NAME = "R"
DEFAULT_AREA_THRESHOLD = 2
DEFAULT_DPI = 300
DEFAULT_FONT_SIZE = 11
TITLE_FONT_SIZE = 14
LABEL_FONT_SIZE = 12
LINE_WIDTH = 2
MARKER_SIZE = 6
CAP_SIZE = 4
CAP_THICKNESS = 1.5
ERROR_BAR_ALPHA = 0.8
HORIZONTAL_LINE_ALPHA = 0.7
HORIZONTAL_LINE_WIDTH = 0.8

# Color palette for consistency
VARIANT_COLORS = {"primary": "#1f77b4", "secondary": "#d62728"}


def check_filepath(filepath_str: str) -> Path:
    """Custom type converter for argparse using pathlib.

    Checks if the given path exists and is a file.

    Args:
        filepath_str: String path to validate

    Returns:
        Path object if valid

    Raises:
        ValueError: If path doesn't exist or isn't a file
    """
    path = Path(filepath_str)
    if not path.exists():
        raise ValueError(f"Error: The path '{filepath_str}' does not exist.")
    if not path.is_file():
        raise ValueError(f"Error: The path '{filepath_str}' is not a file.")
    return path


def check_directory(dir_path_str: str) -> Path:
    """Custom type converter for argparse using pathlib.

    Checks if the given path exists and is a directory.

    Args:
        dir_path_str: String path to validate

    Returns:
        Path object if valid

    Raises:
        ValueError: If path doesn't exist or isn't a directory
    """
    path = Path(dir_path_str)
    if not path.exists():
        raise ValueError(f"Error: The path '{dir_path_str}' does not exist.")
    if not path.is_dir():
        raise ValueError(f"Error: The path '{dir_path_str}' is not a directory.")
    return path


def read_gene_list(file_path: Path) -> List[str]:
    """Read gene names from a text file, one per line.

    Args:
        file_path: Path to gene list file

    Returns:
        List of gene names (stripped of whitespace)
    """
    try:
        with open(file_path, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        logging.info(f"{len(genes)} genes read from gene list")
        return genes
    except FileNotFoundError:
        logging.error(f"Gene list file not found: {file_path}")
        return []
    except PermissionError:
        logging.error(f"Permission denied reading gene list: {file_path}")
        return []
    except UnicodeDecodeError:
        logging.error(f"Unable to decode gene list file: {file_path}")
        return []
    except Exception as e:
        logging.error(f"Unexpected error reading {file_path}: {e}")
        return []


def read_data(path: Path) -> Tuple[pd.DataFrame, float, float]:
    """Read and filter proteomics data file.

    Filters data to include only entries with at least MIN_UNIQUE_PEPTIDES peptides.

    Args:
        path: Path to TSV data file

    Returns:
        Tuple of (filtered_dataframe, max_log2_ratio, min_log2_ratio)
    """
    try:
        data_df = pd.read_csv(path, delimiter="\t")

        # Validate required columns
        required_cols = [
            "# Unique Total Peptides",
            "AVG Log2 Ratio",
            "Genes",
            "Comparison (group1/group2)",
        ]
        missing_cols = [col for col in required_cols if col not in data_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        filtered_df = data_df[data_df["# Unique Total Peptides"] >= MIN_UNIQUE_PEPTIDES]

        if filtered_df.empty:
            raise ValueError(
                f"No data remains after filtering for >= {MIN_UNIQUE_PEPTIDES} peptides"
            )

        highest_fld_chng = filtered_df["AVG Log2 Ratio"].max()
        lowest_fld_chng = filtered_df["AVG Log2 Ratio"].min()
        logging.info(f"{path} read success")
        return filtered_df, highest_fld_chng, lowest_fld_chng
    except FileNotFoundError:
        logging.error(f"Data file not found: {path}")
        raise
    except pd.errors.EmptyDataError:
        logging.error(f"Data file is empty: {path}")
        raise
    except pd.errors.ParserError as e:
        logging.error(f"Error parsing data file {path}: {e}")
        raise
    except ValueError as e:
        logging.error(f"Data validation error for {path}: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error reading {path}: {e}")
        raise


def get_variants(df: pd.DataFrame, single_mode: bool, line_name: str) -> List[str]:
    """Extract compound variants from comparison column.

    Args:
        df: Input dataframe with comparison data
        single_mode: If True, return only the specified line_name
        line_name: Name of specific line to extract (used in single_mode)

    Returns:
        List of variant names

    Raises:
        AssertionError: If expected number of variants not found
    """
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
    df: pd.DataFrame,
    y_min: float,
    y_max: float,
    cell_line: str,
    gene_name: str,
    output_dir: Path,
    variants: List[str],
    error_df: Optional[pd.DataFrame] = None,
) -> None:
    """Create and save thermal proteome profiling plot.

    Generates publication-quality line plots with consistent styling.

    Args:
        df: Plotting dataframe in long format
        y_min: Minimum y-axis value
        y_max: Maximum y-axis value
        cell_line: Cell line name for title
        gene_name: Gene name for title and filename
        output_dir: Directory to save plot
        variants: List of variant names
        error_df: Optional dataframe with error bar data
    """
    try:
        _setup_plot_style()
        variant_palette = _create_variant_palette(variants)
        
        if len(variants) < 1:
            logging.error(f"Insufficient protein variants for plot: {variants}")
            return
            
        _plot_lines(df, variants, variant_palette)
        
        if error_df is not None:
            _add_error_bars(df, error_df, variants, variant_palette)
        
        _format_plot(y_min, y_max, cell_line, gene_name)
        
        plt.savefig(
            f"{output_dir}/{gene_name}_{cell_line}.png",
            dpi=DEFAULT_DPI,
            bbox_inches="tight",
        )
        plt.close()

    except FileNotFoundError:
        logging.error(f"Output directory not found: {output_dir}")
    except PermissionError:
        logging.error(f"Permission denied writing to: {output_dir}")
    except ValueError as e:
        logging.error(f"Invalid plot data for {gene_name}: {e}")
    except Exception as e:
        logging.error(f"Unexpected error plotting {gene_name} from {cell_line}: {e}")


def _setup_plot_style() -> None:
    """Configure plot styling with consistent parameters."""
    sns.set_style("whitegrid")
    plt.rcParams["font.size"] = DEFAULT_FONT_SIZE


def _create_variant_palette(variants: List[str]) -> Dict[str, str]:
    """Create color palette for variants.

    Args:
        variants: List of variant names

    Returns:
        Dictionary mapping variant names to colors
    """
    if len(variants) == 1:
        return {variants[0]: VARIANT_COLORS["secondary"]}

    sorted_variants = sorted(variants, key=lambda x: x[0])
    return {
        sorted_variants[0]: VARIANT_COLORS["primary"],
        sorted_variants[1]: VARIANT_COLORS["secondary"],
    }


def _plot_lines(
    df: pd.DataFrame, variants: List[str], variant_palette: Dict[str, str]
) -> None:
    """Create the main line plot.

    Args:
        df: Plotting dataframe
        variants: List of variant names
        variant_palette: Color mapping for variants
    """
    if len(variants) > 1:
        sorted_variants = sorted(variants, key=lambda x: x[0])
        sns.lineplot(
            data=df,
            x="concentration",
            y="log_fld_change",
            hue="variant",
            palette=variant_palette,
            hue_order=sorted_variants,
            marker="o",
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
        )
    else:
        sns.lineplot(
            data=df,
            x="concentration",
            y="log_fld_change",
            hue="variant",
            palette=variant_palette,
            marker="o",
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
        )


def _add_error_bars(
    df: pd.DataFrame,
    error_df: pd.DataFrame,
    variants: List[str],
    variant_palette: Dict[str, str],
) -> None:
    """Add error bars to the plot.

    Args:
        df: Main plotting dataframe
        error_df: Error data dataframe
        variants: List of variant names
        variant_palette: Color mapping for variants
    """
    ax = plt.gca()
    unique_concentrations = df["concentration"].unique()
    concentration_to_x = {conc: i for i, conc in enumerate(unique_concentrations)}

    sorted_variants = (
        sorted(variants, key=lambda x: x[0]) if len(variants) > 1 else variants
    )

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
            capsize=CAP_SIZE,
            capthick=CAP_THICKNESS,
            alpha=ERROR_BAR_ALPHA,
        )


def _format_plot(y_min: float, y_max: float, cell_line: str, gene_name: str) -> None:
    """Apply final formatting to the plot.

    Args:
        y_min: Minimum y-axis value
        y_max: Maximum y-axis value
        cell_line: Cell line name for title
        gene_name: Gene name for title
    """
    plt.ylim(y_min, y_max)
    plt.xlabel("Concentration", fontsize=LABEL_FONT_SIZE)
    plt.ylabel("AVG Log2 Ratio", fontsize=LABEL_FONT_SIZE)
    plt.title(
        f"{cell_line} {gene_name}", fontsize=TITLE_FONT_SIZE, fontweight="bold", pad=20
    )
    plt.axhline(
        y=0, color="gray", linewidth=HORIZONTAL_LINE_WIDTH, alpha=HORIZONTAL_LINE_ALPHA
    )
    plt.legend(frameon=True, fancybox=True, shadow=True)
    plt.tight_layout()


def get_lines(
    variants: List[str],
    gene_data: pd.DataFrame,
    sorted_concentrations: List[str],
    in_column: str,
    out_column: str,
) -> Tuple[pd.DataFrame, Dict[str, List]]:
    """Extract line data for plotting from gene data.

    Transforms wide-format data into long format suitable for plotting.

    Args:
        variants: List of variant names
        gene_data: Filtered gene data
        sorted_concentrations: Concentrations in desired order
        in_column: Source column name
        out_column: Output column name

    Returns:
        Tuple of (melted_dataframe, raw_lines_dictionary)
    """
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
            except (IndexError, KeyError) as e:
                logging.warning(f"Missing data for {variant} at {concentration}: {e}")
                lines[variant].append(None)
            except Exception as e:
                logging.error(
                    f"Unexpected error getting line for {variant} at {concentration}: {e}"
                )
                lines[variant].append(None)

    plot_df = pd.DataFrame(lines, index=sorted_concentrations).reset_index()
    plot_df = plot_df.melt(id_vars="index", var_name="variant", value_name=out_column)
    plot_df.rename(columns={"index": "concentration"}, inplace=True)

    return plot_df, lines


def plot_gene(
    gene_name: str,
    data_df: pd.DataFrame,
    output_dir: Path,
    variants: List[str],
    highest_y_value: float,
    lowest_y_value: float,
    error_bars: bool = False,
) -> None:
    """Generate thermal proteome profiling plot for a specific gene.

    Creates a line plot showing log2 fold change vs concentration for each variant.

    Args:
        gene_name: Name of gene to plot
        data_df: Full dataset dataframe
        output_dir: Directory to save plot
        variants: List of variant names to include
        highest_y_value: Maximum y-axis value
        lowest_y_value: Minimum y-axis value
        error_bars: Whether to include error bars if available
    """
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
    line_dictionary: Dict[str, List[float]],
    concentrations: List[str],
    variants: List[str],
    threshold: float = DEFAULT_AREA_THRESHOLD,
    flat: bool = False,
) -> Optional[bool]:
    """Filter genes based on area between thermal curves.

    Calculates the area between two curves using trapezoidal integration.
    Optionally checks if control curve is sufficiently flat.

    Args:
        line_dictionary: Dictionary mapping variant names to y-values
        concentrations: List of concentration values (used as x-coordinates)
        variants: List of variant names
        threshold: Minimum area threshold for filtering
        flat: Whether to also check control curve flatness

    Returns:
        True if gene passes filter criteria, False otherwise, None if flat check fails
    """
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
        logging.info(
            f"Standard deviation/Control curve flatness threshold: {FLATNESS_THRESHOLD}"
        )
        logging.info(f"control curve name: {CONTROL_CURVE_NAME}")
        ctrl_curve = np.array(line_dictionary[CONTROL_CURVE_NAME])
        std_dev = np.std(ctrl_curve)
        is_flat = std_dev < FLATNESS_THRESHOLD

        if area > threshold and is_flat:
            return True
    else:
        if area > threshold:
            return True
        return False


def filter_for_flatness() -> None:
    """Filter genes based on control curve flatness.

    Placeholder for future implementation of flatness-based filtering.
    Will identify genes where control curves show minimal variation.
    """
    pass


def filter_against_flatness() -> None:
    """Filter genes that show significant control curve variation.

    Placeholder for future implementation of anti-flatness filtering.
    Will identify genes where control curves show substantial variation.
    """
    pass
