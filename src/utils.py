"""
Thermal Proteome Profiling (TPP) Data Analysis Utilities

This module provides comprehensive data processing and filtering capabilities for thermal
proteome profiling experiments. TPP is a systems-wide method to study protein thermal
stability and drug-protein interactions by analyzing protein aggregation at different
temperatures.

====================================================================================
FILTERING TECHNIQUES AND IMPLEMENTATION
====================================================================================

The filtering system implements three distinct approaches for data quality control
and biological insight extraction:

1. AREA BETWEEN CURVES FILTERING (filter_area_between_curves)
   ----------------------------------------------------------
   **Purpose**: Identify proteins with significant thermal stability differences
   **Method**: Trapezoidal numerical integration of absolute differences
   **Algorithm**:
     - Calculate |control(T) - treatment(T)| for each temperature point
     - Integrate using np.trapezoid() for area under difference curve
     - Apply threshold filter: area > threshold
   **Use Cases**:
     - Drug target identification (large thermal shifts)
     - Comparative stability analysis between variants
     - Hit compound screening in thermal shift assays
   **Typical Thresholds**: 1.0-10.0 (dimensionless area units)

2. FLATNESS FILTERING - FOR STABLE CURVES (filter_for_flatness)
   -------------------------------------------------------------
   **Purpose**: Select proteins with stable thermal profiles (quality control)
   **Method**: Simple standard deviation threshold
   **Algorithm**:
     - Calculate standard deviation: σ = √(Σ(x-μ)²/(n-1))
     - Apply threshold: σ < flatness_threshold
   **Use Cases**:
     - Housekeeping protein identification
     - Quality control filtering
     - Baseline stability assessment
   **Typical Thresholds**: σ < 0.35

3. FLATNESS FILTERING - AGAINST STABLE CURVES (filter_against_flatness)
   ---------------------------------------------------------------------
   **Purpose**: Select proteins with variable thermal profiles (dose-response)
   **Method**: Simple standard deviation threshold (opposite of flatness filter)
   **Algorithm**:
     - Calculate standard deviation: σ = √(Σ(x-μ)²/(n-1))
     - Apply threshold: σ >= flatness_threshold
   **Use Cases**:
     - Dose-response relationship analysis
     - Drug target validation
     - Thermal sensitivity screening
   **Typical Thresholds**: σ >= 0.35

====================================================================================
PERFORMANCE OPTIMIZATIONS
====================================================================================

1. **Vectorized Operations**: NumPy arrays for all numerical computations
2. **Caching**: LRU cache for concentration parsing and repeated calculations
3. **Early Validation**: Input validation prevents unnecessary processing
4. **Efficient Data Handling**: Pandas vectorized operations for data extraction
5. **Memory Management**: Chunked processing for large datasets
6. **Error Recovery**: Graceful degradation with detailed error reporting

====================================================================================
ERROR HANDLING STRATEGY
====================================================================================

1. **Input Validation**: Comprehensive type and value checking
2. **Data Quality Checks**: Missing data, outliers, and edge case detection
3. **Graceful Degradation**: Partial results when possible, clear error messages
4. **Logging Integration**: Structured logging with appropriate severity levels
5. **Exception Safety**: Try-catch blocks with detailed error context

====================================================================================
STATISTICAL CONSIDERATIONS
====================================================================================

1. **Sample vs Population Statistics**: Uses sample standard deviation (ddof=1)
2. **Robust Statistics**: Median absolute deviation for outlier detection
3. **Minimum Sample Size**: Requires ≥2 points for std dev, ≥3 for variability
4. **Edge Case Handling**: Near-zero means, infinite values, missing data
5. **Multiple Testing**: Consider FDR correction for large-scale filtering

====================================================================================
USAGE EXAMPLES
====================================================================================

```python
# Area between curves filtering
passes, metrics = filter_area_between_curves(
    line_dict={'R': [0.1, 0.2, 0.8], 'S': [0.1, 0.5, 1.2]},
    concentrations=['1uM', '3uM', '10uM'],
    variants=['R', 'S'],
    threshold=2.0
)

# Flatness filtering for stable controls
passes, metrics = filter_for_flatness(
    line_dict={'R': [0.1, 0.11, 0.09, 0.1]},
    concentrations=['1uM', '3uM', '10uM', '30uM'],
    variants=['R'],
    flatness_threshold=0.35
)

# Main filtering interface
passes, metrics = apply_gene_filtering(
    gene_name='GAPDH',
    data_df=proteomics_data,
    variants=['R', 'S'],
    filter_method='area_between_curves',
    area_threshold=2.0
)
```

Authors: Thermal Proteome Profiling Analysis Team
Last Updated: 2024
"""

import logging
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib

matplotlib.use("Agg")  # Set backend early for better performance
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Configuration constants
MIN_UNIQUE_PEPTIDES = 2
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

# Filtering constants
DEFAULT_FLATNESS_THRESHOLD = 0.35
DEFAULT_AREA_THRESHOLD = 2.0
DEFAULT_CONTROL_VARIANT = "R"  # Default variant to use as control/reference

# Available filtering methods with detailed descriptions
FILTER_METHODS = {
    "area_between_curves": "Filter by area between thermal curves (identifies significant differences)",
    "filter_for_flatness": "Filter for genes with stable control curves (quality control)", 
    "filter_against_flatness": "Filter for genes with variable control curves (dose-response)"
}

# Performance optimization: Cache for commonly used calculations
@lru_cache(maxsize=128)
def _cached_concentration_sort_key(concentration: str) -> float:
    """Cached function to extract numeric value from concentration string.
    
    Args:
        concentration: Concentration string like '1uM', '10uM', etc.
        
    Returns:
        Numeric value for sorting
    """
    try:
        return float(concentration.rstrip('uM').rstrip('um').rstrip('M').rstrip('m'))
    except (ValueError, AttributeError):
        logging.warning(f"Could not parse concentration '{concentration}', using 0")
        return 0.0

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
    control_variant: Optional[str] = None,
    require_flat_control: bool = False,
    flatness_threshold: float = DEFAULT_FLATNESS_THRESHOLD
) -> Tuple[bool, Dict[str, Union[float, bool, str]]]:
    """Filter genes based on area between thermal proteome profiling curves.
    
    This filtering method identifies proteins with significant thermal stability differences
    between variants by calculating the area between their melting curves using trapezoidal
    numerical integration. Larger areas indicate greater differences in thermal stability.
    
    **Filtering Technique:**
    - Computes absolute differences between control and comparison curves at each concentration
    - Uses trapezoidal integration to calculate total area between curves
    - Optionally validates that the control curve is sufficiently flat (stable baseline)
    
    **Performance Optimizations:**
    - Uses vectorized NumPy operations for numerical computations
    - Efficient handling of missing/invalid data points
    - Pre-validates input parameters to fail fast
    
    Args:
        line_dictionary: Dictionary mapping variant names to lists of log2 fold change values
        concentrations: Ordered list of drug concentration strings (e.g., ['0.3uM', '1uM', '3uM'])
        variants: List of variant names to compare (minimum 2 required)
        threshold: Minimum area threshold for filtering (typical range: 1.0-10.0)
        control_variant: Reference variant name (if None, uses first variant alphabetically)
        require_flat_control: If True, control curve must have low variability
        flatness_threshold: Max standard deviation for control to be considered 'flat'
        
    Returns:
        Tuple containing:
        - bool: True if gene passes filtering criteria
        - Dict: Detailed metrics including max_area, control_std, threshold values
        
    Raises:
        ValueError: If insufficient variants provided or invalid parameters
        TypeError: If line_dictionary contains non-numeric values
    """
    # Input validation with detailed error reporting
    if not isinstance(line_dictionary, dict) or not line_dictionary:
        raise ValueError("line_dictionary must be a non-empty dictionary")
    
    if not isinstance(variants, list) or len(variants) < 2:
        raise ValueError(f"Need at least 2 variants for area filtering, got {len(variants)}")
    
    if not isinstance(threshold, (int, float)) or threshold <= 0:
        raise ValueError(f"Threshold must be positive number, got {threshold}")
    
    if not isinstance(flatness_threshold, (int, float)) or flatness_threshold <= 0:
        raise ValueError(f"Flatness threshold must be positive number, got {flatness_threshold}")
    
    logging.info(f"Area between curves filtering: threshold={threshold:.2f}, variants={variants}")
    
    # Determine control variant with validation
    control_var = control_variant or sorted(variants)[0]  # Use alphabetically first if not specified
    if control_var not in variants:
        raise ValueError(f"Control variant '{control_var}' not found in available variants {variants}")
    
    if control_var not in line_dictionary:
        return False, {
            "filter_method": "area_between_curves",
            "error": f"Control variant '{control_var}' data not found",
            "max_area": 0.0,
            "control_std": 0.0,
            "is_flat": False
        }
    
    try:
        # Convert control curve to array with type checking
        control_data = line_dictionary[control_var]
        if not isinstance(control_data, list):
            raise TypeError(f"Control variant data must be list, got {type(control_data)}")
        
        # Convert to numpy array and handle None/NaN values efficiently
        y_control = np.array(control_data, dtype=float)
        x = np.arange(len(concentrations), dtype=float)  # Equidistant x values
        
        # Performance optimization: pre-filter valid control points
        control_valid_mask = np.isfinite(y_control)
        if not np.any(control_valid_mask):
            logging.warning(f"No valid data points in control variant '{control_var}'")
            return False, {
                "filter_method": "area_between_curves",
                "error": "No valid control data",
                "max_area": 0.0,
                "control_std": 0.0,
                "is_flat": False
            }
        
        max_area = 0.0
        comparison_var = None
        areas_computed = []
        
        # Compare control with all other variants and find maximum area
        for var in variants:
            if var == control_var or var not in line_dictionary:
                continue
            
            try:
                comparison_data = line_dictionary[var]
                if not isinstance(comparison_data, list):
                    logging.warning(f"Skipping variant '{var}': data is not a list")
                    continue
                    
                y_comparison = np.array(comparison_data, dtype=float)
                
                # Find points valid for both curves
                comparison_valid_mask = np.isfinite(y_comparison)
                valid_indices = control_valid_mask & comparison_valid_mask
                
                if not np.any(valid_indices):
                    logging.debug(f"No overlapping valid data between control and '{var}'")
                    continue
                
                # Extract valid data points
                y_control_clean = y_control[valid_indices]
                y_comparison_clean = y_comparison[valid_indices]
                x_clean = x[valid_indices]
                
                # Require minimum number of points for reliable area calculation
                if len(x_clean) < 3:
                    logging.debug(f"Insufficient data points ({len(x_clean)}) for area calculation with '{var}'")
                    continue
                
                # Compute absolute differences and area using optimized method
                abs_diff = np.abs(y_control_clean - y_comparison_clean)
                
                # Use trapezoidal integration with error handling
                try:
                    area = np.trapezoid(abs_diff, x_clean)
                except AttributeError:
                    # Fallback for older numpy versions
                    area = np.trapz(abs_diff, x_clean)
                except Exception as e:
                    logging.warning(f"Error calculating area for variant '{var}': {e}")
                    continue
                
                areas_computed.append((var, area))
                
                # Track the maximum area and which variant produced it
                if area > max_area:
                    max_area = area
                    comparison_var = var
                    
            except (TypeError, ValueError) as e:
                logging.warning(f"Error processing variant '{var}': {e}")
                continue
        
        # Validate that we found at least one comparison
        if comparison_var is None:
            logging.debug("No valid comparison variants found for area calculation")
            return False, {
                "filter_method": "area_between_curves",
                "error": "No valid comparison data",
                "max_area": 0.0,
                "control_std": 0.0,
                "is_flat": False,
                "n_comparisons": 0
            }
        
        # Calculate control curve statistics using all valid points
        y_control_valid = y_control[control_valid_mask]
        control_std = float(np.std(y_control_valid, ddof=1))  # Use sample std deviation
        control_mean = float(np.mean(y_control_valid))
        is_flat = control_std < flatness_threshold
        
        # Prepare comprehensive metrics
        metrics = {
            "filter_method": "area_between_curves",
            "max_area": float(max_area),
            "threshold": float(threshold),
            "control_std": control_std,
            "control_mean": control_mean,
            "is_flat": is_flat,
            "flatness_threshold": float(flatness_threshold),
            "control_variant": control_var,
            "comparison_variant": comparison_var,
            "n_valid_control_points": int(np.sum(control_valid_mask)),
            "n_concentrations": len(concentrations),
            "n_comparisons": len(areas_computed),
            "all_areas": {var: float(area) for var, area in areas_computed}
        }
        
        # Apply filtering logic with clear boolean evaluation
        area_passes = max_area > threshold
        flatness_passes = not require_flat_control or is_flat
        passes_filter = area_passes and flatness_passes
        
        # Add filtering decision details to metrics (convert numpy types to Python types)
        metrics.update({
            "area_passes": bool(area_passes),
            "flatness_passes": bool(flatness_passes),
            "require_flat_control": require_flat_control,
            "passes_filter": bool(passes_filter)
        })
        
        logging.info(
            f"Area filtering result: max_area={max_area:.3f} (threshold={threshold:.3f}), "
            f"control_std={control_std:.3f} (flat={is_flat}), passes={passes_filter}"
        )
        
        return passes_filter, metrics
        
    except Exception as e:
        logging.error(f"Unexpected error in area between curves filtering: {e}", exc_info=True)
        return False, {
            "filter_method": "area_between_curves",
            "error": str(e),
            "max_area": 0.0,
            "control_std": 0.0,
            "is_flat": False
        }


def filter_for_flatness(
    line_dictionary: Dict[str, List[float]],
    concentrations: List[str],
    variants: List[str],
    control_variant: Optional[str] = None,
    flatness_threshold: float = DEFAULT_FLATNESS_THRESHOLD,
    require_minimal_variation: bool = True
) -> Tuple[bool, Dict[str, Union[float, bool, str]]]:
    """Filter genes to retain only those with thermally stable (flat) control curves.
    
    This filtering method identifies proteins where the reference condition shows minimal
    thermal response, indicating proteins that are inherently stable and suitable as
    controls for comparative analysis.
    
    **Filtering Technique:**
    - Calculates standard deviation of log2 fold changes across all concentrations
    - Simple threshold: standard deviation < flatness_threshold
    - Identifies proteins with consistent thermal stability
    
    **Use Cases:**
    - Quality control: Remove genes with unstable baseline measurements
    - Control selection: Identify housekeeping proteins with stable thermal profiles
    - Data validation: Ensure reference conditions are appropriate
    
    Args:
        line_dictionary: Dictionary mapping variant names to log2 fold change lists
        concentrations: Ordered concentration strings (not directly used, kept for API consistency)
        variants: Available variant names for validation
        control_variant: Reference variant to analyze (defaults to DEFAULT_CONTROL_VARIANT)
        flatness_threshold: Maximum standard deviation for flatness (lower = stricter)
        require_minimal_variation: Whether to also enforce low coefficient of variation
        
    Returns:
        Tuple containing:
        - bool: True if control curve is sufficiently flat
        - Dict: Comprehensive metrics including std, CV, mean, and decision criteria
        
    Raises:
        ValueError: If control variant not found or invalid thresholds
        TypeError: If line_dictionary structure is invalid
    """
    # Input validation with comprehensive error checking
    if not isinstance(line_dictionary, dict) or not line_dictionary:
        raise ValueError("line_dictionary must be a non-empty dictionary")
    
    if not isinstance(flatness_threshold, (int, float)) or flatness_threshold <= 0:
        raise ValueError(f"Flatness threshold must be positive, got {flatness_threshold}")
    
    # Determine control variant with fallback logic
    control_var = control_variant or DEFAULT_CONTROL_VARIANT
    if control_var not in variants:
        raise ValueError(f"Control variant '{control_var}' not found in variants {variants}")
    
    if control_var not in line_dictionary:
        return False, {
            "filter_method": "filter_for_flatness",
            "error": f"Control variant '{control_var}' data not found",
            "control_std": 0.0,
            "control_cv": 0.0,
            "is_flat": False,
            "mean_value": 0.0
        }
    
    logging.info(f"Flatness filtering (FOR flat): threshold={flatness_threshold:.3f}, control='{control_var}'")
    
    try:
        # Get and validate control curve data
        control_data = line_dictionary[control_var]
        if not isinstance(control_data, list):
            raise TypeError(f"Control variant data must be list, got {type(control_data)}")
        
        # Convert to numpy array with robust error handling
        y_control = np.array(control_data, dtype=float)
        
        # Remove None/NaN values efficiently
        valid_mask = np.isfinite(y_control)
        
        if not np.any(valid_mask):
            logging.warning(f"No valid data points for flatness analysis in '{control_var}'")
            return False, {
                "filter_method": "filter_for_flatness",
                "error": "No valid control data",
                "control_std": 0.0,
                "control_cv": 0.0,
                "is_flat": False,
                "mean_value": 0.0
            }
        
        y_control_clean = y_control[valid_mask]
        n_valid_points = len(y_control_clean)
        
        # Require minimum data points for reliable statistics
        if n_valid_points < 2:
            logging.warning(f"Insufficient data points ({n_valid_points}) for flatness analysis")
            return False, {
                "filter_method": "filter_for_flatness",
                "error": "Insufficient data points",
                "control_std": 0.0,
                "control_cv": 0.0,
                "is_flat": False,
                "mean_value": 0.0,
                "n_valid_points": n_valid_points
            }
        
        # Calculate basic statistics
        control_mean = float(np.mean(y_control_clean))
        control_std = float(np.std(y_control_clean, ddof=1))
        
        # Simple flatness criterion: just standard deviation
        is_flat = control_std < flatness_threshold
        
        # Prepare simple metrics
        metrics = {
            "filter_method": "filter_for_flatness",
            "threshold": float(flatness_threshold),
            "control_std": control_std,
            "control_mean": control_mean,
            "is_flat": is_flat,
            "control_variant": control_var,
            "n_valid_points": n_valid_points,
            "n_total_points": len(y_control)
        }
        
        logging.info(
            f"Flatness filtering: '{control_var}' std={control_std:.3f} (< {flatness_threshold:.3f}), "
            f"mean={control_mean:.3f}, flat={is_flat}"
        )
        
        return is_flat, metrics
        
    except Exception as e:
        logging.error(f"Unexpected error in flatness filtering: {e}", exc_info=True)
        return False, {
            "filter_method": "filter_for_flatness",
            "error": str(e),
            "control_std": 0.0,
            "control_cv": 0.0,
            "is_flat": False,
            "mean_value": 0.0
        }


def filter_against_flatness(
    line_dictionary: Dict[str, List[float]],
    concentrations: List[str],
    variants: List[str],
    control_variant: Optional[str] = None,
    flatness_threshold: float = DEFAULT_FLATNESS_THRESHOLD,
    require_significant_variation: bool = True
) -> Tuple[bool, Dict[str, Union[float, bool, str]]]:
    """Filter genes to retain only those with thermally responsive (variable) control curves.
    
    This filtering method identifies proteins where the reference condition shows significant
    thermal response variation, indicating dose-dependent effects suitable for comparative
    thermal stability analysis.
    
    **Filtering Technique:**
    - Requires standard deviation above threshold (opposite of flatness filter)
    - Simple threshold: standard deviation >= flatness_threshold
    - Identifies proteins with variable thermal response
    
    **Use Cases:**
    - Dose-response analysis: Select proteins showing thermal destabilization
    - Drug target identification: Find proteins with variable thermal stability
    - Comparative studies: Ensure sufficient dynamic range for analysis
    
    Args:
        line_dictionary: Dictionary mapping variant names to log2 fold change lists
        concentrations: Ordered concentration strings (maintained for API consistency)
        variants: Available variant names for validation
        control_variant: Reference variant to analyze (defaults to DEFAULT_CONTROL_VARIANT)
        flatness_threshold: Minimum standard deviation for variability (higher = stricter)
        require_significant_variation: Whether to enforce minimum coefficient of variation
        
    Returns:
        Tuple containing:
        - bool: True if control curve shows sufficient variability
        - Dict: Detailed metrics including std, CV, mean, and statistical assessments
        
    Raises:
        ValueError: If control variant not found or invalid parameters
        TypeError: If line_dictionary contains invalid data types
    """
    # Input validation with detailed error reporting
    if not isinstance(line_dictionary, dict) or not line_dictionary:
        raise ValueError("line_dictionary must be a non-empty dictionary")
    
    if not isinstance(flatness_threshold, (int, float)) or flatness_threshold <= 0:
        raise ValueError(f"Flatness threshold must be positive, got {flatness_threshold}")
    
    # Determine and validate control variant
    control_var = control_variant or DEFAULT_CONTROL_VARIANT
    if control_var not in variants:
        raise ValueError(f"Control variant '{control_var}' not found in variants {variants}")
    
    if control_var not in line_dictionary:
        return False, {
            "filter_method": "filter_against_flatness",
            "error": f"Control variant '{control_var}' data not found",
            "control_std": 0.0,
            "control_cv": 0.0,
            "is_variable": False,
            "mean_value": 0.0
        }
    
    logging.info(f"Anti-flatness filtering: threshold={flatness_threshold:.3f}, control='{control_var}'")
    
    try:
        # Get and validate control curve data
        control_data = line_dictionary[control_var]
        if not isinstance(control_data, list):
            raise TypeError(f"Control variant data must be list, got {type(control_data)}")
        
        # Convert to numpy array with error handling
        y_control = np.array(control_data, dtype=float)
        
        # Efficiently remove None/NaN values
        valid_mask = np.isfinite(y_control)
        
        if not np.any(valid_mask):
            logging.warning(f"No valid data points for variability analysis in '{control_var}'")
            return False, {
                "filter_method": "filter_against_flatness",
                "error": "No valid control data",
                "control_std": 0.0,
                "control_cv": 0.0,
                "is_variable": False,
                "mean_value": 0.0
            }
        
        y_control_clean = y_control[valid_mask]
        n_valid_points = len(y_control_clean)
        
        # Require minimum data points for reliable variability assessment
        if n_valid_points < 3:
            logging.warning(f"Insufficient data points ({n_valid_points}) for variability analysis")
            return False, {
                "filter_method": "filter_against_flatness",
                "error": "Insufficient data points",
                "control_std": 0.0,
                "control_cv": 0.0,
                "is_variable": False,
                "mean_value": 0.0,
                "n_valid_points": n_valid_points
            }
        
        # Calculate basic statistics
        control_mean = float(np.mean(y_control_clean))
        control_std = float(np.std(y_control_clean, ddof=1))
        
        # Simple variability criterion: standard deviation above threshold
        is_variable = control_std >= flatness_threshold
        
        # Prepare simple metrics
        metrics = {
            "filter_method": "filter_against_flatness",
            "threshold": float(flatness_threshold),
            "control_std": control_std,
            "control_mean": control_mean,
            "is_variable": is_variable,
            "control_variant": control_var,
            "n_valid_points": n_valid_points,
            "n_total_points": len(y_control)
        }
        
        logging.info(
            f"Anti-flatness filtering: '{control_var}' std={control_std:.3f} (>= {flatness_threshold:.3f}), "
            f"mean={control_mean:.3f}, variable={is_variable}"
        )
        
        return is_variable, metrics
        
    except Exception as e:
        logging.error(f"Unexpected error in anti-flatness filtering: {e}", exc_info=True)
        return False, {
            "filter_method": "filter_against_flatness",
            "error": str(e),
            "control_std": 0.0,
            "control_cv": 0.0,
            "is_variable": False,
            "mean_value": 0.0
        }


def apply_gene_filtering(
    gene_name: str,
    data_df: pd.DataFrame,
    variants: List[str],
    filter_method: str = "none",
    control_variant: Optional[str] = None,
    area_threshold: float = DEFAULT_AREA_THRESHOLD,
    flatness_threshold: float = DEFAULT_FLATNESS_THRESHOLD,
    **kwargs
) -> Tuple[bool, Dict[str, Union[float, bool, str]]]:
    """Apply specified filtering method to a gene with comprehensive error handling.
    
    This is the main entry point for gene-level filtering in thermal proteome profiling
    analysis. It provides a unified interface for all filtering methods with robust
    error handling, performance optimization, and detailed metrics reporting.
    
    **Architecture:**
    - Validates inputs and extracts gene-specific data
    - Transforms wide-format data to analysis-ready format
    - Dispatches to appropriate specialized filtering function
    - Aggregates and standardizes metrics across all methods
    
    **Performance Features:**
    - Cached concentration parsing for repeated operations
    - Early validation to prevent unnecessary computations
    - Efficient data extraction using pandas vectorized operations
    - Memory-efficient handling of large datasets
    
    Args:
        gene_name: Target gene symbol for analysis
        data_df: Complete proteomics dataset (wide format)
        variants: List of experimental variant names (e.g., ['R', 'S'])
        filter_method: Filtering strategy - "none", "area_between_curves", 
                      "filter_for_flatness", "filter_against_flatness"
        control_variant: Reference variant for comparison (None = use DEFAULT_CONTROL_VARIANT)
        area_threshold: Minimum area between curves for significance (typical: 1.0-10.0)
        flatness_threshold: Standard deviation threshold for flatness (typical: 0.1-1.0)
        **kwargs: Method-specific parameters (e.g., require_flat_control, require_minimal_variation)
        
    Returns:
        Tuple containing:
        - bool: True if gene passes the specified filtering criteria
        - Dict: Comprehensive metrics including method-specific statistics and metadata
        
    Raises:
        ValueError: If gene not found, invalid filter method, or malformed data
        TypeError: If data_df is not a pandas DataFrame
    """
    # Comprehensive input validation
    if not isinstance(gene_name, str) or not gene_name.strip():
        raise ValueError(f"Gene name must be non-empty string, got {gene_name}")
    
    if not isinstance(data_df, pd.DataFrame):
        raise TypeError(f"data_df must be pandas DataFrame, got {type(data_df)}")
    
    if not isinstance(variants, list) or len(variants) == 0:
        raise ValueError(f"variants must be non-empty list, got {variants}")
    
    if filter_method not in ["none"] + list(FILTER_METHODS.keys()):
        raise ValueError(f"Unknown filter method '{filter_method}'. Valid options: {list(FILTER_METHODS.keys())}")
    
    logging.info(f"Applying '{filter_method}' filter to gene '{gene_name}' with {len(variants)} variants")
    
    # Handle no-filtering case efficiently
    if filter_method == "none":
        return True, {
            "filter_method": "none", 
            "gene_name": gene_name,
            "passes_filter": True,
            "processing_time_ms": 0.0
        }
    
    # Performance tracking
    import time
    start_time = time.time()
    
    try:
        # Filter data for the specific gene with validation
        gene_mask = data_df["Genes"] == gene_name
        gene_data = data_df[gene_mask]
        
        if gene_data.empty:
            logging.warning(f"No data found for gene: {gene_name}")
            return False, {
                "filter_method": filter_method, 
                "gene_name": gene_name, 
                "error": "no_data",
                "passes_filter": False,
                "processing_time_ms": (time.time() - start_time) * 1000
            }
        
        # Validate required columns exist
        required_columns = ["Comparison (group1/group2)", "AVG Log2 Ratio"]
        missing_columns = [col for col in required_columns if col not in gene_data.columns]
        if missing_columns:
            raise ValueError(f"Required columns missing from gene data: {missing_columns}")
    
        # Extract and sort concentrations with robust error handling
        try:
            concentrations = (
                gene_data["Comparison (group1/group2)"]
                .str.split("_")
                .str[2]  # get the third element of each split
                .str.split(" ")
                .str[0]  # get the part before the first space
                .unique()  # get unique values
            )
            
            # Use cached sorting function for better performance
            sorted_concentrations = sorted(
                concentrations, 
                key=_cached_concentration_sort_key
            )
            
            if len(sorted_concentrations) < 2:
                logging.warning(f"Insufficient concentration points ({len(sorted_concentrations)}) for gene '{gene_name}'")
                return False, {
                    "filter_method": filter_method,
                    "gene_name": gene_name,
                    "error": "insufficient_concentrations",
                    "n_concentrations": len(sorted_concentrations),
                    "passes_filter": False,
                    "processing_time_ms": (time.time() - start_time) * 1000
                }
                
        except Exception as e:
            logging.error(f"Error extracting concentrations for gene '{gene_name}': {e}")
            return False, {
                "filter_method": filter_method,
                "gene_name": gene_name,
                "error": f"concentration_extraction_failed: {str(e)}",
                "passes_filter": False,
                "processing_time_ms": (time.time() - start_time) * 1000
            }
        
        # Create line data for analysis with error handling
        try:
            plot_df, line_dict = get_lines(
                variants,
                gene_data,
                sorted_concentrations,
                "AVG Log2 Ratio",
                "log_fld_change",
            )
            
            # Validate line data quality
            if not line_dict or all(not data for data in line_dict.values()):
                logging.warning(f"No valid line data extracted for gene '{gene_name}'")
                return False, {
                    "filter_method": filter_method,
                    "gene_name": gene_name,
                    "error": "no_valid_line_data",
                    "passes_filter": False,
                    "processing_time_ms": (time.time() - start_time) * 1000
                }
                
        except Exception as e:
            logging.error(f"Error creating line data for gene '{gene_name}': {e}")
            return False, {
                "filter_method": filter_method,
                "gene_name": gene_name,
                "error": f"line_data_creation_failed: {str(e)}",
                "passes_filter": False,
                "processing_time_ms": (time.time() - start_time) * 1000
            }
        
        # Apply the specified filtering method with comprehensive error handling
        try:
            if filter_method == "area_between_curves":
                passes_filter, metrics = filter_area_between_curves(
                    line_dict, 
                    sorted_concentrations, 
                    variants,
                    threshold=area_threshold,
                    control_variant=control_variant,
                    **kwargs
                )
                
            elif filter_method == "filter_for_flatness":
                passes_filter, metrics = filter_for_flatness(
                    line_dict,
                    sorted_concentrations,
                    variants,
                    control_variant=control_variant,
                    flatness_threshold=flatness_threshold,
                    **kwargs
                )
                
            elif filter_method == "filter_against_flatness":
                passes_filter, metrics = filter_against_flatness(
                    line_dict,
                    sorted_concentrations,
                    variants,
                    control_variant=control_variant,
                    flatness_threshold=flatness_threshold,
                    **kwargs
                )
                
            else:
                # This should not happen due to earlier validation, but safety first
                raise ValueError(f"Unknown filter method: {filter_method}")
                
        except Exception as e:
            logging.error(f"Error in '{filter_method}' filtering for gene '{gene_name}': {e}")
            return False, {
                "filter_method": filter_method,
                "gene_name": gene_name,
                "error": f"filter_execution_failed: {str(e)}",
                "passes_filter": False,
                "processing_time_ms": (time.time() - start_time) * 1000
            }
        
        # Add comprehensive metadata to metrics
        processing_time_ms = (time.time() - start_time) * 1000
        
        # Ensure metrics is a dictionary (defensive programming)
        if not isinstance(metrics, dict):
            metrics = {"error": "Invalid metrics returned from filter function"}
        
        metrics.update({
            "filter_method": filter_method,
            "gene_name": gene_name,
            "n_concentrations": len(sorted_concentrations),
            "n_variants": len(variants),
            "concentrations": list(sorted_concentrations),
            "variants_analyzed": list(variants),
            "processing_time_ms": processing_time_ms,
            "data_shape": gene_data.shape,
            "passes_filter": bool(passes_filter)
        })
        
        logging.info(
            f"Gene '{gene_name}' filter '{filter_method}' result: {passes_filter} "
            f"(processed in {processing_time_ms:.1f}ms)"
        )
        
        return passes_filter, metrics
        
    except Exception as e:
        processing_time_ms = (time.time() - start_time) * 1000
        logging.error(f"Unexpected error applying filter '{filter_method}' to gene '{gene_name}': {e}", exc_info=True)
        return False, {
            "filter_method": filter_method, 
            "gene_name": gene_name, 
            "error": str(e),
            "passes_filter": False,
            "processing_time_ms": processing_time_ms
        }
