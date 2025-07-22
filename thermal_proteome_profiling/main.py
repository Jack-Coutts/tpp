import sys
import logging
from pathlib import Path
from typing import Tuple, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from PySide6.QtWidgets import (QApplication, QWidget, QLabel, QPushButton, QLineEdit, QCheckBox,
                               QFileDialog, QRadioButton, QButtonGroup, QVBoxLayout, QHBoxLayout, QTextEdit, QMessageBox)
from PySide6.QtCore import Qt

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
    if not path.exists():
        raise ValueError(
            f"Error: The path '{filepath_str}' does not exist."
        )
    if not path.is_file():
        raise ValueError(
            f"Error: The path '{filepath_str}' is not a file."
        )
    return path

# custom type converter for argparse using pathlib
# checks if the given path exists and is a directory. Returns a Path object.
def check_directory(dir_path_str: str) -> Path:
    path = Path(dir_path_str)
    if not path.exists():
        raise ValueError(
            f"Error: The path '{dir_path_str}' does not exist."
        )
    if not path.is_dir():
        raise ValueError(
            f"Error: The path '{dir_path_str}' is not a directory."
        )
    return path

def read_gene_list(file_path: Path) -> List:
    """Read gene names from a text file, one per line"""
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

def plotter(
    df, y_min, y_max, cell_line, gene_name, output_dir, variants, error_df=None
):
    try:
        # Set seaborn style
        sns.set_style("whitegrid")
        plt.rcParams["font.size"] = 11

        # ensure line colour consistency
        try:
            sorted_variants = sorted(variants, key=lambda x: x[0])
            variant_palette = {
                sorted_variants[0]: "#1f77b4",
                sorted_variants[1]: "#d62728",
            }
            logging.info(f"Sorted variants exist: {sorted_variants}")
        except Exception as e:
            logging.error(
                f"Less than 2 protein varainats exist: Variants: {variants}\n Error: {e}"
            )
            return

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
        logging.info(f"{gene_name} from {cell_line} plotted")

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

def filter_or_plot_gene(
    gene_name,
    data_df,
    output_dir,
    highest_y_value,
    lowest_y_value,
    filter=False,
    filter_flat=False,
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

    # Create line df
    plot_df, lines = get_lines(
        variants, gene_data, sorted_concentrations, "AVG Log2 Ratio", "log_fld_change"
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
            error_df,
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

class TPPPlotterGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TPP Plotter")
        self.setMinimumSize(800, 600)
        self.init_ui()

    def init_ui(self):
        # Widgets
        self.data_label = QLabel("Data File:")
        self.data_path = QLineEdit()
        self.data_browse = QPushButton("Browse")

        self.output_label = QLabel("Output Folder:")
        self.output_path = QLineEdit()
        self.output_browse = QPushButton("Browse")

        self.error_bars_checkbox = QCheckBox("Add Error Bars")

        # Radio buttons for modes
        self.mode_group = QButtonGroup(self)
        self.mode_gene = QRadioButton("Plot Single Gene")
        self.mode_gene_list = QRadioButton("Plot from Gene List")
        self.mode_filter = QRadioButton("Filter Proteins")
        self.mode_all = QRadioButton("Plot All Proteins")
        self.mode_gene.setChecked(True)

        self.mode_group.addButton(self.mode_gene)
        self.mode_group.addButton(self.mode_gene_list)
        self.mode_group.addButton(self.mode_filter)
        self.mode_group.addButton(self.mode_all)

        self.gene_label = QLabel("Gene Name:")
        self.gene_input = QLineEdit()

        self.gene_list_label = QLabel("Gene List File:")
        self.gene_list_path = QLineEdit()
        self.gene_list_browse = QPushButton("Browse")

        self.run_button = QPushButton("Run")

        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setFixedHeight(150)

        # Layouts
        main_layout = QVBoxLayout()

        # Data file layout
        data_layout = QHBoxLayout()
        data_layout.addWidget(self.data_label)
        data_layout.addWidget(self.data_path)
        data_layout.addWidget(self.data_browse)
        main_layout.addLayout(data_layout)

        # Output folder layout
        output_layout = QHBoxLayout()
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(self.output_path)
        output_layout.addWidget(self.output_browse)
        main_layout.addLayout(output_layout)

        # Error bars
        main_layout.addWidget(self.error_bars_checkbox)

        # Mode selection
        mode_layout = QVBoxLayout()
        mode_layout.addWidget(self.mode_gene)
        mode_layout.addWidget(self.mode_gene_list)
        mode_layout.addWidget(self.mode_filter)
        mode_layout.addWidget(self.mode_all)
        main_layout.addLayout(mode_layout)

        # Gene input
        gene_layout = QHBoxLayout()
        gene_layout.addWidget(self.gene_label)
        gene_layout.addWidget(self.gene_input)
        main_layout.addLayout(gene_layout)

        # Gene list file layout
        gene_list_layout = QHBoxLayout()
        gene_list_layout.addWidget(self.gene_list_label)
        gene_list_layout.addWidget(self.gene_list_path)
        gene_list_layout.addWidget(self.gene_list_browse)
        main_layout.addLayout(gene_list_layout)

        # Run button
        main_layout.addWidget(self.run_button)

        # Status box
        main_layout.addWidget(self.status_box)

        self.setLayout(main_layout)

        # Connections
        self.data_browse.clicked.connect(self.browse_data_file)
        self.output_browse.clicked.connect(self.browse_output_folder)
        self.gene_list_browse.clicked.connect(self.browse_gene_list_file)
        self.run_button.clicked.connect(self.on_run)
        self.mode_gene.toggled.connect(self.update_mode_inputs)
        self.mode_gene_list.toggled.connect(self.update_mode_inputs)
        self.mode_filter.toggled.connect(self.update_mode_inputs)
        self.mode_all.toggled.connect(self.update_mode_inputs)

        self.update_mode_inputs()

    def log(self, message):
        self.status_box.append(message)

    def browse_data_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Data File", "", "CSV or TSV Files (*.csv *.tsv *.txt);;All Files (*.*)")
        if file_path:
            self.data_path.setText(file_path)

    def browse_output_folder(self):
        folder_path = QFileDialog.getExistingDirectory(self, "Select Output Folder")
        if folder_path:
            self.output_path.setText(folder_path)

    def browse_gene_list_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Gene List File", "", "Text Files (*.txt);;All Files (*.*)")
        if file_path:
            self.gene_list_path.setText(file_path)

    def update_mode_inputs(self):
        # Show/hide input widgets based on selected mode
        self.gene_label.setEnabled(self.mode_gene.isChecked())
        self.gene_input.setEnabled(self.mode_gene.isChecked())

        self.gene_list_label.setEnabled(self.mode_gene_list.isChecked())
        self.gene_list_path.setEnabled(self.mode_gene_list.isChecked())
        self.gene_list_browse.setEnabled(self.mode_gene_list.isChecked())

    def on_run(self):
        try:
            data_file = self.data_path.text()
            output_folder = self.output_path.text()
            error_bars = self.error_bars_checkbox.isChecked()

            if not data_file:
                QMessageBox.warning(self, "Input Error", "Please select a data file.")
                return
            if not output_folder:
                QMessageBox.warning(self, "Input Error", "Please select an output folder.")
                return

            data_file = check_filepath(data_file)
            output_folder = check_directory(output_folder)

            mode = None
            if self.mode_gene.isChecked():
                mode = "gene"
            elif self.mode_gene_list.isChecked():
                mode = "gene_list"
            elif self.mode_filter.isChecked():
                mode = "filter"
            elif self.mode_all.isChecked():
                mode = "all"

            if mode == "gene":
                gene_name = self.gene_input.text().strip()
                if not gene_name:
                    QMessageBox.warning(self, "Input Error", "Enter a gene name for single gene mode.")
                    return
            elif mode == "gene_list":
                gene_list_file = self.gene_list_path.text()
                if not gene_list_file:
                    QMessageBox.warning(self, "Input Error", "Select a gene list file.")
                    return
                gene_list_file = check_filepath(gene_list_file)

            self.log(f"Loading data from {data_file}")
            data_df, y_max, y_min = read_data(data_file)

            if mode == "gene":
                self.log(f"Plotting single gene: {gene_name}")
                filter_or_plot_gene(
                    gene_name, data_df, output_folder, y_max, y_min, error_bars=error_bars
                )
            elif mode == "gene_list":
                self.log(f"Reading gene list from {gene_list_file}")
                gene_list = read_gene_list(gene_list_file)
                available_genes = set(data_df["Genes"].unique())
                found_genes = [gene for gene in gene_list if gene in available_genes]
                missing_genes = [gene for gene in gene_list if gene not in available_genes]
                self.log(f"Found {len(found_genes)} genes, {len(missing_genes)} missing")

                filtered_df = data_df[data_df["Genes"].isin(found_genes)]
                y_max = filtered_df["AVG Log2 Ratio"].max()
                y_min = filtered_df["AVG Log2 Ratio"].min()

                for i, gene_name in enumerate(found_genes, 1):
                    self.log(f"Plotting {i}/{len(found_genes)}: {gene_name}")
                    filter_or_plot_gene(
                        gene_name, data_df, output_folder, y_max, y_min, error_bars=error_bars
                    )

            elif mode == "filter":
                genes = data_df["Genes"].unique()
                self.log("Applying filter and adding genes to gene list...")
                for i, gene_name in enumerate(genes, 1):
                    self.log(f"protein {i}/{len(genes)}")
                    filter_or_plot_gene(gene_name, data_df, output_folder, y_max, y_min, filter=True, filter_flat=True)

            elif mode == "all":
                genes = data_df["Genes"].unique()
                self.log("Plotting all genes...")
                for i, gene_name in enumerate(genes, 1):
                    self.log(f"protein {i}/{len(genes)}")
                    filter_or_plot_gene(gene_name, data_df, output_folder, y_max, y_min, error_bars=error_bars)

            self.log("Processing complete.")
            QMessageBox.information(self, "Success", "Processing complete!")

        except Exception as e:
            self.log(f"Error: {str(e)}")
            logging.error(f"Error in GUI run: {e}")
            QMessageBox.critical(self, "Error", str(e))

def main():
    app = QApplication(sys.argv)
    window = TPPPlotterGUI()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
