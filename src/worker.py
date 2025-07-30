import logging

from PySide6.QtCore import QObject, Signal, Slot

from .utils import (
    filter_or_plot_gene,
    read_data,
    read_gene_list,
)


class Worker(QObject):
    finished = Signal()
    error = Signal(str)

    def __init__(
        self,
        data_file,
        output_folder,
        error_bars,
        mode,
        gene_name=None,
        gene_list_file=None,
        parent=None,
    ):
        super().__init__(parent)
        self.data_file = data_file
        self.output_folder = output_folder
        self.error_bars = error_bars
        self.mode = mode
        self.gene_name = gene_name
        self.gene_list_file = gene_list_file

    @Slot()
    def run(self):
        try:
            logging.info(f"Loading data from {self.data_file}")
            data_df, y_max, y_min = read_data(self.data_file)

            if self.mode == "gene":
                logging.info(f"Plotting single gene: {self.gene_name}")
                filter_or_plot_gene(
                    self.gene_name,
                    data_df,
                    self.output_folder,
                    y_max,
                    y_min,
                    error_bars=self.error_bars,
                )
            elif self.mode == "gene_list":
                logging.info(f"Reading gene list from {self.gene_list_file}")
                gene_list = read_gene_list(self.gene_list_file)
                available_genes = set(data_df["Genes"].unique())
                found_genes = [gene for gene in gene_list if gene in available_genes]
                missing_genes = [
                    gene for gene in gene_list if gene not in available_genes
                ]
                logging.info(
                    f"Found {len(found_genes)} genes, {len(missing_genes)} missing"
                )

                filtered_df = data_df[data_df["Genes"].isin(found_genes)]
                y_max = filtered_df["AVG Log2 Ratio"].max()
                y_min = filtered_df["AVG Log2 Ratio"].min()

                for i, gene_name in enumerate(found_genes, 1):
                    logging.info(f"Plotting {i}/{len(found_genes)}: {gene_name}")
                    filter_or_plot_gene(
                        gene_name,
                        data_df,
                        self.output_folder,
                        y_max,
                        y_min,
                        error_bars=self.error_bars,
                    )
            elif self.mode == "filter":
                genes = data_df["Genes"].unique()
                logging.info("Applying filter and adding genes to gene list...")
                for i, gene_name in enumerate(genes, 1):
                    logging.info(f"protein {i}/{len(genes)}")
                    filter_or_plot_gene(
                        gene_name,
                        data_df,
                        self.output_folder,
                        y_max,
                        y_min,
                        filter=True,
                        filter_flat=True,
                    )
            elif self.mode == "all":
                genes = data_df["Genes"].unique()
                logging.info("Plotting all genes...")
                for i, gene_name in enumerate(genes, 1):
                    logging.info(f"protein {i}/{len(genes)}")
                    filter_or_plot_gene(
                        gene_name,
                        data_df,
                        self.output_folder,
                        y_max,
                        y_min,
                        error_bars=self.error_bars,
                    )

            logging.info("Processing complete.")
            self.finished.emit()
        except Exception as e:
            logging.error(f"Error in worker: {e}")
            self.error.emit(str(e))
