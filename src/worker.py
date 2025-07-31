import logging

from PySide6.QtCore import QObject, Signal, Slot

from .utils import get_variants, plot_gene, read_data, read_gene_list


class Worker(QObject):
    finished = Signal()
    error = Signal(str)
    stopped = Signal()

    def __init__(
        self,
        data_file,
        output_folder,
        single_mode,
        error_bars,
        mode,
        line_name,
        gene_name=None,
        gene_list_file=None,
        parent=None,
    ):
        super().__init__(parent)
        self.data_file = data_file
        self.output_folder = output_folder
        self.single_mode = single_mode
        self.error_bars = error_bars
        self.mode = mode
        self.gene_name = gene_name
        self.gene_list_file = gene_list_file
        self.line_name = line_name
        self._stop_requested = False  # Flag for interruption

    def request_stop(self):
        self._stop_requested = True

    @Slot()
    def run(self):
        try:
            logging.info(f"Loading data from {self.data_file}")
            data_df, y_max, y_min = read_data(self.data_file)

            # the the compound variants (line names)
            print(self.line_name)
            variants = get_variants(data_df, self.single_mode, self.line_name)

            if self.mode == "gene":
                if self._stop_requested:
                    logging.info("Processing stopped")
                    self.stopped.emit()
                    return
                logging.info(f"Plotting single gene: {self.gene_name}")
                plot_gene(
                    self.gene_name,
                    data_df,
                    variants,
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
                    if self._stop_requested:
                        logging.info("Processing stopped")
                        self.stopped.emit()
                        return
                    logging.info(f"Plotting {i}/{len(found_genes)}: {gene_name}")
                    plot_gene(
                        gene_name,
                        data_df,
                        self.output_folder,
                        variants,
                        y_max,
                        y_min,
                        error_bars=self.error_bars,
                    )

            elif self.mode == "all":
                genes = data_df["Genes"].unique()
                logging.info("Plotting all genes...")
                for i, gene_name in enumerate(genes, 1):
                    if self._stop_requested:
                        logging.info("Processing stopped")
                        self.stopped.emit()
                        return
                    logging.info(f"{i}/{len(genes)} plotting {gene_name}")
                    plot_gene(
                        gene_name,
                        data_df,
                        self.output_folder,
                        variants,
                        y_max,
                        y_min,
                        error_bars=self.error_bars,
                    )

            logging.info("Processing complete.")
            self.finished.emit()
        except Exception as e:
            logging.error(f"Error in worker: {e}")
            self.error.emit(str(e))

        finally:
            # THIS IS THE MOST IMPORTANT PART
            # This will run no matter how the 'try' block exits.
            # It allows the user to run proccessing again as the thread is closed
            logging.info("Worker's run method finished, emitting 'finished' signal.")
            self.finished.emit()
