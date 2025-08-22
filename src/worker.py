import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from PySide6.QtCore import QObject, Signal, Slot

from .utils import get_variants, plot_gene, read_data, read_gene_list, apply_gene_filtering


class Worker(QObject):
    finished = Signal()
    error = Signal(str)
    stopped = Signal()

    def __init__(
        self,
        data_file: Path,
        output_folder: Path,
        single_mode: bool,
        error_bars: bool,
        mode: str,
        line_name: str,
        gene_name: Optional[str] = None,
        gene_list_file: Optional[Path] = None,
        filter_method: str = "none",
        control_variant: Optional[str] = None,
        area_threshold: float = 2.0,
        flatness_threshold: float = 0.35,
        require_flat_control: bool = False,
        save_filter_metrics: bool = False,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self.data_file = data_file
        self.output_folder = output_folder
        self.single_mode = single_mode
        self.error_bars = error_bars
        self.mode = mode
        self.gene_name = gene_name
        self.gene_list_file = gene_list_file
        self.line_name = line_name
        
        # Filtering parameters
        self.filter_method = filter_method
        self.control_variant = control_variant
        self.area_threshold = area_threshold
        self.flatness_threshold = flatness_threshold
        self.require_flat_control = require_flat_control
        self.save_filter_metrics = save_filter_metrics
        
        self._stop_requested = False  # Flag for interruption

    def _save_filter_metrics_csv(self, metrics_data: Dict[str, Dict], output_file: Path) -> None:
        """Save filtering metrics to CSV file.
        
        Args:
            metrics_data: Dictionary mapping gene names to their metrics
            output_file: Path where CSV file should be saved
        """
        try:
            if not metrics_data:
                logging.warning("No metrics data to save")
                return
            
            # Flatten the metrics for CSV format
            flattened_data = []
            for gene_name, metrics in metrics_data.items():
                # Ensure gene_name is in the metrics
                if 'gene_name' not in metrics:
                    metrics['gene_name'] = gene_name
                
                # Convert complex data types to strings for CSV compatibility
                flattened_metrics = {}
                for key, value in metrics.items():
                    if isinstance(value, (list, dict)):
                        flattened_metrics[key] = str(value)
                    elif isinstance(value, (tuple,)):
                        flattened_metrics[key] = str(list(value))
                    else:
                        flattened_metrics[key] = value
                
                flattened_data.append(flattened_metrics)
            
            # Create DataFrame and save as CSV
            df = pd.DataFrame(flattened_data)
            
            # Reorder columns to put important ones first
            important_cols = ['gene_name', 'filter_method', 'passes_filter', 'processing_time_ms']
            other_cols = [col for col in df.columns if col not in important_cols]
            ordered_cols = [col for col in important_cols if col in df.columns] + sorted(other_cols)
            df = df[ordered_cols]
            
            df.to_csv(output_file, index=False)
            logging.info(f"Filter metrics saved to: {output_file}")
            
        except Exception as e:
            logging.error(f"Error saving filter metrics to CSV: {e}")

    def request_stop(self) -> None:
        self._stop_requested = True

    def _should_process_gene(self, gene_name: str, data_df, variants) -> Tuple[bool, Dict]:
        """Check if gene should be processed based on filtering criteria.
        
        Returns:
            Tuple of (should_process, filter_metrics)
        """
        if self.filter_method == "none":
            return True, {"filter_method": "none"}
        
        try:
            passes_filter, metrics = apply_gene_filtering(
                gene_name=gene_name,
                data_df=data_df,
                variants=variants,
                filter_method=self.filter_method,
                control_variant=self.control_variant,
                area_threshold=self.area_threshold,
                flatness_threshold=self.flatness_threshold,
                require_flat_control=self.require_flat_control
            )
            
            return passes_filter, metrics
            
        except Exception as e:
            logging.error(f"Error filtering gene {gene_name}: {e}")
            return False, {"error": str(e)}

    @Slot()
    def run(self) -> None:
        try:
            logging.info(f"Loading data from {self.data_file}")
            data_df, y_max, y_min = read_data(self.data_file)

            # Get the compound variants (line names)
            variants = get_variants(data_df, self.single_mode, self.line_name)

            if self.mode == "gene":
                if self._stop_requested:
                    logging.info("Processing stopped")
                    self.stopped.emit()
                    return
                
                should_process, filter_metrics = self._should_process_gene(self.gene_name, data_df, variants)
                
                if should_process:
                    logging.info(f"Plotting single gene: {self.gene_name}")
                    plot_gene(
                        self.gene_name,
                        data_df,
                        self.output_folder,
                        variants,
                        y_max,
                        y_min,
                        error_bars=self.error_bars,
                    )
                else:
                    logging.info(f"Gene {self.gene_name} filtered out: {filter_metrics}")
                
                if self.save_filter_metrics and filter_metrics:
                    metrics_file = self.output_folder / f"{self.gene_name}_filter_metrics.csv"
                    self._save_filter_metrics_csv({self.gene_name: filter_metrics}, metrics_file)
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

                processed_count = 0
                all_filter_metrics = {}
                
                for i, gene_name in enumerate(found_genes, 1):
                    if self._stop_requested:
                        logging.info("Processing stopped")
                        self.stopped.emit()
                        return
                    
                    should_process, filter_metrics = self._should_process_gene(gene_name, data_df, variants)
                    all_filter_metrics[gene_name] = filter_metrics
                    
                    if should_process:
                        processed_count += 1
                        logging.info(f"Plotting {processed_count}/{len(found_genes)}: {gene_name}")
                        plot_gene(
                            gene_name,
                            data_df,
                            self.output_folder,
                            variants,
                            y_max,
                            y_min,
                            error_bars=self.error_bars,
                        )
                    else:
                        logging.info(f"Gene {gene_name} filtered out: {filter_metrics}")
                
                if self.save_filter_metrics:
                    metrics_file = self.output_folder / "gene_list_filter_metrics.csv"
                    self._save_filter_metrics_csv(all_filter_metrics, metrics_file)
                
                logging.info(f"Processed {processed_count} out of {len(found_genes)} genes after filtering")

            elif self.mode == "all":
                genes = data_df["Genes"].unique()
                logging.info("Plotting all genes...")
                processed_count = 0
                all_filter_metrics = {}
                
                for i, gene_name in enumerate(genes, 1):
                    if self._stop_requested:
                        logging.info("Processing stopped")
                        self.stopped.emit()
                        return
                    
                    should_process, filter_metrics = self._should_process_gene(gene_name, data_df, variants)
                    all_filter_metrics[gene_name] = filter_metrics
                    
                    if should_process:
                        processed_count += 1
                        logging.info(f"{processed_count}/{len(genes)} plotting {gene_name}")
                        plot_gene(
                            gene_name,
                            data_df,
                            self.output_folder,
                            variants,
                            y_max,
                            y_min,
                            error_bars=self.error_bars,
                        )
                    else:
                        logging.info(f"Gene {gene_name} filtered out: {filter_metrics}")
                
                if self.save_filter_metrics:
                    metrics_file = self.output_folder / "all_genes_filter_metrics.csv"
                    self._save_filter_metrics_csv(all_filter_metrics, metrics_file)
                
                logging.info(f"Processed {processed_count} out of {len(genes)} genes after filtering")

            logging.info("Processing complete.")
        except FileNotFoundError as e:
            error_msg = f"File not found: {e}"
            logging.error(error_msg)
            self.error.emit(error_msg)
        except PermissionError as e:
            error_msg = f"Permission denied: {e}"
            logging.error(error_msg)
            self.error.emit(error_msg)
        except ValueError as e:
            error_msg = f"Data validation error: {e}"
            logging.error(error_msg)
            self.error.emit(error_msg)
        except Exception as e:
            error_msg = f"Unexpected error in worker: {e}"
            logging.error(error_msg)
            self.error.emit(error_msg)
        finally:
            # Ensure finished signal is always emitted for proper cleanup
            logging.info("Worker's run method finished, emitting 'finished' signal.")
            self.finished.emit()
