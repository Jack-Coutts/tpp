import logging
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QThread, Qt, Slot
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSpinBox,
    QDoubleSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .logging_handler import QtLogHandler
from .utils import check_directory, check_filepath, FILTER_METHODS, DEFAULT_CONTROL_VARIANT  # Import utilities
from .worker import Worker


class TPPPlotterGUI(QWidget):
    """Main GUI window for Thermal Proteome Profiling data analysis and visualization.
    
    This application provides a user-friendly interface for:
    - Loading and processing proteomics data
    - Configuring filtering parameters for data quality control
    - Generating publication-quality thermal shift plots
    - Managing multi-threaded data processing workflows
    """
    
    def __init__(self) -> None:
        """Initialize the main GUI window with default settings and logging."""
        super().__init__()
        self.worker = None
        self.thread = None
        self.setWindowTitle("TPP Plotter - Thermal Proteome Profiling Analysis")
        self.setMinimumSize(900, 700)  # Increased size for better spacing
        self.init_ui()
        self.setup_logging()

    def _create_file_widgets(self) -> None:
        """Create file selection widgets."""
        self.data_label = QLabel("Data File:")
        self.data_path = QLineEdit()
        self.data_browse = QPushButton("Browse")

        self.output_label = QLabel("Output Folder:")
        self.output_path = QLineEdit()
        self.output_browse = QPushButton("Browse")

        self.gene_list_label = QLabel("Gene List File:")
        self.gene_list_path = QLineEdit()
        self.gene_list_browse = QPushButton("Browse")

    def _create_mode_widgets(self) -> None:
        """Create mode selection widgets."""
        # Radio buttons for multi-line plotting
        self.multi_line = QButtonGroup(self)
        self.single_line = QRadioButton("Plot One Line")
        self.two_lines = QRadioButton("Plot Two Lines")
        self.single_line.setChecked(True)
        
        self.multi_line.addButton(self.single_line)
        self.multi_line.addButton(self.two_lines)

        self.line_name_label = QLabel("Compound Name:")
        self.name_of_line = QLineEdit()

        # Radio buttons for plotting modes
        self.mode_group = QButtonGroup(self)
        self.mode_gene = QRadioButton("Plot Single Gene")
        self.mode_gene_list = QRadioButton("Plot from Gene List")
        self.mode_all = QRadioButton("Plot All Proteins")
        self.mode_gene.setChecked(True)

        self.mode_group.addButton(self.mode_gene)
        self.mode_group.addButton(self.mode_gene_list)
        self.mode_group.addButton(self.mode_all)

        self.gene_label = QLabel("Gene Name:")
        self.gene_input = QLineEdit()

    def _create_control_widgets(self) -> None:
        """Create control widgets."""
        self.error_bars_checkbox = QCheckBox("Add Error Bars")
        self.run_button = QPushButton("Run")
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        
        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setFixedHeight(150)

    def _create_filtering_widgets(self) -> None:
        """Create comprehensive filtering configuration widgets with improved organization.
        
        Creates a well-organized group of filtering controls including:
        - Filter method selection with tooltips
        - Dynamic parameter controls that adapt to selected method
        - Clear labeling and logical grouping of related options
        """
        # Main filtering group box with better styling
        self.filtering_group = QGroupBox("Data Quality Filtering Options")
        self.filtering_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 2px solid gray;
                border-radius: 8px;
                margin-top: 1ex;
                padding: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
        """)
        
        # Filter method selection with tooltip
        self.filter_method_label = QLabel("Filter Method:")
        self.filter_method_label.setStyleSheet("font-weight: bold;")
        self.filter_method_combo = QComboBox()
        self.filter_method_combo.addItem("No Filtering", "none")
        for method_key, method_desc in FILTER_METHODS.items():
            self.filter_method_combo.addItem(method_desc, method_key)
        self.filter_method_combo.setToolTip(
            "Select filtering method to apply quality control criteria to protein data"
        )
        
        # Control variant selection with better labeling
        self.control_variant_label = QLabel("Reference Variant:")
        self.control_variant_label.setStyleSheet("font-weight: bold;")
        self.control_variant_combo = QComboBox()
        self.control_variant_combo.addItems(["R", "S", "Auto-detect"])
        self.control_variant_combo.setCurrentText("R")  # Default
        self.control_variant_combo.setToolTip(
            "Select which variant to use as the reference/control for filtering calculations"
        )
        
        # Area threshold (for area between curves) with enhanced controls
        self.area_threshold_label = QLabel("Area Threshold:")
        self.area_threshold_spinbox = QDoubleSpinBox()
        self.area_threshold_spinbox.setRange(0.1, 100.0)
        self.area_threshold_spinbox.setValue(2.0)  # Default
        self.area_threshold_spinbox.setSingleStep(0.1)
        self.area_threshold_spinbox.setDecimals(1)
        self.area_threshold_spinbox.setToolTip(
            "Minimum area between curves required to pass filtering (higher = more stringent)"
        )
        
        # Flatness threshold (for flatness filters) with enhanced controls
        self.flatness_threshold_label = QLabel("Flatness Threshold:")
        self.flatness_threshold_spinbox = QDoubleSpinBox()
        self.flatness_threshold_spinbox.setRange(0.01, 2.0)
        self.flatness_threshold_spinbox.setValue(0.35)  # Default
        self.flatness_threshold_spinbox.setSingleStep(0.05)
        self.flatness_threshold_spinbox.setDecimals(2)
        self.flatness_threshold_spinbox.setToolTip(
            "Standard deviation threshold for determining curve flatness (lower = more strict)"
        )
        
        # Advanced options with better grouping
        self.require_flat_control_checkbox = QCheckBox("Require Flat Control (Area Filter)")
        self.require_flat_control_checkbox.setToolTip(
            "Additionally require the control variant to have a flat curve when using area filtering"
        )
        
        self.save_filter_metrics_checkbox = QCheckBox("Save Detailed Filter Metrics")
        self.save_filter_metrics_checkbox.setToolTip(
            "Export detailed filtering statistics and metrics to CSV files for analysis"
        )

    def _create_layouts(self) -> QVBoxLayout:
        """Create and organize layouts with improved spacing and professional organization.
        
        Returns:
            QVBoxLayout: Main layout with properly organized sections
        """
        # Create scroll area for better usability on smaller screens
        scroll_widget = QWidget()
        main_layout = QVBoxLayout(scroll_widget)
        main_layout.setSpacing(20)  # Increased spacing between sections
        
        # === FILE SELECTION SECTION ===
        file_group = QGroupBox("Data Files")
        file_group.setStyleSheet("QGroupBox { font-weight: bold; border: 2px solid gray; border-radius: 5px; margin-top: 1ex; padding: 8px; }")
        file_layout = QVBoxLayout()
        
        # Data file selection
        data_layout = QHBoxLayout()
        data_layout.addWidget(self.data_label)
        data_layout.addWidget(self.data_path)
        data_layout.addWidget(self.data_browse)
        file_layout.addLayout(data_layout)
        
        # Output folder selection
        output_layout = QHBoxLayout()
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(self.output_path)
        output_layout.addWidget(self.output_browse)
        file_layout.addLayout(output_layout)
        
        # Gene list file selection
        gene_list_layout = QHBoxLayout()
        gene_list_layout.addWidget(self.gene_list_label)
        gene_list_layout.addWidget(self.gene_list_path)
        gene_list_layout.addWidget(self.gene_list_browse)
        file_layout.addLayout(gene_list_layout)
        
        file_group.setLayout(file_layout)
        main_layout.addWidget(file_group)

        # === PLOTTING CONFIGURATION SECTION ===
        plot_config_group = QGroupBox("Plot Configuration")
        plot_config_group.setStyleSheet("QGroupBox { font-weight: bold; border: 2px solid gray; border-radius: 5px; margin-top: 1ex; padding: 8px; }")
        plot_config_layout = QVBoxLayout()
        
        # Multi-line selection
        line_selection_layout = QHBoxLayout()
        line_selection_layout.addWidget(QLabel("Plot Mode:"))
        line_selection_layout.addWidget(self.single_line)
        line_selection_layout.addWidget(self.two_lines)
        line_selection_layout.addStretch()  # Add stretch to prevent spreading
        plot_config_layout.addLayout(line_selection_layout)

        # Compound name
        line_name_layout = QHBoxLayout()
        line_name_layout.addWidget(self.line_name_label)
        line_name_layout.addWidget(self.name_of_line)
        plot_config_layout.addLayout(line_name_layout)
        
        plot_config_group.setLayout(plot_config_layout)
        main_layout.addWidget(plot_config_group)

        # === ANALYSIS MODE SECTION ===
        analysis_group = QGroupBox("Analysis Mode")
        analysis_group.setStyleSheet("QGroupBox { font-weight: bold; border: 2px solid gray; border-radius: 5px; margin-top: 1ex; padding: 8px; }")
        analysis_layout = QVBoxLayout()
        
        # Mode selection radio buttons
        mode_layout = QHBoxLayout()
        mode_layout.addWidget(self.mode_gene)
        mode_layout.addWidget(self.mode_gene_list)
        mode_layout.addWidget(self.mode_all)
        mode_layout.addStretch()
        analysis_layout.addLayout(mode_layout)

        # Gene input
        gene_layout = QHBoxLayout()
        gene_layout.addWidget(self.gene_label)
        gene_layout.addWidget(self.gene_input)
        analysis_layout.addLayout(gene_layout)
        
        analysis_group.setLayout(analysis_layout)
        main_layout.addWidget(analysis_group)

        # === FILTERING SECTION ===
        filtering_group = self._create_filtering_layout()
        main_layout.addWidget(filtering_group)

        # === EXECUTION CONTROLS SECTION ===
        controls_group = QGroupBox("Execution Controls")
        controls_group.setStyleSheet("QGroupBox { font-weight: bold; border: 2px solid gray; border-radius: 5px; margin-top: 1ex; padding: 8px; }")
        controls_layout = QVBoxLayout()
        
        # Error bars option
        controls_layout.addWidget(self.error_bars_checkbox)
        
        # Control buttons
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.stop_button)
        button_layout.addStretch()
        controls_layout.addLayout(button_layout)
        
        # Status output
        controls_layout.addWidget(QLabel("Processing Status:"))
        controls_layout.addWidget(self.status_box)
        
        controls_group.setLayout(controls_layout)
        main_layout.addWidget(controls_group)
        
        # Create scroll area and set the widget
        scroll_area = QScrollArea()
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        
        # Wrap scroll area in main layout
        wrapper_layout = QVBoxLayout()
        wrapper_layout.addWidget(scroll_area)
        
        return wrapper_layout

    def _create_filtering_layout(self) -> QGroupBox:
        """Create well-organized filtering section layout with improved spacing and grouping.
        
        Returns:
            QGroupBox containing all filtering controls with professional layout
        """
        # Use grid layout for better organization and alignment
        filtering_layout = QGridLayout()
        filtering_layout.setSpacing(15)  # Increased spacing for better readability
        
        # Row 0: Filter method selection (spans full width)
        filtering_layout.addWidget(self.filter_method_label, 0, 0)
        filtering_layout.addWidget(self.filter_method_combo, 0, 1, 1, 3)  # Span 3 columns
        
        # Row 1: Control variant selection
        filtering_layout.addWidget(self.control_variant_label, 1, 0)
        filtering_layout.addWidget(self.control_variant_combo, 1, 1)
        
        # Row 2: Threshold parameters organized in columns
        filtering_layout.addWidget(self.area_threshold_label, 2, 0)
        filtering_layout.addWidget(self.area_threshold_spinbox, 2, 1)
        filtering_layout.addWidget(self.flatness_threshold_label, 2, 2)
        filtering_layout.addWidget(self.flatness_threshold_spinbox, 2, 3)
        
        # Row 3 & 4: Advanced options
        filtering_layout.addWidget(self.require_flat_control_checkbox, 3, 0, 1, 4)  # Span full width
        filtering_layout.addWidget(self.save_filter_metrics_checkbox, 4, 0, 1, 4)  # Span full width
        
        # Set column stretch ratios for better proportions
        filtering_layout.setColumnStretch(1, 1)
        filtering_layout.setColumnStretch(3, 1)
        
        # Apply layout to the group box
        self.filtering_group.setLayout(filtering_layout)
        
        return self.filtering_group

    def _create_divider(self) -> QFrame:
        """Create a horizontal divider line."""
        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)
        return divider

    def _connect_signals(self) -> None:
        """Connect widget signals to their handlers."""
        self.data_browse.clicked.connect(self.browse_data_file)
        self.output_browse.clicked.connect(self.browse_output_folder)
        self.gene_list_browse.clicked.connect(self.browse_gene_list_file)
        self.run_button.clicked.connect(self.on_run)
        self.stop_button.clicked.connect(self.on_stop)
        
        self.mode_gene.toggled.connect(self.update_mode_inputs)
        self.mode_gene_list.toggled.connect(self.update_mode_inputs)
        self.mode_all.toggled.connect(self.update_mode_inputs)
        
        self.single_line.toggled.connect(self.update_line_num_inputs)
        self.two_lines.toggled.connect(self.update_line_num_inputs)
        
        # Filtering signal connections with error handling
        self.filter_method_combo.currentTextChanged.connect(self.update_filtering_inputs)
        
        # Add validation signals for parameter inputs
        self.area_threshold_spinbox.valueChanged.connect(self._validate_area_threshold)
        self.flatness_threshold_spinbox.valueChanged.connect(self._validate_flatness_threshold)

    def init_ui(self) -> None:
        """Initialize the user interface with modular components."""
        self._create_file_widgets()
        self._create_mode_widgets()
        self._create_control_widgets()
        self._create_filtering_widgets()
        
        main_layout = self._create_layouts()
        self.setLayout(main_layout)
        
        self._connect_signals()
        self.update_mode_inputs()
        self.update_filtering_inputs()
        
        # Apply initial styling improvements
        self._apply_ui_styling()

    @Slot(str)
    def _append_log_message(self, message: str) -> None:
        self.status_box.append(message)

    def setup_logging(self) -> None:
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.INFO)
        # Console handler
        console_handler = logging.StreamHandler()
        console_formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)
        # GUI handler
        self.gui_handler = QtLogHandler()
        self.gui_handler.new_log_message.connect(self._append_log_message)
        root_logger.addHandler(self.gui_handler)

    def browse_data_file(self) -> None:
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Data File",
            "",
            "CSV or TSV Files (*.csv *.tsv *.txt);;All Files (*.*)",
        )
        if file_path:
            self.data_path.setText(file_path)

    def browse_output_folder(self) -> None:
        folder_path = QFileDialog.getExistingDirectory(self, "Select Output Folder")
        if folder_path:
            self.output_path.setText(folder_path)

    def browse_gene_list_file(self) -> None:
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Gene List File", "", "Text Files (*.txt);;All Files (*.*)"
        )
        if file_path:
            self.gene_list_path.setText(file_path)

    def update_mode_inputs(self) -> None:
        # Show/hide input widgets based on selected mode
        self.gene_label.setEnabled(self.mode_gene.isChecked())
        self.gene_input.setEnabled(self.mode_gene.isChecked())

        self.gene_list_label.setEnabled(self.mode_gene_list.isChecked())
        self.gene_list_path.setEnabled(self.mode_gene_list.isChecked())
        self.gene_list_browse.setEnabled(self.mode_gene_list.isChecked())

    def update_line_num_inputs(self) -> None:
        # Show/hide input widgets based on selected line num
        self.line_name_label.setEnabled(self.single_line.isChecked())
        self.name_of_line.setEnabled(self.single_line.isChecked())

    def update_filtering_inputs(self) -> None:
        """Update filtering input availability and labels based on selected method.
        
        Dynamically adjusts the UI to show only relevant parameters for the chosen
        filtering method, improving user experience and reducing confusion.
        
        The method implements context-sensitive parameter visibility:
        - No filtering: All parameters disabled
        - Area between curves: Area threshold + optional flatness controls
        - Flatness filtering: Only flatness threshold + control variant
        """
        try:
            current_method = self.filter_method_combo.currentData()
            
            # Enable/disable controls based on filter method with error checking
            if current_method == "none":
                # Disable all filtering options
                self._set_filtering_controls_enabled(False, False, False, False)
                self._update_filtering_labels("Area Threshold:", "Flatness Threshold:")
                
            elif current_method == "area_between_curves":
                # Enable area filtering options
                self._set_filtering_controls_enabled(True, True, True, True)
                self._update_filtering_labels(
                    "Area Threshold:", 
                    "Flatness Threshold (optional):"
                )
                
            elif current_method in ["filter_for_flatness", "filter_against_flatness"]:
                # Enable flatness filtering options only
                self._set_filtering_controls_enabled(True, False, True, False)
                
                # Update labels based on specific flatness filter type
                if current_method == "filter_for_flatness":
                    self._update_filtering_labels(
                        "Area Threshold:", 
                        "Max Std Dev for Flatness:"
                    )
                else:  # filter_against_flatness
                    self._update_filtering_labels(
                        "Area Threshold:", 
                        "Min Std Dev for Variability:"
                    )
            else:
                # Unknown method - log warning and use safe defaults
                logging.warning(f"Unknown filter method: {current_method}")
                self._set_filtering_controls_enabled(False, False, False, False)
                self._update_filtering_labels("Area Threshold:", "Flatness Threshold:")
                
        except Exception as e:
            logging.error(f"Error updating filtering inputs: {e}")
            # Fallback to safe state
            self._set_filtering_controls_enabled(False, False, False, False)
    
    def _set_filtering_controls_enabled(self, control_variant: bool, area_threshold: bool, 
                                       flatness_threshold: bool, flat_control: bool) -> None:
        """Helper method to set enabled state of filtering controls.
        
        Args:
            control_variant: Enable control variant selection
            area_threshold: Enable area threshold spinbox
            flatness_threshold: Enable flatness threshold spinbox
            flat_control: Enable require flat control checkbox
        """
        try:
            self.control_variant_combo.setEnabled(control_variant)
            self.area_threshold_spinbox.setEnabled(area_threshold)
            self.flatness_threshold_spinbox.setEnabled(flatness_threshold)
            self.require_flat_control_checkbox.setEnabled(flat_control)
        except AttributeError as e:
            logging.error(f"Error setting control enabled states: {e}")
    
    def _update_filtering_labels(self, area_label: str, flatness_label: str) -> None:
        """Helper method to update filtering parameter labels.
        
        Args:
            area_label: Text for area threshold label
            flatness_label: Text for flatness threshold label
        """
        try:
            self.area_threshold_label.setText(area_label)
            self.flatness_threshold_label.setText(flatness_label)
        except AttributeError as e:
            logging.error(f"Error updating filtering labels: {e}")

    def reset_ui(self) -> None:
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        logging.info("App reset to inactive state.")

    def on_stop(self) -> None:
        if hasattr(self, "worker") and self.worker:
            self.worker.request_stop()
            logging.info("Stop requested... waiting for processing to halt.")
            self.stop_button.setEnabled(
                False
            )  # Disable immediately to prevent multiple clicks

    def on_run(self) -> None:
        if self.thread and self.thread.isRunning():
            logging.warning("A process is already running. Please wait or stop it.")
            return
        try:
            data_file_str = self.data_path.text()
            output_folder_str = self.output_path.text()
            error_bars = self.error_bars_checkbox.isChecked()

            if not data_file_str:
                QMessageBox.warning(self, "Input Error", "Please select a data file.")
                return
            if not output_folder_str:
                QMessageBox.warning(
                    self, "Input Error", "Please select an output folder."
                )
                return

            data_file = check_filepath(data_file_str)
            output_folder = check_directory(output_folder_str)

            # different logic for plotting a single line or two lines
            single_mode = False
            line_name = ""
            if self.single_line.isChecked():
                single_mode = True
                line_name = self.name_of_line.text()

            mode = None
            gene_name = None
            gene_list_file = None

            if self.mode_gene.isChecked():
                mode = "gene"
                gene_name = self.gene_input.text().strip()
                if not gene_name:
                    QMessageBox.warning(
                        self, "Input Error", "Enter a gene name for single gene mode."
                    )
                    return
            elif self.mode_gene_list.isChecked():
                mode = "gene_list"
                gene_list_file_str = self.gene_list_path.text()
                if not gene_list_file_str:
                    QMessageBox.warning(self, "Input Error", "Select a gene list file.")
                    return
                gene_list_file = check_filepath(gene_list_file_str)
            # elif self.mode_filter.isChecked():
            #     mode = "filter"
            elif self.mode_all.isChecked():
                mode = "all"

            # Get filtering parameters
            filter_method = self.filter_method_combo.currentData()
            control_variant_text = self.control_variant_combo.currentText()
            control_variant = None if control_variant_text == "Auto-detect" else control_variant_text
            area_threshold = self.area_threshold_spinbox.value()
            flatness_threshold = self.flatness_threshold_spinbox.value()
            require_flat_control = self.require_flat_control_checkbox.isChecked()
            save_filter_metrics = self.save_filter_metrics_checkbox.isChecked()

            # Create worker and thread
            self.thread = QThread()
            self.worker = Worker(
                data_file,
                output_folder,
                single_mode,
                error_bars,
                mode,
                line_name,
                gene_name,
                gene_list_file,
                filter_method,
                control_variant,
                area_threshold,
                flatness_threshold,
                require_flat_control,
                save_filter_metrics,
            )

            self.worker.moveToThread(self.thread)

            # Connect signals
            self.thread.started.connect(self.worker.run)
            self.worker.finished.connect(self.thread.quit)
            self.worker.finished.connect(self.worker.deleteLater)

            # Connect thread.finished to safe cleanup
            self.thread.finished.connect(self._on_cleanup_finished)

            # Start the thread
            self.thread.start()

            # update ui state
            self.run_button.setEnabled(False)
            self.stop_button.setEnabled(True)

        except Exception as e:
            logging.error(f"Error starting worker: {e}")
            self.reset_ui()
            # Ensure stale references are cleared on setup error
            self.thread = None
            self.worker = None

    # safe cleanup method
    @Slot()
    def _on_cleanup_finished(self) -> None:
        """
        This slot runs only after the thread has completely finished.
        It's the safest place to reset the UI and clear references.
        """
        logging.info("Thread finished, performing final cleanup.")
        if self.thread:  # Check if the object still exists before deleting
            # thread deleted when current tasks are finished
            self.thread.deleteLater()

        # Nullify references to ensure a clean state for the next run
        self.thread = None
        self.worker = None

        # Now it's safe to reset the UI
        self.reset_ui()

    def _validate_area_threshold(self, value: float) -> None:
        """Validate area threshold parameter and provide user feedback.
        
        Args:
            value: New threshold value to validate
        """
        try:
            if value < 0.1:
                logging.warning("Area threshold below recommended minimum (0.1)")
            elif value > 50.0:
                logging.warning("Area threshold above typical range (>50)")
        except Exception as e:
            logging.error(f"Error validating area threshold: {e}")
    
    def _validate_flatness_threshold(self, value: float) -> None:
        """Validate flatness threshold parameter and provide user feedback.
        
        Args:
            value: New threshold value to validate
        """
        try:
            if value < 0.01:
                logging.warning("Flatness threshold very low (<0.01)")
            elif value > 1.0:
                logging.warning("Flatness threshold high (>1.0) - may be too permissive")
        except Exception as e:
            logging.error(f"Error validating flatness threshold: {e}")
    
    def _apply_ui_styling(self) -> None:
        """Apply consistent styling across the UI components.
        
        Loads styles from external QSS file for better maintainability.
        Sets up professional styling including:
        - Button styling with hover effects
        - Consistent color scheme
        - Improved spacing and padding
        - Windows-compatible checkbox indicators
        """
        try:
            import sys
            import os
            
            # Determine if we're running as a bundled executable
            if getattr(sys, 'frozen', False):
                # Running as bundled executable (PyInstaller)
                # sys._MEIPASS is the temp folder where PyInstaller extracts files
                bundle_dir = Path(sys._MEIPASS)
            else:
                # Running in normal Python environment
                bundle_dir = Path(__file__).parent
            
            # Try to load the QSS file
            qss_path = bundle_dir / "styles.qss"
            
            if qss_path.exists():
                with open(qss_path, 'r', encoding='utf-8') as f:
                    stylesheet = f.read()
                    self.setStyleSheet(stylesheet)
                    logging.info(f"Loaded styles from {qss_path}")
            else:
                # Fallback to basic styling if QSS file not found
                logging.warning(f"Style sheet not found at {qss_path}, using basic styling")
                self.setStyleSheet("""
                    QWidget {
                        background-color: #f8f9fa;
                        color: #343a40;
                    }
                    QPushButton {
                        background-color: #007bff;
                        color: white;
                        padding: 8px;
                        border-radius: 4px;
                    }
                """)
            
        except Exception as e:
            logging.error(f"Error applying UI styling: {e}")
