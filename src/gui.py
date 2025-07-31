import logging

from PySide6.QtCore import QThread, Slot
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .logging_handler import QtLogHandler
from .utils import check_directory, check_filepath  # Import utilities
from .worker import Worker


class TPPPlotterGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.worker = None
        self.thread = None
        self.setWindowTitle("TPP Plotter")
        self.setMinimumSize(800, 600)
        self.init_ui()
        self.setup_logging()

    def init_ui(self):
        # Widgets
        self.data_label = QLabel("Data File:")
        self.data_path = QLineEdit()
        self.data_browse = QPushButton("Browse")

        self.output_label = QLabel("Output Folder:")
        self.output_path = QLineEdit()
        self.output_browse = QPushButton("Browse")

        # radio buttons for multi-line plotting
        self.multi_line = QButtonGroup(self)
        self.single_line = QRadioButton("Plot One Line")
        self.two_lines = QRadioButton("Plot Two Lines")
        self.single_line.setChecked(True)

        self.multi_line.addButton(self.single_line)
        self.multi_line.addButton(self.two_lines)

        self.line_name_label = QLabel("Compound Name:")
        self.name_of_line = QLineEdit()

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

        # error bars checkbox
        self.error_bars_checkbox = QCheckBox("Add Error Bars")

        self.run_button = QPushButton("Run")

        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setFixedHeight(150)

        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)  # Disabled until processing starts

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

        # multi line
        multi_line = QVBoxLayout()
        multi_line.addWidget(self.single_line)
        multi_line.addWidget(self.two_lines)
        main_layout.addLayout(multi_line)

        # single line name
        line_name = QHBoxLayout()
        line_name.addWidget(self.line_name_label)
        line_name.addWidget(self.name_of_line)
        main_layout.addLayout(line_name)

        # horizonatal divider
        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.HLine)  # Set the shape to a horizontal line
        divider.setFrameShadow(
            QFrame.Shadow.Sunken
        )  # Set the shadow to make it visible
        main_layout.addWidget(divider)

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

        # horzontal divider
        divider1 = QFrame()
        divider1.setFrameShape(QFrame.Shape.HLine)  # Set the shape to a horizontal line
        divider1.setFrameShadow(
            QFrame.Shadow.Sunken
        )  # Set the shadow to make it visible
        main_layout.addWidget(divider1)

        # Error bars
        main_layout.addWidget(self.error_bars_checkbox)

        # Run button
        main_layout.addWidget(self.run_button)

        # logging box
        main_layout.addWidget(self.status_box)

        # stop button
        main_layout.addWidget(self.stop_button)

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
        self.single_line.toggled.connect(self.update_line_num_inputs)
        self.two_lines.toggled.connect(self.update_line_num_inputs)

        self.stop_button.clicked.connect(self.on_stop)

        self.update_mode_inputs()

    @Slot(str)
    def _append_log_message(self, message):
        self.status_box.append(message)

    def setup_logging(self):
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

    def browse_data_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Data File",
            "",
            "CSV or TSV Files (*.csv *.tsv *.txt);;All Files (*.*)",
        )
        if file_path:
            self.data_path.setText(file_path)

    def browse_output_folder(self):
        folder_path = QFileDialog.getExistingDirectory(self, "Select Output Folder")
        if folder_path:
            self.output_path.setText(folder_path)

    def browse_gene_list_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Gene List File", "", "Text Files (*.txt);;All Files (*.*)"
        )
        if file_path:
            self.gene_list_path.setText(file_path)

    def update_mode_inputs(self):
        # Show/hide input widgets based on selected mode
        self.gene_label.setEnabled(self.mode_gene.isChecked())
        self.gene_input.setEnabled(self.mode_gene.isChecked())

        self.gene_list_label.setEnabled(self.mode_gene_list.isChecked())
        self.gene_list_path.setEnabled(self.mode_gene_list.isChecked())
        self.gene_list_browse.setEnabled(self.mode_gene_list.isChecked())

    def update_line_num_inputs(self):
        # Show/hide input widgets based on selected line num
        self.line_name_label.setEnabled(self.single_line.isChecked())
        self.line.setEnabled(self.single_line.isChecked())

    def reset_ui(self):
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        logging.info("App reset to inactive state.")

    def on_stop(self):
        if hasattr(self, "worker") and self.worker:
            self.worker.request_stop()
            logging.info("Stop requested... waiting for processing to halt.")
            self.stop_button.setEnabled(
                False
            )  # Disable immediately to prevent multiple clicks

    def on_run(self):
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
                print(line_name)
            print(line_name)

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
            elif self.mode_filter.isChecked():
                mode = "filter"
            elif self.mode_all.isChecked():
                mode = "all"

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
    def _on_cleanup_finished(self):
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
