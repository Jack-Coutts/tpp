import sys

from PySide6.QtWidgets import QApplication

from src.gui import TPPPlotterGUI


def main() -> None:
    app = QApplication(sys.argv)
    window = TPPPlotterGUI()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
