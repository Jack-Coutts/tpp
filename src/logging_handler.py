import logging

from PySide6.QtCore import QObject, Signal


class QtLogHandler(logging.Handler, QObject):
    new_log_message = Signal(str)

    def __init__(self) -> None:
        logging.Handler.__init__(self)
        QObject.__init__(self)
        self.setFormatter(logging.Formatter("%(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        msg = self.format(record)
        self.new_log_message.emit(msg)
