import logging
import sys

class MainOnlyInfoFilter(logging.Filter):
    def filter(self, record):
        return not (record.levelno == logging.INFO and record.name != "__main__")

def setup_logging(level=logging.WARN):
    handler = logging.StreamHandler(sys.stdout)
    handler.addFilter(MainOnlyInfoFilter())
    formatter = logging.Formatter("[%(levelname)s] %(name)s: %(message)s")
    handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()
    logger.addHandler(handler)
    logger.propagate = False
