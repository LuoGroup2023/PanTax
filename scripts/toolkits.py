import logging
from logging import handlers

class Logger(object):
    level_relations = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "crit": logging.CRITICAL
    }

    def __init__(self, level="debug",
                 fmt="%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s"):
        self.logger = logging.getLogger()
        format_str = logging.Formatter(fmt)
        self.logger.setLevel(self.level_relations.get(level))
        sh = logging.StreamHandler()
        sh.setFormatter(format_str)
        self.logger.addHandler(sh)