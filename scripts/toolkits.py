import logging, shutil
from logging import handlers
from time import time
from pathlib import Path

# runtime
def timeit(func):
    def func_wrapper(*args, **kwargs):
        start = time()
        ret = func(*args, **kwargs)
        end = time()
        spend = end - start
        print("{} cost time: {:.3f} s".format(func.__name__, spend))
        return ret

    return func_wrapper

# logger
class Logger(object):
    level_relations = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "crit": logging.CRITICAL
    }

    def __init__(self, level="debug",
                 fmt="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s"):
        self.logger = logging.getLogger()
        format_str = logging.Formatter(fmt)
        self.logger.setLevel(self.level_relations.get(level))
        sh = logging.StreamHandler()
        sh.setFormatter(format_str)
        self.logger.addHandler(sh)

# a class can accept all para
class Config:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return str(self.__dict__) 

# file process
def delete_directory_contents(directory_path):
    path = Path(directory_path)
    if path.is_dir():
        shutil.rmtree(path)  
        path.mkdir()  
    else:
        raise IOError(f"{directory_path} is not valid.")
    
def count_lines(filepath):
    with open(filepath, 'r') as f:
        return sum(1 for line in f)
    
def is_file_non_empty(file_path):
    return file_path.exists() and file_path.stat().st_size > 0