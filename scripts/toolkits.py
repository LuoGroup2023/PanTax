import logging, shutil, re, sys
# from logging import handlers
# from multiprocessing import Queue
from time import time
from pathlib import Path
import numpy as np
from typing import Dict, List
# from logging.handlers import QueueHandler, QueueListener

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

# def set_stream_handler(formatter: logging.Formatter):
#     handler = logging.StreamHandler(sys.stdout)
#     handler.setLevel(logging.WARNING)
#     handler.setFormatter(formatter)
#     return handler

# def set_queue_handler(log_queue):
#     handler = QueueHandler(log_queue)
#     handler.setLevel(logging.DEBUG)
#     return handler

# def get_logger(level: int = logging.DEBUG):
#     log_queue = Queue(-1)
#     logger = logging.getLogger()
#     logger.setLevel(level)
#     formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
#     stream_handler = set_stream_handler(formatter)
#     queue_handler = set_queue_handler(log_queue)
#     logger.addHandler(queue_handler)
#     global queue_listener
#     queue_listener = QueueListener(log_queue, stream_handler, respect_handler_level=True)
#     print(f"queue_listener_ori: {queue_listener}")
#     queue_listener.start()
#     return logger
# multi_logger = get_logger()

# def close_log_queue():
#     if queue_listener:
#         queue_listener.stop()

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

def check_file_avail(file_path):
    absolut_file_path = Path(file_path).resolve()
    if is_file_non_empty(absolut_file_path):
        return absolut_file_path
    else:
        raise FileNotFoundError(f"{absolut_file_path} does not exist, please check the path.")

def check_dir_avail(dir):
    dir = Path(dir).resolve()
    dir.mkdir(exist_ok=True)
    return dir

# genome process
def extract_genome_name(genome, gtdb):
    if gtdb:
        genome_ID = "_".join(str(Path(genome).name).split("_")[:2])
    else:
        genome_ID = re.sub(r"_genomic\.fna(\.gz)?$", "", str(Path(genome).name))
    return genome_ID