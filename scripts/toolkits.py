import logging, shutil, h5py, re, sys
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

def read_gfa(graph_name: str, previous: int = 0):
    """
    Reads a graph from a GFA-file and returns graph in gt-format.
    """
    nodes_len_list: List[int] = []
    paths: Dict[str, List[int]] = {} # key = haplotype_name, vlaue = [node1, node2, ...]
    with open(graph_name, "r") as f:
        i = 0
        for line in f:
            if line.startswith("S"):
                line = line.rstrip('\n').split('\t')              
                # vertex
                node_id = int(line[1]) - 1 - previous
                assert i == node_id
                i += 1
                node_len = len(line[2])
                assert node_len != 0, "Error: Node length 0 appears in the GFA!"
                nodes_len_list.append(node_len)
            elif line[0] == "W" or line[0] == "P":
                reverse_flag = None
                haplotype_id = None
                line = line.rstrip('\n').split('\t')
                # W line               
                if line[0] == "W":
                    # haplotype_id = line[1] + "_" + line[3]
                    haplotype_id = line[1]
                    if line[-1].startswith("<"):
                        reverse_flag = True
                    else:
                        reverse_flag = False
                    path = [int(match.group())-1-previous for match in re.finditer(r'-?\d+', line[-1])]
                # P line
                elif line[0] == "P":
                    haplotype_id = line[1].split("#", 1)[0]
                    if line[2].split(",", 1)[0].endswith("-"):
                        reverse_flag = True
                    else:
                        reverse_flag = False
                    path = [int(match.group())-1-previous for match in re.finditer(r'\d+', line[2])]
                else:
                    continue
                if reverse_flag:
                    path = list(reversed(path))
                # Multiple chromosomes from the same genome merge into one path. 
                # Although this would introduce two non-existent trio nodes for each additional chromosome, 
                # it is not worth mentioning that only two nodes are relative to all nodes in the entire graph.
                if haplotype_id not in paths:
                    paths[haplotype_id] = path 
                else:
                    paths[haplotype_id].extend(path)             
    return paths, np.array(nodes_len_list)

def write_h5py_file(nodes_len_npy, paths, file_name):
    """
    Save the node length and path information of the graph in h5py format 
    """
    with h5py.File(file_name, "w") as hf:
        path_grp = hf.create_group("paths")
        for key,value in paths.items():
            path_grp.create_dataset(key, data=np.array(value))
        node_len_grp = hf.create_group("node_len")
        node_len_grp.create_dataset("node_len", data=nodes_len_npy)

def read_h5py_file(file_name):
    """
    Get the node length and path information of the graph from h5py format if exists
    """
    with h5py.File(file_name, "r") as hf:
            paths_grp = hf["paths"]
            paths = {key: list(value[:]) for key, value in paths_grp.items()}
            node_len_grp = hf["node_len"]
            nodes_len_npy = node_len_grp["node_len"][:]
    return paths, nodes_len_npy