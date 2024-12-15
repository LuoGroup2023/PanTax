#!/usr/bin/env python3
from toolkits import read_gfa, write_h5py_file
import sys

def main():
    gfa_file = sys.argv[1]
    h5_out_file = sys.argv[2]
    paths, nodes_len_npy  = read_gfa(gfa_file)
    write_h5py_file(nodes_len_npy, paths, h5_out_file)    

if __name__ == "__main__":
    sys.exit(main())
