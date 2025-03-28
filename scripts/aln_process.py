#!/usr/bin/env python3
import argparse, sys
import pandas as pd
from typing import Dict, List
from toolkits import timeit

usage = "alignment process"

def main():
    parser = argparse.ArgumentParser(prog="strain_abundance_cal.py", description=usage)
    parser.add_argument("read_cls_file", type=str, help="read classification file")
    parser.add_argument("aln_file", type=str, help="alignment file(GAF format)")
    args = parser.parse_args()
    read_group_data = read_gaf(args.read_cls_file, args.aln_file)
    # print(sys.getsizeof(read_group_data))
    
# @timeit
def read_gaf(read_cls_file: str, gaf_file: str) -> Dict[str, List[List[str]]]:
    read_cls = pd.read_csv(read_cls_file, sep="\t", header=None, dtype={0:str,1:int,2:str})
    try:
        read_cls.columns = ["read_id", "mapq", "species_taxid"]
    except:
        read_cls.columns = ["read_id", "mapq", "species_taxid", "read_len"]
    read_cls = read_cls.dropna(subset=["species_taxid"])
    # read_cls["species_taxid"] = read_cls["species_taxid"].astype(int).astype(str)
    alignment_not_paired_read_id = []
    grouped = read_cls.groupby("read_id")
    for group_id, group in grouped:
        if len(group["species_taxid"].unique()) != 1:  
            alignment_not_paired_read_id.append(group_id)
    read_id2species_taxids = read_cls.set_index("read_id")["species_taxid"].to_dict()
    species_taxids = set(read_cls["species_taxid"].tolist())
    read_group_data: Dict[str, List[List[str]]] = {} # key = species_taxid, value = [[read_id, read_path, read_path_len, read_start, read_end], ...]
    for species_taxid in species_taxids:
        if species_taxid != "0":
            read_group_data[species_taxid] = []    
    reads_aln_info: Dict[str, List[str]] = {} # key = read_id, value = [read_id, read_path, read_path_len, read_start, read_end]
    # first method: line by line
    dup_flag = False
    with open(gaf_file, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")
            if tokens[8] != "*" and tokens[9] != "*":
                if tokens[0] in reads_aln_info:
                    read_id = tokens[0] + "_2"
                    reads_aln_info[read_id] = [tokens[0]] + tokens[5:9]
                    dup_flag = True
                else:    
                    reads_aln_info[tokens[0]] = [tokens[0]] + tokens[5:9] #[read_id, read_path, read_path_len, read_start, read_end]
    for read_id in read_cls["read_id"].tolist():
        if read_id in alignment_not_paired_read_id: continue
        species_taxid = read_id2species_taxids[read_id]
        try:
            read_info = reads_aln_info[read_id]
            read_group_data[species_taxid].append(read_info)
            if dup_flag:
                second_read_id = read_id + "_2"
                read_info = reads_aln_info.get(second_read_id, None)
                if read_info: read_group_data[species_taxid].append(read_info)
        except:
            print(f"{read_id} does not exist in GAF file.")
            continue 
    # print(f"read_gaf read_group_data size: {sys.getsizeof(read_group_data)} Bytes")
    return read_group_data   
   

if __name__ == "__main__":
    sys.exit(main())