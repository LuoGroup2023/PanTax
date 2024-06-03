#!/usr/bin/env python3
import argparse, sys, json, h5py
import pandas as pd
from tqdm import tqdm

usage = "alignment json process"

def main():
    parser = argparse.ArgumentParser(prog="strain_abundance_cal.py", description=usage)
    parser.add_argument("read_cls_file", type=str, help="read classification file")
    parser.add_argument("aln_json_file", type=str, help="alignment json file")
    args = parser.parse_args()
    read_group_data = read_group(args.read_cls_file, args.aln_json_file)
    # print(sys.getsizeof(read_group_data))
    

def read_group(read_cls_file, aln_json_file):
    read_cls = pd.read_csv(read_cls_file, sep="\t", header=None)
    try:
        read_cls.columns = ["read_id", "mapq", "species_taxid"]
    except:
        read_cls.columns = ["read_id", "mapq", "species_taxid", "read_len"]
    read_cls = read_cls.dropna(subset=["species_taxid"])
    read_cls["species_taxid"] = read_cls["species_taxid"].astype(int).astype(str)
    read_id2species_taxids = read_cls.set_index("read_id")["species_taxid"].to_dict()
    species_taxids = set(read_cls["species_taxid"].tolist())
    read_group_data = {}
    for species_taxid in species_taxids:
        if species_taxid != "0":
            read_group_data[species_taxid] = []
    reads_aln_info = read_process(aln_json_file)
    for read_id in read_cls["read_id"].tolist():
        species_taxid = read_id2species_taxids[read_id]
        try:
            read_info = reads_aln_info[read_id]
            read_group_data[species_taxid].append(read_info)
        except:
            # print(f"{read_id} does not exist in json file.")
            continue
    return read_group_data

def read_process(aln_json_file):
    print("Processing alignments...")
    reads_aln_info = {}
    with open(aln_json_file, 'r') as aln_json:
        for line in tqdm(aln_json):
            aln = json.loads(line)
            seq_id = aln["name"]           
            try:
                # mapq = aln["mapping_quality"]
                path = aln["path"]
                mapping = path["mapping"]
            except KeyError:
                # read unmapped
                continue
            each_read_aln_info = {}
            for node_info in mapping:
                position = node_info["position"]
                node_id = int(position["name"])
                aln_len = 0
                edit = node_info["edit"]
                for aln_piece in edit:
                    try:
                        from_len = int(aln_piece["from_length"])
                    except KeyError:
                        from_len = 0
                    try:
                        to_len = int(aln_piece["to_length"])
                    except KeyError:
                        to_len = 0
                    aln_len += min(from_len, to_len)
                each_read_aln_info[node_id] = aln_len
            reads_aln_info[seq_id] = each_read_aln_info
    return reads_aln_info


if __name__ == "__main__":
    sys.exit(main())