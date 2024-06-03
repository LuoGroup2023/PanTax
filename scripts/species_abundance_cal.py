#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
from staticsData import statics_and_write
import concurrent.futures
import subprocess

usage = "Compute species abundance"

def main():
    parser = argparse.ArgumentParser(prog="species_abundance_cal.py", description=usage)
    parser.add_argument("genomes_info", type=str, help="Genomes information file.Format: genomeID\tstrain_taxid\tspecies_taxid\ttorganism_name\tid.")
    parser.add_argument("read_cls_file", type=str, help="Read classification file")
    parser.add_argument("read_type", type=str, help="long/short")
    parser.add_argument("-ft", "--isfilter", dest="isfilter", type=int, default=1, help="MAPQ-based filter")
    parser.add_argument("wd", type=str, default=None, help="Work directory")
    args = parser.parse_args()
    global wd
    wd = args.wd
    species_genomes_len = "species_genomes_stats.txt"
    if not os.path.exists(f"{species_genomes_len}"):
        species_genomes_len_cal(f"{args.genomes_info}")
    genomes_stats = pd.read_csv(f"{species_genomes_len}", sep="\t", header=None, usecols=[0,1], dtype={0: str, 1: float, 2:float})
    genomes_stats.columns = ["species_taxid", "avg_len"]
    if args.read_type == "short":
        if args.isfilter:
            counts, read_len = short_filter_read_reads_cls(args.read_cls_file)
        else:
            counts, read_len = short_read_reads_cls(args.read_cls_file)
        abundance_set = short_abundance_cal(counts, genomes_stats, read_len)
    elif args.read_type == "long":
        abundance_set = long_read_reads_cls(args.read_cls_file, genomes_stats, args.isfilter)
    write(abundance_set)

def short_read_reads_cls(cls_file_path, sep="\t", header=None, usecol=2):
    reads_cls_data = pd.read_csv(cls_file_path, sep=sep, header=header, usecols=[usecol],dtype=object)
    reads_cls_data.iloc[:,0] = reads_cls_data.iloc[:,0].fillna("0")
    counts = reads_cls_data.iloc[:,0].value_counts().sort_values(ascending=False)
    read_len = pd.read_csv(cls_file_path, header=None, sep="\t", nrows=1).iloc[0,3]
    return counts, read_len

def short_filter_read_reads_cls(cls_file_path):
    reads_cls_data = pd.read_csv(cls_file_path, sep="\t", header=None, usecols=[1, 2],dtype={0:str, 1:int, 2:str, 3:int})
    reads_cls_data.columns = ["mapq", "species"]
    reads_cls_data = reads_cls_data.dropna(subset=["species"])
    grouped = reads_cls_data.groupby("species")
    index = []
    value_counts = [] 
    for group_name, group_data in grouped:
        read_count = len(group_data)
        # muti_match_read_count = len(group_data[group_data["mapq"] <= 1])
        less_muti_match_read_count = len(group_data[(group_data["mapq"] >= 3) & (group_data["mapq"] <= 60)])
        uniq_read_count = len(group_data[group_data["mapq"] == 60])
        if uniq_read_count == 0 or less_muti_match_read_count <= read_count/10:
            continue
        else:
            index.append(group_name)
            value_counts.append(read_count)
    counts = pd.Series(value_counts, index=index)
    read_len = pd.read_csv(cls_file_path, header=None, sep="\t", nrows=1).iloc[0,3]
    return counts, read_len

def long_read_reads_cls(cls_file_path, genomes_stats, isfilter):
    reads_cls_data = pd.read_csv(cls_file_path, sep="\t", header=None, dtype={0:str, 1:int, 2:str, 3:int})
    reads_cls_data.columns = ["readID", "mapq", "species", "read_len"]
    reads_cls_data = reads_cls_data.drop_duplicates(subset='readID', keep='first')
    reads_cls_data = reads_cls_data.dropna(subset=["species"])
    grouped = reads_cls_data.groupby("species")
    counts = {}
    for group_name, group_data in grouped:
        if isfilter:
            read_count = len(group_data)
            # muti_match_read_count = len(group_data[group_data["mapq"] <= 1])
            less_muti_match_read_count = len(group_data[(group_data["mapq"] >= 3) & (group_data["mapq"] <= 60)])
            uniq_read_count = len(group_data[group_data["mapq"] == 60])
            if uniq_read_count == 0 or less_muti_match_read_count <= read_count/10:
                continue
        read_size = group_data["read_len"].sum()
        counts[group_name] = read_size            
    abundance_set = {}
    for species_taxid, read_size in counts.items():
        if species_taxid == "0":
            continue
        genome_len = statics_species_genome_len(species_taxid, genomes_stats)
        if not genome_len:
            continue
        coverage = read_size/genome_len
        abundance_set[species_taxid] = coverage
    suml = sum(abundance_set.values())
    for key, value in abundance_set.items():
        abundance_set[key] = [value/suml, value]
    abundance_set = dict(sorted(abundance_set.items(), key=lambda x: x[1][0], reverse=True))
    return abundance_set

def statics_species_genome_len(species_taxid, genomes_stats):
    genome_len = None
    find_genome = genomes_stats.loc[genomes_stats["species_taxid"] == species_taxid]
    if not find_genome.empty:
        genome_len = find_genome["avg_len"].values[0]   
    return genome_len

def short_abundance_cal(counts, genomes_stats, each_read_size):
    abundance_set = {}
    for i in range(len(counts)):
        species_taxid = counts.index[i]
        if species_taxid == "0":
            abundance_set[species_taxid] = 0
            continue
        count = counts.values[i]
        genome_len = statics_species_genome_len(species_taxid, genomes_stats)
        if not genome_len:
            continue
        coverage = count*each_read_size/genome_len
        abundance_set[species_taxid] = coverage
    suml = sum(abundance_set.values())
    for key, value in abundance_set.items():
        abundance_set[key] = [value/suml, value]
    abundance_set = dict(sorted(abundance_set.items(), key=lambda x: x[1][0], reverse=True))
    return abundance_set

def species_genomes_len_cal(genomes_info_file):
    genomes_info = pd.read_csv(genomes_info_file, sep="\t", usecols=[2,4], dtype=object)       
    if not os.path.exists("genome_statics.txt"):
        genomes = genomes_info["id"].tolist()
        genomes = [os.path.join(wd, genome) if not os.path.isabs(genome) and wd else genome for genome in genomes]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(statics_and_write, genomes)    
    genomes_info["genomes"] = genomes_info["id"].apply(os.path.basename)
    genomes_statics = pd.read_csv("genome_statics.txt",sep="\t", header=None, usecols=[0,2])
    genomes_statics.drop_duplicates(inplace=True)
    genomes_statics.columns = ["genomes", "genome_length"]
    merge_df = pd.merge(genomes_info, genomes_statics, on="genomes", how="left")
    grouped = merge_df.groupby("species_taxid")
    species_taxid2avg_len = []
    for species_taxid, group_data in grouped:
        avg_len = group_data["genome_length"].mean()
        species_taxid2avg_len.append((species_taxid, avg_len))
    species_taxid2avg_len_df = pd.DataFrame(species_taxid2avg_len)
    species_taxid2avg_len_df.columns = ["species_taxid", "avg_len"]
    species_taxid2avg_len_df.to_csv(f"species_genomes_stats.txt", sep="\t", header=False, index=False)
    subprocess.run("rm genome_statics.txt", shell=True)

def write(abundance_set):
    with open("species_abundance.txt", "w") as f:
        f.write("species_taxid\tpredicted_abundance\tpredicted_coverage\n")
        for key, value in abundance_set.items():
            f.write(f"{key}\t{value[0]}\t{value[1]}\n")

if __name__ == "__main__":
    sys.exit(main())