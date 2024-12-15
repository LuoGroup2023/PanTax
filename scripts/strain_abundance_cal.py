#!/usr/bin/env python3
import sys, os, re, argparse, time, random, logging
import numpy as np
import pandas as pd
from gurobipy import *
from aln_process import read_gaf
import concurrent.futures 
from pathlib import Path
from toolkits import timeit, Logger, read_gfa, read_h5py_file
from typing import Dict, List, Tuple, Union
from functools import partial
from collections import deque

random.seed(42)

usage = "Compute strain abundance"

def main():
    parser = argparse.ArgumentParser(prog="strain_abundance_cal.py", description=usage)
    parser.add_argument("db", type=str, help="pantax database directory")
    parser.add_argument("genomes_info", type=str, help="Genomes information file")
    parser.add_argument("read_cls", type=str, help="Reads classification file")
    parser.add_argument("aln_file", type=str, help="Reads alignment file(GAF format)")
    parser.add_argument("otu_range_file", type=str, help="OTUs range in graph")
    parser.add_argument("species_abundance_file", type=str, help="Species abundance file")
    parser.add_argument("-sq", "--sylph_query", dest="sylph_query_file", type=str, help="Sylph query file")
    parser.add_argument("-c", "--min_cov", dest="min_cov", type=float, default=0, help="Minimum coverage required per strain")
    parser.add_argument("-d", "--min_depth", dest="min_depth", type=int, default=0, help="Output a list of nodes with sequence depth less than <min_depth>.")
    parser.add_argument("-fr", "--unique_trio_nodes_fraction", dest="unique_trio_nodes_fraction", type=float, default=0.3, help="Unique trio nodes fraction")
    parser.add_argument("-fc", "--unique_trio_nodes_count", dest="unique_trio_nodes_count", type=float, default=0.45, help="Unique trio nodes mean count fraction")
    parser.add_argument("-a", "--min_species_abundance", dest="min_species_abundance", type=float, default=1e-04, help="Minimum species abundance(filter low abundance species)")
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=64, help="Set number of threads used for species.")
    parser.add_argument("-gt", "--gurobi_threads", dest="gurobi_threads", type=int, default=1, help="Set number of threads used for Gurobi.")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="Show log")
    parser.add_argument("--sample", dest="sample", type=int, default=None, help="Sampling nodes(500000) are used for the graph with too many nodes.")
    parser.add_argument("--sample_test", dest="sample_test", type=int, default=0, help="Sampling nodes are used for small model testing")
    # for test
    parser.add_argument("-sd", "--single_cov_diff", dest="single_cov_diff", type=float, default=0.2, help="Coverage difference for the species with single strain.")
    parser.add_argument("-m", "--mode", dest="mode", type=str, default="all", help="Select output pangenome result.")
    parser.add_argument("-p", "--parallel", dest="parallel", type=str, default="true", help="Whether parallel POA.")
    parser.add_argument("-ds", "--designated_species", dest="designated_species", type=str, default=None, help="Only return designated species result.")
    args = parser.parse_args()

    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()

    strain_abundance_est = StrainAbundanceEst(**vars(args))
    strain_abundance_est.run()

class StrainAbundanceEst():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.set_paras()

    def set_paras(self):
        global log
        if self.verbose:
            log = Logger(level="debug").logger
        else:
            log = Logger(level="info").logger
        global mode, verbose, sample_test, sample
        mode = self.mode    
        verbose = self.verbose
        sample_test = self.sample_test    
        sample = self.sample
        if self.parallel.lower() == "true":
            self.parallel = True
        else:
            self.parallel = False
        self.minimization_min_cov = 0
        if self.designated_species != "0" and self.designated_species:
            self.designated_species = self.designated_species.strip().split(",")

    def prepare_before_poa(self):
        self.otu_to_range: Dict[str, List[str]] = {} # [key] = species_taxid, [value] = [start,end] -> the start and end position in reference_pangenome.gfa
        self.pangenome_eq1 = []
        self.pangenome_ge2 = []
        self.read_cls_file()
        self.worker_num = self.threads // self.gurobi_threads
        if self.worker_num > self.gurobi_threads and self.worker_num > len(self.otu_to_range):
            factor = self.worker_num // len(self.otu_to_range)
            self.gurobi_threads = self.gurobi_threads * factor
            self.worker_num = self.threads // self.gurobi_threads
        else:
            self.worker_num = self.threads
        log.info("Reading GAF file...")
        global read_group_data
        read_group_data = self.read_group_data = read_gaf(self.read_cls, self.aln_file)
        log.info("Reading GAF file...done")
        self.otu_cov= [] # [(species_taxid, strain_record_name(i.e. GCF...), coverage)...]

    def run_poa(self):
        if self.mode == "all":
            otu_list = list(self.otu_to_range.keys())
        elif self.mode == "single":
            otu_list = self.pangenome_eq1
        elif self.mode == "double":
            otu_list = self.pangenome_ge2
        if self.designated_species:
            ds_otu_range = {}
            for _ds_otu in self.designated_species:
                if _ds_otu in self.otu_to_range:
                    ds_otu_range[_ds_otu] = self.otu_to_range[_ds_otu]
                else:
                    log.info(f"{_ds_otu} does not exist in graph.")
            if ds_otu_range:
                self.otu_to_range = ds_otu_range
        count = 0
        current_percentage = last_percentage = 0
        all_species = len(otu_list)
        partial_process = partial(optimize_otu, 
                                db=self.db, 
                                min_depth=self.min_depth, 
                                minimization_min_cov=self.minimization_min_cov, 
                                min_cov=self.min_cov, 
                                unique_trio_nodes_fraction=self.unique_trio_nodes_fraction, 
                                unique_trio_nodes_count=self.unique_trio_nodes_count,
                                gurobi_threads=self.gurobi_threads)
        if self.parallel:
            log.info(f"Parallel execution of strain optimization(threads {self.threads}, worker num {self.worker_num}, gurobipy threads {self.gurobi_threads})...")
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.worker_num) as executor:
                futures = {executor.submit(partial_process, key, value): key for key, value in self.otu_to_range.items()}
                for future in concurrent.futures.as_completed(futures): 
                    result = future.result()
                    if result:
                        result = self.abundace_constraint(result)
                        # print(result)
                        self.otu_cov.extend(result)
                    count += 1
                    current_percentage = count * 100 / all_species
                    if current_percentage > last_percentage + 5:
                        log.info(f"Strain abundance estimate percentage: {round(current_percentage, 1)}%")
                        last_percentage = current_percentage
        else:
            log.info("Sequential execution of strain optimization...")
            for key, value in self.otu_to_range.items():
                result = partial_process(key, value)
                if result:
                    result = self.abundace_constraint(result)
                    self.otu_cov.extend(result)
                    count += 1
                    current_percentage = count * 100 / all_species
                    if current_percentage > last_percentage + 5:
                        log.info(f"Strain abundance estimate percentage: {round(current_percentage, 1)}%")
                        last_percentage = current_percentage

    def run(self):
        global_t1_start = time.perf_counter()
        global_t2_start = time.process_time()
        self.prepare_before_poa()
        self.run_poa()
        self.abundance_cal()        
        t3_stop = time.perf_counter()
        t4_stop = time.process_time()

        log.info("All strains optimize completed")
        log.info("Elapsed time: {:.1f} seconds".format(t3_stop-global_t1_start))
        log.info("CPU process time: {:.1f} seconds".format(t4_stop-global_t2_start))
        print()

    def read_cls_file(self):
        """
        Obtain species filtered for abundance, and then obtain their location information in the reference pangenome.
        """
        species_abundance = pd.read_csv(self.species_abundance_file, sep="\t", dtype={0:str, 1:float, 2:float})
        species_abundance.columns = ["species_taxid", "abundance", "coverage"]
        species_abundance = species_abundance[species_abundance["abundance"] > self.min_species_abundance]
        otu_list = species_abundance["species_taxid"].tolist()
        # otu_list = ["1589", "630"]
        # print(len(otu_list))
        with open(self.otu_range_file, "r") as f:
            for line in f:
                tokens = line.strip().split("\t")
                otu = tokens[0]
                if otu in otu_list:
                    self.otu_to_range[otu] = (tokens[1], tokens[2])
                    if len(tokens) == 4 and tokens[3] == "0":
                        self.pangenome_eq1.append(otu)
                    elif len(tokens) == 4 and tokens[3] == "1":
                        self.pangenome_ge2.append(otu)


    def abundace_constraint(self, otu_mapping_list):
        species_abundance_df = pd.read_csv(self.species_abundance_file, sep="\t")
        species_abundance_df.columns = ["species_taxid", "abundance", "avg_coverage"]
        species_abundance_df["species_taxid"] = species_abundance_df["species_taxid"].astype(str)
        species_taxid = otu_mapping_list[0][0]
        strains: List[str] = [] 
        strains_coverage: List[Union[int, float]] = []
        for mapping in otu_mapping_list:
            strains.append(mapping[1])
            strains_coverage.append(mapping[7])
        species_coverage = species_abundance_df[species_abundance_df["species_taxid"] == species_taxid]["avg_coverage"].tolist()[0]
        total_cov_diff = abs(sum(strains_coverage)-species_coverage)/((sum(strains_coverage)+species_coverage)/2)
        for i in range(len(otu_mapping_list)):
            otu_mapping_list[i].append(total_cov_diff)        
        if max(strains_coverage) > 1.05*species_coverage:
            factor = species_coverage/sum(strains_coverage)
            strains_coverage = np.array(strains_coverage)*factor
            for i in range(len(otu_mapping_list)):
                otu_mapping_list[i][7] = strains_coverage[i]
        return otu_mapping_list

    def abundance_cal(self) -> None:
        genomes_info = pd.read_csv(self.genomes_info, sep="\t",usecols=[0,1],dtype={"genome_ID": str, "strain_taxid": str})
        genomes_info["hap_id"] = genomes_info["genome_ID"].str.split("_").str[:2].str.join("_")
        otu_cov_df = pd.DataFrame(self.otu_cov, columns=["species_taxid", "hap_id", "unique_trio_fraction", "uniq_trio_cov_mean", "path_base_cov", "first_sol", "strain_cov_diff", "predicted_coverage", "total_cov_diff"])
        otu_cov_df = pd.merge(otu_cov_df, genomes_info, on="hap_id", how="left")
        otu_cov_df["predicted_abundance"] = otu_cov_df["predicted_coverage"] / otu_cov_df["predicted_coverage"].sum()
        new_order = ["species_taxid", "strain_taxid", "genome_ID", "predicted_coverage", "predicted_abundance", "path_base_cov", "unique_trio_fraction", "uniq_trio_cov_mean", "first_sol","strain_cov_diff", "total_cov_diff"]
        ori_otu_cov_df = otu_cov_df[new_order]
        ori_otu_cov_df.to_csv("ori_strain_abundance.txt", sep="\t", index=False)
        otu_cov_df = otu_cov_df.groupby("species_taxid").filter(
            lambda group: not (len(group) == 1 and group["total_cov_diff"].iloc[0] > self.single_cov_diff)
        )
        otu_cov_df = otu_cov_df[(otu_cov_df["predicted_coverage"] >= self.min_cov) & (otu_cov_df["predicted_coverage"] != 0)]
        otu_cov_df["predicted_abundance"] = otu_cov_df["predicted_coverage"] / otu_cov_df["predicted_coverage"].sum()
        new_otu_cov_df = otu_cov_df.drop("hap_id", axis=1)
        new_otu_cov_df = new_otu_cov_df[new_order]
        new_otu_cov_df = new_otu_cov_df.sort_values(by="predicted_abundance", ascending=False)
        # new_otu_cov_df = new_otu_cov_df.round(2)
        new_otu_cov_df.to_csv("strain_abundance.txt", sep="\t", index=False)


def optimize_otu(otu: str, otu_range: List[str], db: str, min_depth: int, 
                          minimization_min_cov: int, min_cov: float, unique_trio_nodes_fraction: float, 
                          unique_trio_nodes_count: float, gurobi_threads: int):
    # otu_range = self.otu_to_range[otu]
    start = int(otu_range[0])-1
    end = int(otu_range[1])-1
    log.debug(f"Reading {otu} information...\n")        

    species_gfa_dir = Path(db) / "species_gfa"
    species_graph_info_dir = Path(db) / "species_graph_info"
    if species_gfa_dir.exists():
        paths, nodes_len_npy  = read_gfa(f"{db}/species_gfa/{otu}.gfa")
    elif species_graph_info_dir.exists():
        paths, nodes_len_npy = read_h5py_file(f"{db}/species_graph_info/{otu}.h5")
    else:
        raise FileNotFoundError(f"{str(species_gfa_dir)} and {str(species_graph_info_dir)} are not found. No single speices GFA file information.")
    t1_start = time.perf_counter()
    unique_trio_nodes, unique_trio_nodes_len, hap2unique_trio_nodes_m = trio_nodes_info(paths, nodes_len_npy) 
    t1_stop = time.perf_counter()
    log.debug(f"{otu} obtain unique trio node information elapsed time spend: {round(t1_stop-t1_start, 1)} seconds")
    nodes_len = {node:node_len for node, node_len in enumerate(nodes_len_npy)}
    t1_start = time.perf_counter()
    node_abundances, unique_trio_node_abundances, node_base_cov, node_base_all = get_node_abundances(otu, nodes_len, unique_trio_nodes, unique_trio_nodes_len, start) 
    t1_stop = time.perf_counter()
    log.debug(f"{otu} obtain all node information elapsed time spend: {round(t1_stop-t1_start, 1)} seconds")
    non_zero_count = sum(1 for elem in node_abundances if elem != 0)
    log.debug(f"{otu} species node abundance > 0 number:{non_zero_count}\n")
    haps_id = list(paths.keys())
    nvert = end-start+1      
    otu_abundance_list = [x if x > min_depth else 0 for x in node_abundances]
    if non_zero_count == 0 or max(otu_abundance_list) == 0 :
        return [[otu, hap_id]+["-"]*4+[0] for hap_id in haps_id]
    # Now solve the minimization problem  
    result = optimize(otu_abundance_list, nvert, paths, hap2unique_trio_nodes_m, unique_trio_node_abundances, node_base_cov, node_base_all, minimization_min_cov, unique_trio_nodes_fraction, min_cov, unique_trio_nodes_count, gurobi_threads)
    if result:
        hap_metrics, objVal = result
        return [[otu] + _hap_metrics for _hap_metrics in hap_metrics]
    else:
        return None
    
def trio_nodes_info(paths: Dict[str, List[int]], nodes_len_npy: np.ndarray):
    trio_nodes: List[Tuple[int, int, int]] = []
    hap_trio_paths: Dict[str, List[Tuple[int, int, int]]] = {}
    # chop trio nodes for all paths
    for hap, path in paths.items():
        trio_path = [(path[i], path[i+1], path[i+2]) for i in range(len(path)-2)]
        hap_trio_paths[hap] = trio_path
        trio_nodes.extend(trio_path)
    # get unique trio nodes
    trio_nodes = list(set(trio_nodes))
    # save trio nodes whether it exists in every path
    trio_nodes_dict: Dict[Tuple[int, int, int], List[List[np.ndarray, int, int]]] = {}
    for idx, trio_node in enumerate(trio_nodes):
        trio_nodes_dict[trio_node] = [np.array([0] * len(hap_trio_paths)), 0, idx]

    i = 0
    haps = []
    for hap, trio_path in hap_trio_paths.items():
        haps.append(hap)
        for trio_node in trio_path:
            try:
                trio_nodes_dict[trio_node][0][i] = 1
                trio_nodes_dict[trio_node][1] += 1
                if i >= 2 and trio_nodes_dict[trio_node][1] == 2:
                    del trio_nodes_dict[trio_node]
            except:
                continue
        i += 1

    hap2unique_trio_nodes_v: List[np.ndarray] = []
    unique_trio_nodes_idx: List[int] = []
    for trio_node, values in trio_nodes_dict.items():
        hap2unique_trio_nodes_v.append(values[0])
        unique_trio_nodes_idx.append(values[2])
    try:
        hap2unique_trio_nodes_m = np.stack(hap2unique_trio_nodes_v)
        # trio nodes length
        unique_trio_nodes_len = [sum(nodes_len_npy[idx] for idx in trio_nodes[i]) for i in unique_trio_nodes_idx]
        unique_trio_nodes = {trio_nodes[i]: idx for idx, i in enumerate(unique_trio_nodes_idx)}
    except:
        print(f"Warnings: {paths} have the same path or the total number of nodes per path is less than 3.")
        unique_trio_nodes_len = []
        unique_trio_nodes = {}
        hap2unique_trio_nodes_m = np.array([])
    return unique_trio_nodes, unique_trio_nodes_len, hap2unique_trio_nodes_m

def get_node_abundances(otu: str, nodes_len: Dict[int, int], trio_nodes: Dict[Tuple[int, int, int], int], 
                            trio_nodes_len: List[int], start: int) -> Tuple[List[int], List[int]]:
    """
    For now we do not consider the insertions and deletions of the nodes. For the start node 
    and end node of a read, we get matched base which is equal to the node length minus the offset.
    For the intermediate nodes, we keep the node length as matched base. In contract, by parsing 
    json format or gam format alignment, the intermediate nodes' matched base will consider indels.
    In previous version, the matched base number is the node length minus the number of deletions, 
    but don't add the number of insertions(keep node length). 
    Actually, we can consider indels by the cs tag or cg tag in gaf file. 
    """
    bases_per_node = {} # map node IDs to read alignments
    node_base_cov_info = {}
    for i in range(len(nodes_len)):
        bases_per_node[i] = 0
        node_base_cov_info[i] = [0, 0, [0]*nodes_len[i]] #start, end, base_cov, is_all_cov
    trio_nodes_bases_count = {}
    for i in range(len(trio_nodes)):
        trio_nodes_bases_count[i] = 0
    log.debug("Processing alignments...")
    count = 0
    for read_info in read_group_data[otu]:
        # read_info -> [read_id, read_path, read_path_len, read_start, read_end]
        read_nodes = [int(match.group())-1-start for match in re.finditer(r'-?\d+', read_info[1])]
        # should keep the insertion order. Python version > 3.7 !!!
        read_nodes_len_in_graph = {node:nodes_len[node] for node in read_nodes}
        start_node = read_nodes[0]
        end_node = read_nodes[-1]
        # read_path_len = int(read_info[2])
        read_start = int(read_info[3])
        read_end = int(read_info[4])
        # if no indels, target_len = read_length
        target_len = read_end - read_start
        seen = 0
        read_nodes_len = {node:0 for node in read_nodes}
        undup_read_nodes = set()
        if start_node == end_node and len(read_nodes) == 1:
            count += 1
            read_nodes_len[start_node] += target_len
            bases_per_node[start_node] += target_len
            if node_base_cov_info[start_node][0] == 0:
                for i in range(read_start, read_end):
                    node_base_cov_info[start_node][2][i] = 1
                node_cov = sum(node_base_cov_info[start_node][2])
                node_base_cov_info[start_node][1] = node_cov
                if node_cov == nodes_len[start_node]:
                    is_all_cov = 1
                else:
                    is_all_cov = 0
                node_base_cov_info[start_node][0] = is_all_cov

        else:
            i = 1
            for node in read_nodes:
                node_len = read_nodes_len_in_graph[node]
                if i == 1:
                    if node_len < read_start:
                        print(f"OTU: {otu}, read_id: {read_info[0]}, i: {i}, node: {node}, node_len: {node_len}, read_start:{read_start}")
                        print(read_nodes)
                        print(read_nodes_len_in_graph)
                        print(read_nodes_len)                        
                    assert node_len >= read_start
                    node_aln_len = node_len - read_start
                    if node_base_cov_info[node][0] == 0:
                        for j in range(read_start, node_len):
                            node_base_cov_info[node][2][j] =  1

                elif i == len(read_nodes):
                    if target_len < seen: target_len = seen
                    node_aln_len = target_len - seen
                    if node_base_cov_info[node][0] == 0:
                        for j in range(0, node_aln_len - 1):
                            try:
                                node_base_cov_info[node][2][j] =  1
                            except:
                                print(len(node_base_cov_info[node][2]), node_aln_len, j)
                                print(f"OTU: {otu}, read_id: {read_info[0]}, i: {i}, node: {node}, node_len: {node_len}, read_start:{read_start}")
                                print(read_nodes)
                                print(read_nodes_len_in_graph)
                                print(read_nodes_len)
                                raise IndexError
                else:
                    node_aln_len = node_len
                    if node_base_cov_info[node][0] == 0:
                        for j in range(0, node_aln_len):
                            node_base_cov_info[node][2][j] =  1
                node_cov = sum(node_base_cov_info[node][2])
                node_base_cov_info[node][1] = node_cov
                if node_cov == nodes_len[node]:
                    is_all_cov = 1
                else:
                    is_all_cov = 0
                node_base_cov_info[node][0] = is_all_cov
                seen += node_aln_len
                if node not in undup_read_nodes:
                    read_nodes_len[node] += node_aln_len
                    bases_per_node[node] += node_aln_len
                    undup_read_nodes.add(node)
                i += 1
        if len(read_nodes) < 3:
            continue
        read_trio_nodes = [(read_nodes[i], read_nodes[i+1], read_nodes[i+2]) for i in range(len(read_nodes)-2)]
        try:
            read_trio_nodes_len = [sum(read_nodes_len[node] for node in trio_node) for trio_node in read_trio_nodes]
        except:
            print(f"OTU: {otu}, read_id: {read_info[0]}")
            print(read_nodes)
            print(read_nodes_len)
            print(read_trio_nodes)
            sys.exit(1)
        for i, read_trio_node in enumerate(read_trio_nodes):
            idx = trio_nodes.get(read_trio_node)
            if not idx:
                read_trio_node_reversed = tuple(reversed(read_trio_node))
                idx = trio_nodes.get(read_trio_node_reversed)
                if not idx:
                    continue
            trio_nodes_bases_count[idx] += read_trio_nodes_len[i]       
    node_abundance_list: List[int] = []
    log.debug("Computing node abundance rates...")
    for node, node_len in nodes_len.items():
        aligned_len = bases_per_node[node]
        if node_len > 0:
            node_abundance = aligned_len / node_len
        else:
            raise ZeroDivisionError(f"Node length 0 for node {node}")
        node_abundance_list.append(node_abundance)
    trio_node_abundance_list: List[int] = []
    log.debug("Computing trio node abundance rates...")
    for i, trio_node_len in enumerate(trio_nodes_len):
        aligned_len = trio_nodes_bases_count[i]
        if trio_node_len > 0:
            node_abundance = aligned_len / trio_node_len
        else:
            raise ZeroDivisionError(f"Trio node length 0 for node {node}")
        trio_node_abundance_list.append(node_abundance)

    node_base_cov = np.array([value[1] for value in node_base_cov_info.values()])
    node_base_all = np.array(list(nodes_len.values()))
    return node_abundance_list, trio_node_abundance_list, node_base_cov, node_base_all

def optimize(
        a: List[int], nvert: int, paths: Dict[str, List[int]], 
        hap2trio_nodes_m: np.ndarray, trio_node_abundances: List[int], node_base_cov, node_base_all,
        minimization_min_cov, unique_trio_nodes_fraction, min_cov, unique_trio_nodes_count, gurobi_threads
        ):
        """
        Defines Gurobi minimization problem and then applies the LP solver.
        Returns the solution values, the corresponding objective value, and a matrix
        counting pairwise differences between haplotypes.
        """
        max_strains = origin_paths_len = len(paths)
        if mode == "single" and origin_paths_len != 1:
            return None
        elif mode == "double" and origin_paths_len == 1:
            return None

        t1_start = time.perf_counter()
        t2_start = time.process_time()
        # set default value
        otu_paths = list(paths.values())
        haps_id = list(paths.keys())    
        hap_metrics = [[_hap]+["-"]*5+[0] for _hap in haps_id]
        trio_node_abundances = np.array(trio_node_abundances)
        size = hap2trio_nodes_m.size
        same_path_flag = False
        possible_strains_idx = []
        possible_strains_frequencies_mean = []
        
        if origin_paths_len != 1 and size:
            row_sums = np.sum(hap2trio_nodes_m, axis=1)
            for idx in range(len(haps_id)):
                selected_indices = np.where((hap2trio_nodes_m[:, idx] == 1) & (row_sums == 1))[0]
                if len(selected_indices) == 0: continue
                frequencies = np.zeros(len(trio_node_abundances))
                frequencies[selected_indices] = trio_node_abundances[selected_indices]
                count_greater_than_zero = np.sum(frequencies > 0)
                otu_unique_trio_nodes_fraction = count_greater_than_zero/len(selected_indices)
                if verbose: print(f"\t\t{haps_id[idx]} unique trio node abundance > 0 ratio: {otu_unique_trio_nodes_fraction}")
                hap_metrics[idx][1] = round(otu_unique_trio_nodes_fraction, 2)
                if otu_unique_trio_nodes_fraction < unique_trio_nodes_fraction: continue
                possible_strains_idx.append(idx)
                non_zero_frequencies = frequencies[frequencies > 0]
                non_zero_frequencies = filter_outliers(non_zero_frequencies, method='zscore')
                # np.save(f"{haps_id[idx]}_non_zero_frequencies.npy", non_zero_frequencies)
                frequencies_mean = np.mean(non_zero_frequencies) if len(non_zero_frequencies) > 0 else 0
                possible_strains_frequencies_mean.append(frequencies_mean)
                hap_metrics[idx][2] = round(frequencies_mean, 2)
            if verbose: print("\t\tFisrt filter #strains / #paths = {} / {}".format(len(possible_strains_idx), origin_paths_len))
            otu_paths = [otu_paths[idx] for idx in possible_strains_idx]
        elif origin_paths_len != 1 and size == 0:
            if all(x == otu_paths[0] for x in otu_paths):
                otu_paths = [otu_paths[0]]
                same_path_flag = True
                node_abundance = np.array(a)
                non_zero_frequencies = node_abundance[node_abundance > 0]
                frequencies_mean = np.mean(non_zero_frequencies) if len(non_zero_frequencies) > 0 else 0
                possible_strains_frequencies_mean.append(frequencies_mean)                
                possible_strains_idx = [0]
            else:
                possible_strains_idx = [_i for _i in range(origin_paths_len)]
                # raise ValueError(f"{origin_paths_len} more than 1, but have trio node is None.")
        elif origin_paths_len == 1:
            node_abundance = np.array(a)
            non_zero_frequencies = node_abundance[node_abundance > 0]
            frequencies_mean = np.mean(non_zero_frequencies) if len(non_zero_frequencies) > 0 else 0
            possible_strains_frequencies_mean.append(frequencies_mean)
            possible_strains_idx = [0]

        # t_stop = time.perf_counter()
        # print("     [PAO first filter based on fr] - Elapsed time: {:.1f} seconds".format(t_stop-t1_start))

        # absolute error
        m = Model('lp')
        obj = LinExpr()
        npaths = len(otu_paths)          #The number of feasible paths
        P = np.zeros((nvert,npaths))
        x = m.addVars(list(range(npaths)), lb=0, ub=1.05*max(a), vtype=GRB.CONTINUOUS, name='x')
        X = [x[i] for i in range(npaths)]
        X = np.array([X]).reshape(npaths,1)    #Set x in an array for multiplication
        nvert_list = list(range(nvert))

        # Store paths in P: p_ij = 1 if node i contains path j
        # print('\nSave for every node which paths are passing through:')
        # node_abundance = []
        for i in range(npaths):
            node_abundance = []
            for v in otu_paths[i]:
                P[int(v),i] = 1
                node_abundance.append(a[v])
            # np.save(f"{haps_id[i]}_node_abundance.npy", node_abundance)
                # if self.mode == "single":
                #     node_abundance.append(a[v])
        # log.debug(f"pangenome eq1 node abundance: {node_abundance}")

        path_cov = np.dot(node_base_cov, P)
        path_all_cov = np.dot(node_base_all, P)
        for _i, _x in enumerate(path_cov):
            ori_idx = possible_strains_idx[_i]
            hap_metrics[ori_idx][3] = _x / path_all_cov[_i]
        if sample_test:
            if nvert > 500:
                sample_nodes = random.sample(nvert_list, 500)
                sample_nodes = sorted(sample_nodes)
                nvert_list = sample_nodes
        if not sample_test and sample:
            if nvert > sample:
                log.debug(f"The graph has too many nodes. Subsample {sample}.")
                sample_nodes = random.sample(nvert_list, sample)
                sample_nodes = sorted(sample_nodes)
                nvert_list = sample_nodes

        # If objective involves absolute values, add extra variables
        y = m.addVars(nvert_list, lb=0, vtype=GRB.CONTINUOUS, name='y')
        # log.debug(f"nvert len:{len(nvert_list)}\ty len:{len(y)}")
        # add indicator variables to count strains
        if max_strains > 0:
            # add indicator variables for counting strains
            x_binary = m.addVars(list(range(npaths)), vtype=GRB.BINARY, name='strain_indicator')
            for i in range(npaths):
                max_x = 2*max(a)
                m.addConstr(x_binary[i] >= (x[i]-minimization_min_cov)/max_x)
            # define total strain count
            sum_x_binary = LinExpr()
            for i in range(npaths):
                sum_x_binary += x_binary[i]
            # bound the number of strains
            m.addConstr(sum_x_binary <= max_strains)

        del otu_paths, hap2trio_nodes_m, trio_node_abundances, path_cov, path_all_cov
        m.update()
        # Define the objective function
        # print('\nDefine the objective function:')
        n_eval = 0
        # for v in tqdm(range(nvert)):
        # nvert_list = deque(nvert_list)
        # for v in range(nvert):
        #     if (sample_test or sample) and nvert_list:
        #         # if v not in nvert_list: continue
        #         if v == nvert_list[0]:
        #             nvert_list.popleft()
        #         else:
        #             continue
        for v in nvert_list:
            # sum the calculated abundances of strains through v
            # sum_xv = 0
            # for idx, i in enumerate(P[v]):
            #     if i == 1:
            #         sum_xv += x[idx]

            # memory may be expensive!!!
            sum_xv = np.dot(P[v,:],X)[0]
            abundance = a[v]
            if abundance == 0:
                continue
            # absolute difference
            obj += y[v] #abs(abundance - sum_xv)
            # set constraints on y[v] to obtain absolute value
            m.addConstr(y[v] >= sum_xv - abundance, "y_{}_-".format(v))
            m.addConstr(y[v] >= abundance - sum_xv, "y_{}_+".format(v))
            n_eval += 1

        assert n_eval > 0
        
        obj *= (1/n_eval)
        # set objective and minimize
        m.setObjective(obj, GRB.MINIMIZE)
        # print("Objective function ready, starting Gurobi optimization:")
        # print()
        m.update()

        m.Params.LogToConsole = 0
        m.Params.Threads = gurobi_threads
        m.Params.NumericFocus = 0
        m.Params.PoolSearchMode = 0
        m.Params.PoolSolutions = 10
        m.Params.Method = 4

        #Minimize the model for the given objective function and constraints
        if verbose: print("\t\t*** Phase 1 optimization***")
        m.optimize()

        # print("Number of solutions = {}".format(m.solcount))
        if m.status == GRB.Status.OPTIMAL:
            x_sol = []
            for v in m.getVars():
                if 'x' in v.varName:
                    x_sol.append(v.x)
            if verbose: print(f"\t\tObjective value: {m.objVal}")
            objVal = m.objVal
            if verbose: print(f"\t\tFirst sol:{x_sol}\n")
            selected_strains = [1 if cov > 0 else 0 for cov in x_sol]
            nstrains = sum(selected_strains)
            for _i, _x in enumerate(x_sol):
                ori_idx = possible_strains_idx[_i]
                hap_metrics[ori_idx][4] = _x

            if verbose: print("\t\tFirst optimization #strains / #paths = {} / {}".format(nstrains, npaths))
            
            if origin_paths_len != 1 and size:            
                for idx, frequencies_mean in enumerate(possible_strains_frequencies_mean):
                    if frequencies_mean == 0:
                        selected_strains[idx] = 0
                        continue
                    f = abs(x_sol[idx]-frequencies_mean)/(x_sol[idx]+frequencies_mean)
                    ori_idx = possible_strains_idx[idx]
                    hap_metrics[ori_idx][5] = f = round(f, 2)
                    if verbose: print(f"\t\t{haps_id[ori_idx]}\tfrequencies_mean:{frequencies_mean}\txsol:{x_sol[idx]}\tf:{f}")
                    if f > unique_trio_nodes_count: selected_strains[idx] = 0
            elif (origin_paths_len != 1 and size == 0 and same_path_flag) or origin_paths_len == 1: 
                frequencies_mean = possible_strains_frequencies_mean[0]
                if frequencies_mean:   
                    f = abs(x_sol[0]-frequencies_mean)/(x_sol[0]+frequencies_mean)
                    if verbose: print(f"\t\tidx:{haps_id[0]}\tfrequencies_mean:{frequencies_mean}\txsol:{x_sol[0]}\tf:{f}")
                    hap_metrics[0][5] = round(f, 2)
                    hap_metrics[0][6] = hap_metrics[0][4]
                return hap_metrics, objVal
            elif origin_paths_len != 1 and size == 0 and not same_path_flag:
                for _i in range(origin_paths_len):
                    ori_idx = possible_strains_idx[_i]
                    hap_metrics[ori_idx][6] = hap_metrics[ori_idx][4]
                return hap_metrics, objVal                
        
            nstrains = sum(selected_strains)
            if verbose: print("\t\tSecond filter #strains / #paths = {} / {}".format(nstrains, npaths))

            # run phase 2 optimization:
            if verbose: print("\t\t*** Phase 2 optimization***")
            m.reset()
            for i in range(npaths):
                if selected_strains[i] == 0:
                    m.addConstr(x[i] == 0)
            m.optimize()

            x_final = []
            for v in m.getVars():
                if 'x' in v.varName:
                    x_final.append(v.x)
            if verbose: print(f"\t\tObjective value: {m.objVal}")
            objVal = m.objVal

            # selected_strains = [1 if cov > self.min_cov else 0 for cov in x_final]
            nstrains = sum(selected_strains)
            if verbose: print("\t\tSecond optimization #strains / #paths = {} / {}".format(nstrains, npaths))

            for _i, _x in enumerate(x_final):
                ori_idx = possible_strains_idx[_i]
                hap_metrics[ori_idx][6] = _x

            t1_stop = time.perf_counter()
            t2_stop = time.process_time()
            if verbose: 
                print("\nStrain path optimize completed")
                print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
                print("CPU process time: {:.1f} seconds\n".format(t2_stop-t2_start))
            return hap_metrics, objVal

        else:
            try:
                m.computeIIS()
                # Print the names of all of the constraints in the IIS set.
                print("IIS constraints:")
                for c in m.GetConstrs():
                    if c.Get(GRB.IntAttr.IISConstr) > 0:
                        print(c.Get(GRB.StringAttr.ConstrName))
                # Print the names of all of the variables in the IIS set.
                print("IIS variables:")
                for v in m.GetVars():
                    if v.Get(GRB.IntAttr.IISLB) > 0 or v.Get(GRB.IntAttr.IISUB) > 0:
                        print(v.Get(GRB.StringAttr.VarName))
                print("ERROR: Infeasible model.")
            except:
                #print(m.getAttr(GRB.Attr.UnbdRay, m.getVars()))
                print("ERROR: Unbounded model.")

            print('\nNo optimal solution found, exiting.')
            sys.exit(1)

def filter_outliers(data, method='iqr', threshold=3):
    if method == 'zscore':
        mean = np.mean(data)
        std = np.std(data)
        return [x for x in data if abs((x - mean) / std) < threshold]
    elif method == 'iqr':
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        return [x for x in data if lower_bound <= x <= upper_bound]
    else:
        raise ValueError("Invalid method. Choose 'zscore' or 'iqr'.")


if __name__ == "__main__":
    sys.exit(main())
