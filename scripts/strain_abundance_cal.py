#!/usr/bin/env python3
import sys, os, re, argparse, time, h5py
import numpy as np
import pandas as pd
from gurobipy import *
from tqdm import tqdm # progress tracker
from functools import partial
from aln_process import read_gaf
import concurrent.futures 
from toolkits import timeit, Logger
from typing import Dict, List, Tuple, Union
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

usage = "Compute strain abundance"

def main():
    parser = argparse.ArgumentParser(prog="strain_abundance_cal.py", description=usage)
    parser.add_argument("db", type=str, help="pantax database directory")
    parser.add_argument("genomes_info", type=str, help="Genomes information file")
    parser.add_argument("read_cls", type=str, help="Reads classification file")
    parser.add_argument("aln_file", type=str, help="Reads alignment file(GAF format)")
    parser.add_argument("otu_range_file", type=str, help="OTUs range in graph")
    parser.add_argument("species_abundance_file", type=str, help="Species abundancefile")
    parser.add_argument("-c", "--min_cov", dest="min_cov", type=float, default=0, help="Minimum coverage required per strain")
    parser.add_argument("-d", "--min_depth", dest="min_depth", type=int, default=0, help="Output a list of nodes with sequence depth less than <min_depth>.")
    parser.add_argument("-fr", "--unique_trio_nodes_fraction", dest="fr", type=float, default=0.3, help="Unique trio nodes fraction")
    parser.add_argument("-fc", "--unique_trio_nodes_count", dest="fc", type=float, default=0.45, help="Unique trio nodes mean count fraction")
    parser.add_argument("-a", "--min_species_abundance", dest="min_species_abundance", type=float, default=1e-04, help="Minimum species abundance(filter low abundance species)")
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=64, help="Set number of threads used for species.")
    parser.add_argument("-gt", "--gurobi_threads", dest="gurobi_threads", type=int, default=1, help="Set number of threads used for Gurobi.")
    parser.add_argument("-s", "--save_graph_info", dest="s", type=int, default=0, help="Save graph information")
    parser.add_argument("-v", "--verbose", dest="is_show_log", action="store_true", help="Show log")
    args = parser.parse_args()

    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()
    log = Logger()
    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()
    species_graph_info = os.path.join(args.db, "species_graph_info")
    if not os.path.exists(f"{species_graph_info}") and args.s:
        os.mkdir("species_graph_info")
    minimization_min_cov = 0
    otu_to_range = read_cls(args.otu_range_file, args.species_abundance_file, args.min_species_abundance)
    if len(otu_to_range) == 1:
        args.gurobi_threads = max(args.threads, args.gurobi_threads)
        args.threads = 1    
    global read_group_data, is_show_log
    is_show_log = args.is_show_log
    read_group_data = read_gaf(args.read_cls, args.aln_file)
    log.logger.info("Parallel estimation of abundance for each strain of every species(this step takes a long time)...")
    otu_cov: List[Tuple[str, str, Union[int, float]]] = [] # [(species_taxid, strain_record_name(i.e. GCF...), coverage)...]
    partial_process = partial(parallel_optimize_otu, 
                              pantax_db=args.db, 
                              min_depth=args.min_depth, 
                              minimization_min_cov=minimization_min_cov, 
                              min_cov=args.min_cov, 
                              unique_trio_nodes_fraction=args.fr, 
                              unique_trio_nodes_mean_count_fraction=args.fc,
                              s=args.s, 
                              threads=args.gurobi_threads)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(partial_process, key, value): key for key, value in otu_to_range.items()}
        for future in concurrent.futures.as_completed(futures): 
            result = future.result()
            result = abundace_constraint(result, species_abundance_file = args.species_abundance_file)
            # print(result)
            otu_cov.extend(result)
            # try:
            #     result = future.result()
            #     result = abundace_constraint(result, species_abundance_file = args.species_abundance_file)
            #     # print(result)
            #     otu_cov.extend(result)
            # except Exception as e:
            #     print(f"An error occurred in strain recall process: {e}")
    abundance_cal(otu_cov, args.genomes_info)

    t3_stop = time.perf_counter()
    t4_stop = time.process_time()

    log.logger.info("All strains optimize completed")
    log.logger.info("Elapsed time: {:.1f} seconds".format(t3_stop-global_t1_start))
    log.logger.info("CPU process time: {:.1f} seconds".format(t4_stop-global_t2_start))
    print()
    return

################################################################################

def write_h5py_file(nodes_len_npy: np.ndarray, paths: Dict[str, List[int]], file_name: str):
    """
    Save the node length and path information of the graph in h5py format 
    """
    with h5py.File(file_name, "w") as hf:
        path_grp = hf.create_group("paths")
        for key,value in paths.items():
            path_grp.create_dataset(key, data=np.array(value))
        node_len_grp = hf.create_group("node_len")
        node_len_grp.create_dataset("node_len", data=nodes_len_npy)

def read_h5py_file(file_name: str):
    """
    Get the node length and path information of the graph from h5py format if exists
    """
    with h5py.File(file_name, 'r') as hf:
            paths_grp = hf['paths']
            paths = {key: list(value[:]) for key, value in paths_grp.items()}
            node_len_grp = hf['node_len']
            nodes_len_npy = node_len_grp['node_len'][:]
    return paths, nodes_len_npy

def parallel_optimize_otu(otu: str, otu_range: List[str], pantax_db: str, min_depth: int, 
                          minimization_min_cov: int, min_cov: float, unique_trio_nodes_fraction: float, 
                          unique_trio_nodes_mean_count_fraction: float, s: int, threads: int) -> List[Tuple[str, str, Union[int, float]]]:
    start = int(otu_range[0])-1
    end = int(otu_range[1])-1
    if is_show_log: print(f"Reading {otu}.gfa file...\n")        
    if not s:
        if is_show_log: print("Skipping save graph information")
        if os.path.exists(f"{pantax_db}/species_gfa"):
            paths, nodes_len_npy  = read_gfa(f"{pantax_db}/species_gfa/{otu}.gfa")
        elif os.path.exists(f"{pantax_db}/species_graph_info"):
            paths, nodes_len_npy = read_h5py_file(f"{pantax_db}/species_graph_info/{otu}.h5")
    else:
        if os.path.exists(f"species_graph_info/{otu}.h5"):
            paths, nodes_len_npy = read_h5py_file(f"species_graph_info/{otu}.h5")
        else:
            paths, nodes_len_npy  = read_gfa(f"{pantax_db}/species_gfa/{otu}.gfa")
            write_h5py_file(nodes_len_npy, paths, f"species_graph_info/{otu}.h5")
    unique_trio_nodes, unique_trio_nodes_len, hap2unique_trio_nodes_m = trio_nodes_info(paths, nodes_len_npy) 
    nodes_len = {node:node_len for node, node_len in enumerate(nodes_len_npy)}
    node_abundances, unique_trio_node_abundances = get_node_abundances_dev(otu, read_group_data[otu], nodes_len, unique_trio_nodes, unique_trio_nodes_len, start) 
    non_zero_count = sum(1 for elem in node_abundances if elem != 0)
    if is_show_log: print(f"{otu} species node abundance > 0 number:{non_zero_count}\n")
    otu_paths = list(paths.values())
    haps_id = list(paths.keys())
    nvert = end-start+1      
    otu_abundance_list = [x if x > min_depth else 0 for x in node_abundances]
    if non_zero_count == 0 or max(otu_abundance_list) == 0 :
        return [(otu, hap_id, 0) for hap_id in haps_id]
    # Now solve the minimization problem  
    x, objVal = optimize(otu_abundance_list, nvert, otu_paths, minimization_min_cov, min_cov, 
                         unique_trio_nodes_fraction, unique_trio_nodes_mean_count_fraction,
                         threads, hap2unique_trio_nodes_m, unique_trio_node_abundances)
    return [(otu, hap_id, x[i]) for i, hap_id in enumerate(haps_id)]

def read_cls(otu_range_file: str, species_abundance_file: str, min_species_abundance: int = 1e-04) -> Dict[str, List[str]]:
    """
    Obtain species filtered for abundance, and then obtain their location information in the reference pangenome.
    """
    species_abundance = pd.read_csv(species_abundance_file, sep="\t", dtype={0:str, 1:float, 2:float})
    species_abundance.columns = ["species_taxid", "abundance", "coverage"]
    species_abundance = species_abundance[species_abundance["abundance"] > min_species_abundance]
    otu_list = species_abundance["species_taxid"].tolist()
    # otu_list = ["1589", "630"]
    # print(len(otu_list))
    otu_to_range: Dict[str, List[str]] = {} # [key] = species_taxid, [value] = [start,end] -> the start and end position in reference_pangenome.gfa
    with open(otu_range_file, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")
            otu = tokens[0]
            if otu in otu_list:
                otu_to_range[otu] = [tokens[1], tokens[2]]
    return otu_to_range

# @timeit
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
                haplotype_id=None
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
                if reverse_flag:
                    path = list(reversed(path))
                # Multiple chromosomes from the same genome merge into one path. 
                try:
                    paths[haplotype_id] = path 
                except:
                    paths[haplotype_id].extend(path)
                    paths[haplotype_id] = list(set(paths[haplotype_id]))               
    return paths, np.array(nodes_len_list)

def abundance_cal(otu_cov: List[Tuple[str, str, Union[int, float]]], genomes_info_file: str) -> None:
    otu_cov_df = pd.DataFrame(otu_cov, columns=["species_taxid", "hap_id", "predicted_coverage"])
    sum_coverage = otu_cov_df["predicted_coverage"].sum()
    otu_cov_df["predicted_abundance"] = otu_cov_df["predicted_coverage"] / sum_coverage
    dtypes = {"genome_ID": str, "strain_taxid": str}
    genomes_info = pd.read_csv(genomes_info_file, sep="\t",usecols=[0,1],dtype=dtypes)
    genomes_info["hap_id"] = genomes_info["genome_ID"].str.split("_").str[:2].str.join("_")
    new_otu_cov_df = pd.merge(otu_cov_df, genomes_info, on="hap_id", how="left")
    new_otu_cov_df = new_otu_cov_df.drop("hap_id", axis=1)
    new_order = ["species_taxid", "strain_taxid", "genome_ID", "predicted_coverage", "predicted_abundance"]
    new_otu_cov_df = new_otu_cov_df[new_order]
    new_otu_cov_df = new_otu_cov_df.sort_values(by="predicted_abundance", ascending=False)
    new_otu_cov_df = new_otu_cov_df[new_otu_cov_df["predicted_abundance"] != 0]
    new_otu_cov_df.to_csv("strain_abundance.txt", sep="\t", index=False)

def abundace_constraint(otu_mapping_list: List[Tuple[str, str, Union[int, float]]], species_abundance_file: str) -> List[Tuple[str, str, Union[int, float]]]:
    species_abundance_df = pd.read_csv(species_abundance_file, sep="\t")
    species_abundance_df.columns = ["species_taxid", "abundance", "avg_coverage"]
    species_abundance_df["species_taxid"] = species_abundance_df["species_taxid"].astype(str)
    species_taxid = otu_mapping_list[0][0]
    strains: List[str] = [] 
    strains_coverage: List[Union[int, float]] = []
    for mapping in otu_mapping_list:
        strains.append(mapping[1])
        strains_coverage.append(mapping[2])
    species_coverage = species_abundance_df[species_abundance_df["species_taxid"] == species_taxid]["avg_coverage"].tolist()[0]
    if max(strains_coverage) > 1.05*species_coverage:
        factor = species_coverage/sum(strains_coverage)
        strains_coverage = np.array(strains_coverage)*factor
        new_otu_mapping_list: List[Tuple[str, str, Union[int, float]]] = []
        for i in range(len(strains)):
            new_otu_mapping_list.append((species_taxid, strains[i], strains_coverage[i]))
        return new_otu_mapping_list
    else:
        return otu_mapping_list

# @timeit
def trio_nodes_info(paths: Dict[str, List[int]], nodes_len_npy: np.ndarray) -> Tuple[Dict[Tuple[int, int, int], int], List[int], np.ndarray]:
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
    hap2unique_trio_nodes_m = np.stack(hap2unique_trio_nodes_v)
    # trio nodes length
    unique_trio_nodes_len = [sum(nodes_len_npy[idx] for idx in trio_nodes[i]) for i in unique_trio_nodes_idx]
    unique_trio_nodes = {trio_nodes[i]: idx for idx, i in enumerate(unique_trio_nodes_idx)}
    return unique_trio_nodes, unique_trio_nodes_len, hap2unique_trio_nodes_m

# @timeit
def get_node_abundances_dev(otu: str, reads_info: List[List[str]], nodes_len: Dict[int, int], trio_nodes: Dict[Tuple[int, int, int], int], 
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
    for i in range(len(nodes_len)):
        bases_per_node[i] = 0
    trio_nodes_bases_count = {}
    for i in range(len(trio_nodes)):
        trio_nodes_bases_count[i] = 0
    if is_show_log: ("Processing alignments...")
    for read_info in reads_info:
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
        if start_node == end_node and len(read_nodes) == 1:
            read_nodes_len[start_node] += target_len
            bases_per_node[start_node] += target_len
        else:
            i = 1
            for node, node_len in read_nodes_len_in_graph.items():
                if i == 1:
                    if node_len < read_start:
                        print(f"OTU: {otu}, read_id: {read_info[0]}, node_len: {node_len}, read_start:{read_start}")
                        print(read_nodes)
                        print(read_nodes_len_in_graph)
                        print(read_nodes_len)                        
                    assert node_len >= read_start
                    node_aln_len = node_len - read_start
                elif i == len(read_nodes_len_in_graph):
                    if target_len < seen: target_len = seen
                    node_aln_len = target_len - seen
                else:
                    node_aln_len = node_len
                seen += node_aln_len
                read_nodes_len[node] += node_aln_len
                bases_per_node[node] += node_aln_len
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
    if is_show_log: print("Computing node abundance rates...")
    for node, node_len in nodes_len.items():
        aligned_len = bases_per_node[node]
        if node_len > 0:
            node_abundance = aligned_len / node_len
        else:
            raise ZeroDivisionError(f"Node length 0 for node {node}")
        node_abundance_list.append(node_abundance)
    trio_node_abundance_list: List[int] = []
    if is_show_log: print("Computing trio node abundance rates...")
    for i, trio_node_len in enumerate(trio_nodes_len):
        aligned_len = trio_nodes_bases_count[i]
        if trio_node_len > 0:
            node_abundance = aligned_len / trio_node_len
        else:
            raise ZeroDivisionError(f"Trio node length 0 for node {node}")
        trio_node_abundance_list.append(node_abundance)
    return node_abundance_list, trio_node_abundance_list

# @timeit
def optimize(a: List[int], nvert: int, paths: List[List[int]], min_cov: int, min_cov_final: float, 
             unique_trio_nodes_fraction: float, unique_trio_nodes_mean_count_fraction: float, 
             threads: int, hap2trio_nodes_m: np.ndarray, trio_node_abundances: List[int]) -> Tuple[List[Union[int, float]], Union[int, float]]:
    """
    Defines Gurobi minimization problem and then applies the LP solver.
    Returns the solution values, the corresponding objective value, and a matrix
    counting pairwise differences between haplotypes.
    """
    t1_start = time.perf_counter()
    t2_start = time.process_time()
    origin_paths_len = len(paths)
    max_strains = origin_paths_len
    trio_node_abundances = np.array(trio_node_abundances)
    row_sums = np.sum(hap2trio_nodes_m, axis=1)
    if origin_paths_len != 1:
        possible_strains_frequencies_mean = []
        possible_strains_idx = []
        for idx in range(len(paths)):
            selected_indices = np.where((hap2trio_nodes_m[:, idx] == 1) & (row_sums == 1))[0]
            if len(selected_indices) == 0: continue
            frequencies = np.zeros(len(trio_node_abundances))
            frequencies[selected_indices] = trio_node_abundances[selected_indices]
            count_greater_than_zero = np.sum(frequencies > 0)
            if is_show_log: print(f"\t\t{count_greater_than_zero}/{len(selected_indices)}")
            if count_greater_than_zero/len(selected_indices) < unique_trio_nodes_fraction: continue
            possible_strains_idx.append(idx)
            non_zero_frequencies = frequencies[frequencies > 0]
            frequencies_mean = np.mean(non_zero_frequencies) if len(non_zero_frequencies) > 0 else 0
            possible_strains_frequencies_mean.append(frequencies_mean)
        if is_show_log: print("\t\tFisrt filter #strains / #paths = {} / {}".format(len(possible_strains_idx), origin_paths_len))
        paths = [paths[idx] for idx in possible_strains_idx]
    else:
        possible_strains_idx = [0]
    # t_stop = time.perf_counter()
    # print("     [PAO first filter based on fr] - Elapsed time: {:.1f} seconds".format(t_stop-t1_start))

    # absolute error
    m = Model('lp')
    obj = LinExpr()

    npaths = len(paths)          #The number of feasible paths
    P = np.zeros((nvert,npaths))
    x = m.addVars(list(range(npaths)), lb=0, ub=1.05*max(a), vtype=GRB.CONTINUOUS, name='x')
    X = [x[i] for i in range(npaths)]
    X = np.array([X]).reshape(npaths,1)    #Set x in an array for multiplication

    # If objective involves absolute values, add extra variables
    y = m.addVars(list(range(nvert)), lb=0, vtype=GRB.CONTINUOUS, name='y')
    # add indicator variables to count strains
    if max_strains > 0:
        # add indicator variables for counting strains
        x_binary = m.addVars(list(range(npaths)), vtype=GRB.BINARY, name='strain_indicator')
        for i in range(npaths):
            max_x = 2*max(a)
            m.addConstr(x_binary[i] >= (x[i]-min_cov)/max_x)
        # define total strain count
        sum_x_binary = LinExpr()
        for i in range(npaths):
            sum_x_binary += x_binary[i]
        # bound the number of strains
        m.addConstr(sum_x_binary <= max_strains)

    # Store paths in P: p_ij = 1 if node i contains path j
    # print('\nSave for every node which paths are passing through:')
    # for i in tqdm(range(npaths)):
    for i in range(npaths):
        for v in paths[i]:
            P[int(v),i] = 1
    npaths = len(paths)
    del paths
    m.update()

    # Define the objective function
    # print('\nDefine the objective function:')
    n_eval = 0
    # for v in tqdm(range(nvert)):
    for v in range(nvert):
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
    m.Params.Threads = threads
    m.Params.NumericFocus = 0
    m.Params.PoolSearchMode = 0
    m.Params.PoolSolutions = 10
    m.Params.Method = 4

    #Minimize the model for the given objective function and constraints
    if is_show_log: print("\t\t*** Phase 1 optimization***")
    m.optimize()

    # print("Number of solutions = {}".format(m.solcount))
    if m.status == GRB.Status.OPTIMAL:
        x_sol = []
        for v in m.getVars():
            if 'x' in v.varName:
                x_sol.append(v.x)
        if is_show_log: print(f"\t\tObjective value: {m.objVal}")
        objVal = m.objVal
        if is_show_log: print(f"\t\tfirst sol:{x_sol}\n")
        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_sol]
        nstrains = sum(selected_strains)
        if is_show_log: print("\t\tFirst optimization #strains / #paths = {} / {}".format(nstrains, npaths))

        if origin_paths_len == 1:
            if x_sol[0] < min_cov_final:
                x_sol[0] = 0
            return x_sol, objVal
        
        for idx, frequencies_mean in enumerate(possible_strains_frequencies_mean):
            if frequencies_mean == 0:
                selected_strains[idx] = 0
                continue
            f = abs(x_sol[idx]-frequencies_mean)/(x_sol[idx]+frequencies_mean)
            if is_show_log: print(f"\t\tidx:{idx}\tfrequencies_mean:{frequencies_mean}\txsol:{x_sol[idx]}\tf:{f}")
            if f > unique_trio_nodes_mean_count_fraction:
                selected_strains[idx] = 0
       
        nstrains = sum(selected_strains)
        if is_show_log: print("\t\tSecond filter #strains / #paths = {} / {}".format(nstrains, npaths))

        # run phase 2 optimization:
        if is_show_log: print("\t\t*** Phase 2 optimization***")
        m.reset()
        for i in range(npaths):
            if selected_strains[i] == 0:
                m.addConstr(x[i] == 0)
        m.optimize()

        x_final = []
        for v in m.getVars():
            if 'x' in v.varName:
                x_final.append(v.x)
        if is_show_log: print(f"\t\tObjective value: {m.objVal}")
        objVal = m.objVal

        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_final]
        nstrains = sum(selected_strains)
        if is_show_log: print("\t\tSecond optimization #strains / #paths = {} / {}".format(nstrains, npaths))

        sol = [0] * origin_paths_len
        i = 0
        for idx in range(len(sol)):
            if idx in possible_strains_idx:
                sol[idx] = x_final[i]
                i += 1

        t1_stop = time.perf_counter()
        t2_stop = time.process_time()
        if is_show_log: 
            print("\nStrain path optimize completed")
            print("Elapsed time: {:.1f} seconds".format(t1_stop-t1_start))
            print("CPU process time: {:.1f} seconds\n".format(t2_stop-t2_start))
        return sol, objVal

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


if __name__ == "__main__":
    sys.exit(main())
