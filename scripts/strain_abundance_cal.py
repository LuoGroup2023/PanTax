#!/usr/bin/env python3
import sys, os, re, argparse, time, h5py, random
import numpy as np
import pandas as pd
from gurobipy import *
from tqdm import tqdm # progress tracker
from functools import partial
from aln_json_process import read_group
import concurrent.futures 
from toolkits import timeit, Logger
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

usage = "Compute strain abundance"

def main():
    parser = argparse.ArgumentParser(prog="strain_abundance_cal.py", description=usage)
    parser.add_argument("db", type=str, help="pantax database directory")
    parser.add_argument("genomes_info", type=str, help="Genomes information file")
    parser.add_argument("read_cls", type=str, help="Reads classification file")
    parser.add_argument("aln_file", type=str, help="Reads alignment file(json format)")
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
    parser.add_argument("--sample", dest="sample", type=int, default=0, help="Sampling nodes are used for small model testing")
    args = parser.parse_args()

    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()
    log = Logger().logger
    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()
    global sample
    sample = args.sample
    species_graph_info = os.path.join(args.db, "species_graph_info")
    if not os.path.exists(f"{species_graph_info}") and args.s:
        os.mkdir("species_graph_info")
    minimization_min_cov = 0
    otu_to_range = read_cls(args.otu_range_file, args.species_abundance_file, args.min_species_abundance)
    otu_cov = []
    if len(otu_to_range) == 1:
        args.gurobi_threads = max(args.threads, args.gurobi_threads)
        args.threads = 1
    global read_group_data, is_show_log
    is_show_log = args.is_show_log
    read_group_data = read_group(args.read_cls, args.aln_file) 
    log.info("Parallel estimation of abundance for each strain of every species(this step takes a long time)...")
    # for key, value in otu_to_range.items():
    #     result = parallel_optimize_otu(key, value, pantax_db=args.db, aln_file=args.aln_file, min_depth=args.min_depth, reduce_obj=args.reduce_obj, 
    #                           minimization_min_cov=minimization_min_cov, min_cov=args.min_cov, unique_trio_nodes_fraction=args.fr, unique_trio_nodes_mean_count_fraction=args.fc,
    #                           s=args.s, threads=args.gurobi_threads)
    #     result = abundace_constraint(result, species_abundance_file = args.species_abundance_file)
    #     otu_cov.extend(result)
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

    log.info("All strains optimize completed")
    log.info("Elapsed time: {:.1f} seconds".format(t3_stop-global_t1_start))
    log.info("CPU process time: {:.1f} seconds".format(t4_stop-global_t2_start))
    print()
    return

################################################################################

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

def parallel_optimize_otu(otu, otu_range, pantax_db, min_depth, minimization_min_cov, min_cov, unique_trio_nodes_fraction, unique_trio_nodes_mean_count_fraction, s, threads):
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
    node_abundances, unique_trio_node_abundances = get_node_abundances(read_group_data[otu], list(nodes_len_npy), unique_trio_nodes, unique_trio_nodes_len, start)  
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

def read_cls(otu_range_file, species_abundance_file, min_species_abundance=1e-04):
    """
    Obtain species filtered for abundance, and then obtain their location information in the reference pangenome.
    """
    species_abundance = pd.read_csv(species_abundance_file, sep="\t", dtype={0:str, 1:float, 2:float})
    species_abundance.columns = ["species_taxid", "abundance", "coverage"]
    species_abundance = species_abundance[species_abundance["abundance"] > min_species_abundance]
    otu_list = species_abundance["species_taxid"].tolist()
    # print(len(otu_list))
    otu_to_range = {}
    with open(otu_range_file, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")
            otu = tokens[0]
            if otu in otu_list:
                otu_to_range[otu] = [tokens[1], tokens[2]]
    return otu_to_range

def read_gfa(graph_name, previous = 0):
    """
    Reads a graph from a GFA-file and returns graph in gt-format.
    """
    nodes_len_list = []
    paths = {}
    with open(graph_name, 'r') as f:
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

def abundance_cal(otu_cov, genomes_info_file):
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

def abundace_constraint(otu_mapping_list, species_abundance_file):
    species_abundance_df = pd.read_csv(species_abundance_file, sep="\t")
    species_abundance_df.columns = ["species_taxid", "abundance", "avg_coverage"]
    species_abundance_df["species_taxid"] = species_abundance_df["species_taxid"].astype(str)
    species_taxid = otu_mapping_list[0][0]
    strains = []
    strains_coverage = []
    for mapping in otu_mapping_list:
        strains.append(mapping[1])
        strains_coverage.append(mapping[2])
    species_coverage = species_abundance_df[species_abundance_df["species_taxid"] == species_taxid]["avg_coverage"].tolist()[0]
    if max(strains_coverage) > 1.05*species_coverage:
        factor = species_coverage/sum(strains_coverage)
        strains_coverage = np.array(strains_coverage)*factor
        new_otu_mapping_list = []
        for i in range(len(strains)):
            new_otu_mapping_list.append((species_taxid, strains[i], strains_coverage[i]))
        return new_otu_mapping_list
    else:
        return otu_mapping_list

def trio_nodes_info(paths, nodes_len_npy):
    trio_nodes = []
    hap_trio_paths = {}
    # chop trio nodes for all paths
    for hap, path in paths.items():
        trio_path = [(path[i], path[i+1], path[i+2]) for i in range(len(path)-2)]
        hap_trio_paths[hap] = trio_path
        trio_nodes.extend(trio_path)
    # get unique trio nodes
    trio_nodes = list(set(trio_nodes))
    # save trio nodes whether it exists in every path
    trio_nodes_dict = {}
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

    hap2unique_trio_nodes_v = []
    unique_trio_nodes_idx = []
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

def get_node_abundances(reads_info, nodes_len, trio_nodes, trio_nodes_len, start):
    bases_per_node = {} # map node IDs to read alignments
    for i in range(len(nodes_len)):
        bases_per_node[i] = 0
    trio_nodes_bases_count = {}
    for i in range(len(trio_nodes)):
        trio_nodes_bases_count[i] = 0
    if is_show_log: print("Processing alignments...")
    for read_info in reads_info:
        read_nodes = []
        read_nodes_len = {}
        for node_id, aln_len in read_info.items():
            origin_node_id = node_id-1-start
            read_nodes.append(origin_node_id)
            try:
                bases_per_node[origin_node_id] += aln_len
            except:
                print(f"node_id:{node_id}, start:{start}")
                continue
                # sys.exit(1)
            read_nodes_len[origin_node_id] = aln_len
        if len(read_nodes) < 3:
            continue
        read_trio_nodes = [(read_nodes[i], read_nodes[i+1], read_nodes[i+2]) for i in range(len(read_nodes)-2)]
        read_trio_nodes_len = [sum(read_nodes_len[idx] for idx in trio_node) for trio_node in read_trio_nodes]
        for i, read_trio_node in enumerate(read_trio_nodes):
            idx = trio_nodes.get(read_trio_node)
            if not idx:
                read_trio_node_reversed = tuple(reversed(read_trio_node))
                idx = trio_nodes.get(read_trio_node_reversed)
                if not idx:
                    continue
            trio_nodes_bases_count[idx] += read_trio_nodes_len[i]
    node_abundance_list = []
    if is_show_log: print("Computing node abundance rates...")
    for node, node_len in enumerate(nodes_len):
        aligned_len = bases_per_node[node]
        if node_len > 0:
            node_abundance = aligned_len / node_len
        else:
            raise ZeroDivisionError(f"Node length 0 for node {node}")
        node_abundance_list.append(node_abundance)
    trio_node_abundance_list = []
    if is_show_log: print("Computing trio node abundance rates...")
    for i, trio_node_len in enumerate(trio_nodes_len):
        aligned_len = trio_nodes_bases_count[i]
        if trio_node_len > 0:
            node_abundance = aligned_len / trio_node_len
        else:
            raise ZeroDivisionError(f"Trio node length 0 for node {node}")
        trio_node_abundance_list.append(node_abundance)
    return node_abundance_list, trio_node_abundance_list

def optimize(a, nvert, paths, min_cov, min_cov_final, 
             unique_trio_nodes_fraction, unique_trio_nodes_mean_count_fraction,
             threads, hap2trio_nodes_m, trio_node_abundances):
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
    same_path_flag = False
    if origin_paths_len != 1 and hap2trio_nodes_m.size:
        possible_strains_frequencies_mean = []
        possible_strains_idx = []
        row_sums = np.sum(hap2trio_nodes_m, axis=1)
        for idx in range(len(paths)):
            selected_indices = np.where((hap2trio_nodes_m[:, idx] == 1) & (row_sums == 1))[0]
            if len(selected_indices) == 0: continue
            frequencies = np.zeros(len(trio_node_abundances))
            frequencies[selected_indices] = trio_node_abundances[selected_indices]
            count_greater_than_zero = np.sum(frequencies > 0)
            if is_show_log: print(f"\t\tTrio node abundance > 0 ratio: {count_greater_than_zero}/{len(selected_indices)}")
            if count_greater_than_zero/len(selected_indices) < unique_trio_nodes_fraction: continue
            possible_strains_idx.append(idx)
            non_zero_frequencies = frequencies[frequencies > 0]
            frequencies_mean = np.mean(non_zero_frequencies) if len(non_zero_frequencies) > 0 else 0
            possible_strains_frequencies_mean.append(frequencies_mean)
        if is_show_log: print("\t\tFisrt filter #strains / #paths = {} / {}".format(len(possible_strains_idx), origin_paths_len))
        paths = [paths[idx] for idx in possible_strains_idx]
    # elif origin_paths_len == 1:
    #     possible_strains_idx = [0]
    elif origin_paths_len != 1 and hap2trio_nodes_m.size == 0:
        if all(x == paths[0] for x in paths):
            paths = [paths[0]]
            same_path_flag = True
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
    nvert_list = list(range(nvert))
    
    if sample:
        if len(nvert_list) > 500:
            sample_nodes = random.sample(nvert_list, 500)
            sample_nodes = sorted(sample_nodes)
            nvert_list = sample_nodes

    # If objective involves absolute values, add extra variables
    y = m.addVars(nvert_list, lb=0, vtype=GRB.CONTINUOUS, name='y')
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
        if sample:
            if v not in nvert_list: continue
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
    # print('\nObjective function ready, starting Gurobi optimization:\n')
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
        if is_show_log: print(f"\t\tFirst sol:{x_sol}\n")
        selected_strains = [1 if cov > min_cov_final else 0 for cov in x_sol]
        nstrains = sum(selected_strains)
        if is_show_log: print("\t\tFirst optimization #strains / #paths = {} / {}".format(nstrains, npaths))

        if origin_paths_len == 1:
            if x_sol[0] < min_cov_final:
                x_sol[0] = 0
            return x_sol, objVal
        elif origin_paths_len != 1 and hap2trio_nodes_m.size == 0:
            if same_path_flag:
                if x_sol[0] < min_cov_final:
                    x_sol[0] = 0  
                sol = [0] * origin_paths_len 
                sol[0] = x_sol[0]    
                x_sol = sol         
            else:
                for i in range(len(x_sol)):
                    if x_sol[i] < min_cov_final: x_sol[i] = 0
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
            print("CPU process time: {:.1f} seconds".format(t2_stop-t2_start))
        return(sol, objVal)

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
