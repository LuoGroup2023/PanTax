#!/usr/bin/env python3
import sys, os, subprocess, argparse, re
import concurrent.futures
import numpy as np
import pandas as pd
from functools import partial
import networkx as nx
from itertools import combinations
from toolkits import Logger

usage = "Using graph-based clustering algorithms to reduce genome redundancy and obtain non-redundant genomes."

class GenomesCluster:

    def __init__(self, summary_file_path, genomes_info, genome_stat, database, output_cluster, gtdb, m, n, p, j):
        self.summary_file_path = summary_file_path
        self.genomes_info = genomes_info
        self.genome_stat = genome_stat
        self.database = database
        self.output_cluster = output_cluster
        self.gtdb_accession = gtdb
        self.m = m
        self.n = n
        self.p = p
        self.j = j
        self.genome_statics_data = self.genome_statics_data_read(genome_stat)
        if os.path.exists(self.summary_file_path):
            self.reference_or_respresentative_species_set = self.preprocess()
        else:
            self.reference_or_respresentative_species_set = {}
        

    def preprocess(self):
        """
        Preprocessing. Obtaining a representative genome or reference genome of each species.
        """
        high_quality_species_genome_data_set, refer_species_taxid = self.reference_or_respresentative_genome()
        reference_or_respresentative_species_set = self.reference_or_respresentative_genome_duplication_remove_and_restruct(high_quality_species_genome_data_set, species_taxid=refer_species_taxid)
        return reference_or_respresentative_species_set

    def genome_statics_data_read(self, genome_statics_file_path):
        """"read genome_statics_data"""
        genome_statics_data = pd.read_csv(genome_statics_file_path, usecols= [0, 3, 5], delimiter='\t', header=None)
        genome_statics_data.columns = ["bacteria", "sca_gap_length", "sca_N50"]
        sca_N50_threshold = 0
        sca_gap_length_threshold = 10000000
        genome_statics_data = genome_statics_data[genome_statics_data["sca_N50"] > sca_N50_threshold]
        genome_statics_data = genome_statics_data[genome_statics_data["sca_gap_length"] < sca_gap_length_threshold]
        return genome_statics_data

    def genome_cluster(self):
        """For every strain, if meet the conditions(eg. N50), added to species cluster"""
        genome_statics_data = self.genome_statics_data
        filter_genomes = genome_statics_data["bacteria"].tolist()
        if not self.gtdb_accession:
            filter_genomes = [re.sub(r"_genomic\.fna(\.gz)?$", "", file) for file in filter_genomes]
        else:
            filter_genomes = ["_".join(file.split("_")[:2]) for file in filter_genomes]
        filter_genomes_df = pd.DataFrame({"genome_ID": filter_genomes})
        provided_genomes = pd.read_csv(self.genomes_info, sep="\t")
        filter_genomes_info = pd.merge(filter_genomes_df, provided_genomes, on="genome_ID", how="left")
        filter_genomes_info["species_taxid"] = filter_genomes_info["species_taxid"].astype(str)
        grouped_genomes_info = filter_genomes_info.groupby("species_taxid")
        cluster_species = {}
        for species_taxid, group in grouped_genomes_info:
            grouped_genomes = group["id"].tolist()
            grouped_genomes = [os.path.basename(genome) for genome in grouped_genomes]
            if self.gtdb_accession: species_taxid = species_taxid.replace(" ", "_")
            cluster_species[species_taxid] = grouped_genomes
        return cluster_species

    def reference_or_respresentative_genome(self):
        """get reference genome(only one) or respresentative genomes(possibly many) for some species(not all species have)"""
        species_taxid = []
        genome_statics_data = self.genome_statics_data
        bacteria_genome_name = genome_statics_data["bacteria"].tolist()
        fn_map = {}
        with open(self.summary_file_path, "r") as f:
            f.readline()
            f.readline()
            for line in f:
                line = line.strip()
                tokens = line.split("\t")
                if (tokens[4] == "reference genome" or tokens[4] == "representative genome") and tokens[11] == "Complete Genome":
                    fn = os.path.basename(tokens[19]) + "_genomic.fna"
                    if fn not in bacteria_genome_name or fn in fn_map:
                        continue
                    fn_map[fn] = [tokens[6], tokens[7]]
                    species_taxid.append(tokens[6])
        return fn_map, species_taxid

    def reference_or_respresentative_genome_duplication_remove_and_restruct(self, high_quality_species_genome_data_set, species_taxid):
        """If some species has more than one respresentative genomes, choose a optimal genome(max N50) for the species. If reference genome, only one"""
        genome_statics_data = self.genome_statics_data
        species_taxid_duplication = [x for x in species_taxid if species_taxid.count(x) > 1]
        species_taxid_duplication = list(np.unique(species_taxid_duplication))
        species_taxid_delete_index = []
        genome_data_set_delete_key = []
        for taxid in species_taxid_duplication:
            taxid_index = list(np.where(np.array(species_taxid) == taxid)[0])
            bacteria_files = [file for file, species_data in high_quality_species_genome_data_set.items() if species_data[0] == taxid]
            matched_data_frames = []
            for bacteria_file in bacteria_files:
                matched_row = genome_statics_data[genome_statics_data["bacteria"] == bacteria_file]
                if not matched_row.empty:
                    matched_data_frames.append(matched_row)
            N50_max_row = max(matched_data_frames, key=lambda df: df["sca_N50"].iloc[0])
            opt_bacteria_file = N50_max_row["bacteria"].iloc[0]
            bacteria_files.remove(opt_bacteria_file)
            genome_data_set_delete_key.append(bacteria_files)
            opt_bacteria_file_index = {key: index for index, key in enumerate(high_quality_species_genome_data_set)}.get(opt_bacteria_file)
            taxid_index.remove(opt_bacteria_file_index)
            species_taxid_delete_index.append(taxid_index)
        for bacteria_files in genome_data_set_delete_key:
            for bacteria_file in bacteria_files:
                del high_quality_species_genome_data_set[bacteria_file]
        species_taxid_delete_index = [index for sublist in species_taxid_delete_index for index in sublist]
        species_taxid = [item for index, item in enumerate(species_taxid) if index not in species_taxid_delete_index]

        reference_or_respresentative_species_set = {}
        for key, value in high_quality_species_genome_data_set.items():
            taxid = value[0]
            reference_or_respresentative_species_set[taxid] = key
        return reference_or_respresentative_species_set

    def gernerate_map_file(self, species_taxid, species_cluster, reference_or_respresentative_species_set):    
        """Select genomes used to cluster."""
        # if not os.path.exists(f"{self.output_cluster}/{species_taxid}/{species_taxid}_result.txt"):
        #     print(f"{self.output_cluster}/{species_taxid}")
        #     subprocess.run(f"rm -rf {self.output_cluster}/{species_taxid}", shell=True)
        if not os.path.exists(f"{self.output_cluster}/{species_taxid}"):
            os.mkdir(f"{self.output_cluster}/{species_taxid}")
            species_cluster = list(set(species_cluster))
            if self.m == -1: self.m = len(species_cluster)
            if len(species_cluster) >= self.m:
                genome_statics_data = self.genome_statics_data
                strain_subset = genome_statics_data[genome_statics_data["bacteria"].isin(species_cluster)].reset_index()
                strain_subset.sort_values(by='sca_N50', ascending=False, inplace=True)
                strain_subset100 = strain_subset.iloc[0:self.m, [1]]["bacteria"].tolist()
                if species_taxid in reference_or_respresentative_species_set and reference_or_respresentative_species_set[species_taxid] not in strain_subset:
                    strain_subset100.append(reference_or_respresentative_species_set[species_taxid])
                species_cluster = strain_subset100
            species_cluster = [os.path.join(self.database, file) for file in species_cluster]
            with open(f"{self.output_cluster}/{species_taxid}/{species_taxid}.txt", "w") as f:
                f.write("\n".join(species_cluster) + "\n")
    
    def fastANI_compute(self, species_taxid):
        """fastANI used to compute ANI"""
        input_file = f"{self.output_cluster}/{species_taxid}/{species_taxid}.txt"
        result_file = f"{self.output_cluster}/{species_taxid}/{species_taxid}_result.txt"
        if not os.path.exists(result_file):
            command = f"fastANI --ql {input_file} --rl {input_file} -o {result_file} --threads {self.j} > /dev/null 2>&1"
            subprocess.run(command, shell=True)

    def parallel_compute(self, species_taxid, species_clusters):
        partial_process = partial(self.gernerate_map_file, reference_or_respresentative_species_set = self.reference_or_respresentative_species_set)
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.j) as executor:
            futures = {key: executor.submit(partial_process, key, paths) for key, paths in species_clusters.items()}
            concurrent.futures.wait(futures.values())
        # for taxid in species_taxid:
        #     fastANI_compute(taxid)
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.p) as executor:
            executor.map(self.fastANI_compute, species_taxid)

    def generate_softlink_and_strain_file(self, genome_data_set, rep_to_cluster, species_taxid):
        """The genomes of each species used to construct the pangenome are saved in the diff_strains file"""
        source_dir = self.database
        failure_file = f"{self.output_cluster}/{species_taxid}/failure_file.txt"
        if not os.path.exists(failure_file):
            with open(failure_file, "w") as failure_f:
                if self.n == -1: self.n = len(genome_data_set)
                for bacteria_file in genome_data_set[:self.n]:
                    bacteria_file = os.path.basename(bacteria_file)
                    source_path = os.path.join(source_dir, bacteria_file)
                    target_path = os.path.join(self.output_cluster, species_taxid, bacteria_file)
                    try:
                        os.symlink(source_path, target_path)
                    except Exception as e:
                        failure_f.write(bacteria_file + "\n")
        diff_strains = f"{self.output_cluster}/{species_taxid}/{species_taxid}_species_with_diff_strains.txt"
        if not os.path.exists(diff_strains):
            with open(diff_strains, "w") as f:
                for bacteria_file in genome_data_set:
                    bacteria_file = os.path.basename(bacteria_file)
                    bacteria_file = os.path.join(self.output_cluster, species_taxid, bacteria_file)
                    f.write(bacteria_file + "\n")
        rep_to_cluster_file = f"{self.output_cluster}/{species_taxid}/{species_taxid}_rep_to_cluster.txt"
        if not os.path.exists(rep_to_cluster_file):
            with open(rep_to_cluster_file, "w") as f:
                for rep, cluster in rep_to_cluster.items():
                    cluster = ",".join(cluster)
                    f.write(f"{rep}\t{cluster}\n")

    def process_among_clusters(self, species_taxid, reference_or_respresentative_species_set, genome_statics_data):
        """Genome cluster based connected graph. Select genome whose ANI between 95 to 99.9."""
        if os.path.exists(f"{self.output_cluster}/{species_taxid}/{species_taxid}_result.txt") and not os.path.exists(f"{self.output_cluster}/{species_taxid}/failure_file.txt"):
            ANI_data = pd.read_csv(f"{self.output_cluster}/{species_taxid}/{species_taxid}_result.txt", usecols= [0, 1, 2], delimiter='\t', header=None)
            ANI_data.columns = ["query", "refer", "ANI"]
            result_strains_set = []
            rep_to_cluster = {}
            if len(ANI_data) == 1 and (ANI_data["query"] == ANI_data["refer"]).item():
                result_strains_set = ANI_data["query"].tolist()
                rep_to_cluster[result_strains_set[0]] = result_strains_set[0]
            else:
                # due to computing, some compared with themselves appear to have 99.9999 - very close to 100 
                ANI_data_origin = ANI_data
                ANI_data = ANI_data.drop(index=ANI_data[ANI_data["query"] == ANI_data["refer"]].index)
                ANI_data = ANI_data[ANI_data["ANI"] >= 95].copy()
                if not ANI_data.empty:
                    all_file_nodes_in_species_cluster = list(set(ANI_data['query'].tolist() + ANI_data['refer'].tolist()))
                    # ANI_data['qrPair'] = ANI_data.apply(lambda row: tuple(sorted([row['query'], row['refer']])), axis=1)
                    # if using parallel, pandas apply function may cause some problems, it can't process some data, but no error returns
                    qrPairs = []
                    for _ , row in ANI_data.iterrows():
                        qrPair = tuple(sorted([row['query'], row['refer']]))
                        qrPairs.append(qrPair)
                    ANI_data['qrPair'] = qrPairs
                    ANI_data = ANI_data.groupby('qrPair')['ANI'].max().reset_index()
                    # ANI_data_ge99 = ANI_data[(ANI_data["ANI"] >= 99.9) & (ANI_data["ANI"] <= 100)].copy()
                    ANI_data_ge99 = ANI_data[ANI_data["ANI"] >= 99.9].copy()
                    edges = ANI_data_ge99["qrPair"].tolist()
                    
                    G = nx.Graph()
                    G.add_nodes_from(all_file_nodes_in_species_cluster)
                    G.add_edges_from(edges)
                    connected_components_set = nx.connected_components(G)
                    if species_taxid in reference_or_respresentative_species_set:
                        refer_node = os.path.join(self.database, reference_or_respresentative_species_set[species_taxid])
                        # print(refer_node)
                    else:
                        refer_node = None
                    for connected_component in sorted(connected_components_set, key=len, reverse=True):
                        # if len(connected_component) == 1:
                        #     continue
                        if refer_node in connected_component:
                            result_strains_set.append(refer_node)
                        else:
                            list_strains = [os.path.basename(i) for i in list(connected_component)]
                            subset_strains = genome_statics_data[genome_statics_data["bacteria"].isin(list_strains)]
                            # max_value = subset_strains["sca_N50"].max()
                            # single_component_result = subset_strains[subset_strains["sca_N50"] == max_value]["bacteria"].values[0]
                            # result_strains_set.append(single_component_result)
                            if not subset_strains.empty:
                                best_strain = subset_strains.loc[subset_strains["sca_N50"].idxmax(), "bacteria"]
                                result_strains_set.append(best_strain)
                                rep_to_cluster[best_strain] = list_strains
                    result_strains_set = [os.path.join(self.database, os.path.basename(strain)) for strain in result_strains_set]            
                    sorted_pairs = [tuple(sorted(pair)) for pair in combinations(result_strains_set, 2)]
                    qrPairs = ANI_data["qrPair"].tolist()
                    filtered_pair = [pair for pair in sorted_pairs if pair in qrPairs]
                    nodes = result_strains_set
                    edges = filtered_pair
                    G = nx.Graph()
                    G.add_nodes_from(nodes)
                    G.add_edges_from(edges)
                    cliques = nx.find_cliques(G)
                    cliques = sorted(cliques, key=len, reverse=True)
                    result_strains_set = cliques[0]
                    for clique in cliques:
                        if refer_node in clique:
                            result_strains_set = clique
                            break
            # some special condition, the ANI of all strains(more than two strains) are lower than 95, it will filter 
            if len(result_strains_set) == 0:
                if species_taxid in reference_or_respresentative_species_set:
                    refer_node = os.path.join(self.database, reference_or_respresentative_species_set[species_taxid])
                else:
                    refer_node = None
                list_strains = list(set(ANI_data_origin["query"].tolist()))
                if refer_node in list_strains:
                    result_strains_set.append(refer_node)
                else:
                    list_strains = [os.path.basename(i) for i in list(list_strains)]
                    subset_strains = genome_statics_data[genome_statics_data["bacteria"].isin(list_strains)]
                    max_value = subset_strains["sca_N50"].max()
                    better_strain = subset_strains[subset_strains["sca_N50"] == max_value]["bacteria"].values[0]
                    result_strains_set.append(better_strain)  
                rep_to_cluster[result_strains_set[0]] = result_strains_set[0]    
            self.generate_softlink_and_strain_file(result_strains_set, rep_to_cluster, species_taxid)  

    def parallel_cls(self, species_taxid):        
        partial_process = partial(self.process_among_clusters, reference_or_respresentative_species_set = self.reference_or_respresentative_species_set, genome_statics_data = self.genome_statics_data)
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(partial_process, species_taxid)  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="python genomes_cluster.py", description=usage)
    parser.add_argument("-m", default=100, type=int, help="Max genomes number used to cluster every species(default:100)")
    parser.add_argument("-n", default=10, type=int, help="Max genomes number used to build pangenome every species(default:10)")
    parser.add_argument("-s", "--summary_file", dest="summary_file_path", default="assembly_summary_bacteria.txt", type=str, help="Assembly summary file path")
    parser.add_argument("-i", "--genomes_info", dest="genomes_info", default="genomes_info_provided_origin.txt", type=str, help="Provided genomes information. Same with custom option in the last step")
    parser.add_argument("-gs", "--genome_stat", dest="genome_stat", default="genome_stat.txt", type=str, help="Genome statistics file path")
    parser.add_argument("-d", "--database", dest="database", default="reference_genomes_database", type=str, help="All genomes database path(absolute path better)")
    parser.add_argument("-o", "--output_cluster", dest="output_cluster", default="output_cluster", type=str, help="Output genomes cluster path")
    parser.add_argument("-p", default=2, type=int, help="Number of parallel processes used for ANI calculation(default:2)")
    parser.add_argument("-j", default=32, type=int, help="Number of parallel processes used for fastANI(default:32)")
    parser.add_argument("--gtdb", dest="gtdb", action="store_true", help="GTDB genome accession.")
    parser.add_argument("--compute", dest="compute", action="store_true", help="ANI calculation.")
    parser.add_argument("--cluster", dest="cluster", action="store_true", help="Genomes cluster based on ANI.")
    parser.add_argument("--test", dest="test", action="store_true", help="Test.")
    args = parser.parse_args()
    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    log = Logger()
    args.database = os.path.abspath(args.database)
    genomes_cluster = GenomesCluster(args.summary_file_path, args.genomes_info, args.genome_stat, args.database, args.output_cluster, args.gtdb, args.m, args.n, args.p, args.j)
    species_clusters = genomes_cluster.genome_cluster()
    log.logger.info(f"Cluster species number:{len(species_clusters)}")
    species_taxid = list(species_clusters.keys())
    if args.test:
        log.logger.info("Executing test compute steps")
        if not os.path.exists(args.output_cluster):
            os.mkdir(args.output_cluster)
        first_10_keys = list(species_clusters.keys())[:10]
        first_10_values = [species_clusters[key] for key in first_10_keys]
        small_dict = {key: value for key, value in zip(first_10_keys, first_10_values)}
        genomes_cluster.parallel_compute(first_10_keys, small_dict)
        log.logger.info("Executing test cluster steps")
        genomes_cluster.parallel_cls(first_10_keys)
        sys.exit(0)
    if args.compute:
        log.logger.info("Executing compute steps")
        if not os.path.exists(args.output_cluster):
            os.mkdir(args.output_cluster)
        genomes_cluster.parallel_compute(species_taxid, species_clusters)
    else:
        log.logger.info("Warning: skipping computing steps")
    if args.cluster:
        log.logger.info("Executing cluster steps")
        genomes_cluster.parallel_cls(species_taxid)
    else:
        log.logger.info("Warning: skipping cluster steps")
    log.logger.info("Completely.")

    