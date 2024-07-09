#!/usr/bin/env python3
# process read classfication file with multiple processes and threads, but it need gaf(it only use in the script, so it need some time to produce)
# And it use pandas to read whole file once, so it need too many memory.
# 2G data need 10 min, may be less. Some step can be optimized.
import argparse, re
import numpy as np
import pandas as pd
from typing import List
import concurrent.futures

class ReadClassification:
    def __init__(self, *args):
        self.species_range_file, self.mapped_gaf_file = args
        self.pair_range = None
    
    def read_species_range_file(self) -> List[str]:
        """
        Obtain species and its range information.
        """
        species_range = pd.read_csv(self.species_range_file, sep="\t", header=None, dtype={0: str, 1: int, 2:int})
        species = species_range.iloc[:, 0].tolist()
        self.pair_range = list(zip(species_range.iloc[:, 1], species_range.iloc[:, 2]))
        return species

    def read_mapped_gaf_file(self):
        """
        Read GAF file with pandas. One read per line.
        """
        mapped_gaf = pd.read_csv(self.mapped_gaf_file, sep="\t", header=None, usecols=[0, 1, 5, 11])
        # reads = list(zip(mapped_gaf.iloc[:, 0], mapped_gaf.iloc[:, 1]))
        reads = [tuple(row) for row in mapped_gaf.values]
        return reads

    def sigle_read_classification(self, read):
        """
        Process single read.
        """
        read_name = read[0]
        read_nodes = np.array([int(match.group()) for match in re.finditer(r'\d+', read[2])])
        mapq = read[3]
        if len(read_nodes) > 1:
            read_nodes = np.array([np.min(read_nodes), np.max(read_nodes)])
        species_range = np.array(self.pair_range)
        is_in_range = np.any((species_range[:, 0] <= read_nodes[:, None]) & (read_nodes[:, None] <= species_range[:, 1]), axis=0)
        if np.size(np.where(is_in_range)):
            species_index = np.where(is_in_range)[0][0]
        else:
            species_index = None
            # print(read_name)
        return (read_name, species_index, read_nodes.tolist(), mapq, read[1])


    def thread_parallel_batch_read_classification(self, reads):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(self.sigle_read_classification, read) for read in reads]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        return results

    def process_parallel_batch_read_classification(self, reads, species, batch_size=100000):       
        with concurrent.futures.ProcessPoolExecutor(max_workers=64) as executor:
            results = []
            for i in range(0, len(reads), batch_size):
                batch_reads = reads[i:i+batch_size]
                future = executor.submit(self.thread_parallel_batch_read_classification, batch_reads)
                results.append(future)
            final_results = []
            for future in concurrent.futures.as_completed(results):
                final_results.extend(future.result())
        df1 = pd.DataFrame(final_results, columns=["read_name", "species_index", "read_nodes_range", "mapq", "read_len"])
        df2 = pd.DataFrame(species, columns=["species"]).reset_index()
        df2.columns = ["species_index", "species"]
        merge = pd.merge(df1, df2, on="species_index", how="left")
        merge = merge.drop(columns=["species_index", "read_nodes_range"])
        merge = merge.reindex(columns=["read_name", "mapq", "species", "read_len"])
        merge.to_csv("reads_classification.csv", index=False, header=False, sep="\t", quoting=0)
# awk -F '\t' '{print $3}' reads_classification.csv | sort | uniq -c | sort -nr > species_match.txt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="python read_classification.py")
    parser.add_argument("-s", "--species_range_file", dest="species_range_file", type=str, help="Species range file")
    parser.add_argument("-m", "--mapped_gaf_file", dest="mapped_gaf_file", type=str, help="Mapped gfa file")
    args = parser.parse_args()
    if not args.mapped_gaf_file:
        args.mapped_gaf_file = "gfa_mapped.gaf"
    read_cls = ReadClassification(args.species_range_file, args.mapped_gaf_file)
    species = read_cls.read_species_range_file()
    reads = read_cls.read_mapped_gaf_file()
    read_cls.process_parallel_batch_read_classification(reads, species)