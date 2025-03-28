#!/usr/bin/env python3


import sys,argparse
import pandas as pd
from toolkits import check_file_avail, extract_genome_name

usage = """
    Get the filtered genomes information from hierarchical clustering. The genomes information of all the genomes of the cluster 
    corresponding to these representative genomes is obtained according to the representative genomes obtained
    from hierarchical clustering.
"""


def main():
    parser = argparse.ArgumentParser(prog="python get_rep_cluster_genomes_info.py", usage=usage)
    parser.add_argument("-cls", "--hcls_file", dest="hcls_file", type=str, help="Hierarchical clustering result file.")
    parser.add_argument("-f", "--genomes_info", dest="genomes_info", type=str, help="Genomes information file.")
    parser.add_argument("-i", "--strain_abund", dest="strain_abund", type=str, help="Strain abundance file.")
    parser.add_argument("--gtdb", dest="gtdb", action="store_true", help="GTDB genome accession.")
    parser.add_argument("-o", "--out", dest="out", default="rep_cluster_genomes_info.txt", type=str, help="Output file path.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if args.gtdb:
        gtdb = True
    else:
        gtdb = False

    hcls_file = check_file_avail(args.hcls_file)
    genomes_info_file = check_file_avail(args.genomes_info)
    strain_abund_file = check_file_avail(args.strain_abund)

    strain_abund = pd.read_csv(strain_abund_file, sep="\t")
    represent_genomes = strain_abund["genome_ID"].tolist()
            
    genome_cluster = []
    with open(hcls_file, "r") as f:
        for line in f:
            tokens = line.strip().split("\t")
            genomes = tokens[2].split(",")
            genomes = [extract_genome_name(genome, gtdb) for genome in genomes]
            for genome in represent_genomes:
                if genome in genomes:
                    genome_cluster.extend(genomes)

    genomes_info = pd.read_csv(genomes_info_file, sep="\t")
    filtered_genomes_info = genomes_info[genomes_info["genome_ID"].isin(genome_cluster)]
    filtered_genomes_info.to_csv(f"{args.out}", sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())

