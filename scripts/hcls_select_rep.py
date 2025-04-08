#!/usr/bin/env python3
import sys, argparse, warnings, subprocess
from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from toolkits import check_file_avail, check_dir_avail, is_file_non_empty, extract_genome_name

# This script is based on:
# https://github.com/liaoherui/StrainScan/blob/main/library/select_rep.py

usage = "Using hierarchical clustering algorithms to reduce genome redundancy and obtain non-redundant genomes."



def filter_matrix(genomes_info_file, matrix_file, out_dir):
    matrix_file_dir = Path(matrix_file).parent
    matrix_file_name = Path(matrix_file).name
    # filter_matrix_file_path = matrix_file_dir / f"filter_{matrix_file_name}"
    filter_matrix_file_path = Path(out_dir) / f"filter_{matrix_file_name}"

    if is_file_non_empty(filter_matrix_file_path):
        genomes_id = pd.read_csv(genomes_info_file, sep="\t")["id"].tolist()
        strings = []
        idx = []
        with open(matrix_file, "r") as f:
            f.readline()
            for i, line in enumerate(f):
                tokens = line.strip().split("\t")
                genome_id = tokens.pop(0)
                if genome_id in genomes_id:
                    if tokens:
                        ani = [tokens[_idx] for _idx in idx]
                        ani = "\t".join(ani)
                        string = f"{genome_id}\t{ani}"
                        strings.append(string)
                    idx.append(i)
        
        with open(filter_matrix_file_path, "w") as f:
            f.write(f"{len(strings)}\n")
            f.write("\n".join(strings) + "\n")
    else:
        warnings.warn(f"{filter_matrix_file_path} exists. Skip filtering matrix.")
    return filter_matrix_file_path

def rebuild_matrix(matrix_file, out_dir):
    matrix_file_dir = Path(matrix_file).parent
    matrix_file_name = Path(matrix_file).name
    # rebuild_matrix_file = matrix_file_dir / f"rebuild_{matrix_file_name}"
    rebuild_matrix_file = Path(out_dir) / f"rebuild_{matrix_file_name}"
    if is_file_non_empty(rebuild_matrix_file):
        return rebuild_matrix_file
    
    with open(matrix_file, "r") as f:
        num = int(f.readline().strip())
        nn = np.zeros((num, num))
        genomes = []
        i = 0
        for line in f:
            tokens = line.strip().split("\t")
            genomes.append(tokens[0])
            len_n = len(tokens)
            ani = tokens[1:len_n]
            if ani:
                for j in range(len(ani)):
                    nn[i][j] = 100 - float(ani[j])
                    nn[j][i] = 100 - float(ani[j])
            i += 1

    with open(rebuild_matrix_file, "w") as f:
        f.write("\t" + "\t".join(genomes) + "\n")
        for i in range(len(genomes)):
            row = list(nn[i])
            row = [str(a) for a in row]
            f.write(genomes[i] + "\t" + "\t".join(row) + "\n")
    return rebuild_matrix_file    

def hcls(matrix_file, method, cutoff, out_dir):
    with open(f'{out_dir}/tem_hcls.R', 'w+') as f:
        f.write(f"""
            x <- read.table("{matrix_file}", header=T, row.names=1)
            d <- as.dist(as(x, "matrix"))
            hc <- hclust(d, method="{method}")
            res <- sort(cutree(hc, h={cutoff}))
            res
        """)

    subprocess.run(f'Rscript {out_dir}/tem_hcls.R > {out_dir}/hcls_res.txt', shell=True)
    with open(f'{out_dir}/hcls_res.txt', 'r') as f:
        lines = f.readlines()

    cluster = {}
    cluster_num = None
    for idx, line in enumerate(lines[::-1]):
        line = line.strip()
        if idx % 2 == 0:
            tokens = line.split()
            if len(tokens) == 1:
                cluster_num = tokens[0]
                if cluster_num not in cluster:
                    cluster[cluster_num] = {}
            else:
                cluster_num = tokens
                for token in tokens:
                    if token not in cluster:
                        cluster[token] = {}  
        else:
            tokens = line.split()
            if len(tokens) == 1:
                genome = tokens[0]
                genome = extract_genome_name(genome, gtdb)
                cluster[cluster_num][genome] = ''
            else:
                for i, genome in enumerate(tokens):
                    cluster[cluster_num[i]][genome] = ''
    cls_file = f"{out_dir}/hclsMap_{100 - float(cutoff)}.txt"
    with open(cls_file, 'w+') as f:
        for cluster_number in cluster:
            cls_genomes = ','.join(cluster[cluster_number])
            f.write(f"{cluster_number}\t{len(cluster[cluster_number])}\t{cls_genomes}\n")

    Path(f"{out_dir}/tem_hcls.R").unlink()
    Path(f"{out_dir}/hcls_res.txt").unlink()
    return cls_file

def pick_rep(matrix_file, cls_file, cutoff, out_dir):
    """ 
    Select representative strains from the FastANI similarity matrix and classify them based on the hierarchical clustering file. 
    """
    
    # Read the distance matrix
    strain_index = {}  # Strain name -> Index
    strain_path = {}   # Strain name -> Full file path
    strain_dist = {}   # Strain name -> Distance array
    with open(matrix_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for idx, genome in enumerate(header):
            strain_name = extract_genome_name(genome, gtdb)  # Extract strain name (remove path and extension)
            strain_index[strain_name] = idx
            strain_path[strain_name] = genome

        for line in f:
            tokens = line.strip().split('\t')
            genome = tokens[0]
            strain_name = extract_genome_name(genome, gtdb)
            strain_dist[strain_name] = np.array(tokens[1:], dtype=float)
    

    # Read cluster file
    clusters = defaultdict(dict)  # Cluster ID -> {Representative strain file path}
    strain_to_rep = {}            # Strain -> Representative strain
    final_clusters = {}           # Representative strain -> Cluster ID

    cls_order = []  # Track cluster ID order
    cls_outfile = Path(out_dir) / f"hclsMap_{cutoff}_Rep.txt"

    with open(cls_file, 'r') as f, open(cls_outfile, 'w') as out:
        for line in f:
            cluster_id, count, *strains = line.strip().split('\t')
            cluster_id, count = int(cluster_id), int(count)
            cls_order.append(cluster_id)

            strain_list = strains[-1].split(',')
            
            if count == 1:
                rep = strain_list[0]
                final_clusters[rep] = cluster_id
                strain_to_rep[rep] = rep
                clusters[cluster_id][strain_path[rep]] = ''
                out.write(f"{line.strip()}\t0\n")

            elif count == 2:
                rep = strain_list[0]
                final_clusters[rep] = cluster_id
                clusters[cluster_id][strain_path[rep]] = ''
                for strain in strain_list:
                    strain_to_rep[strain] = rep
                out.write(f"{cluster_id}\t{count}\t{rep}\t0\n")

            else:
                # Calculate average distance for each strain and select the most suitable representative strain
                distances = {
                    strain: np.mean(strain_dist[strain][[strain_index[s] for s in strain_list if s != strain]])
                    for strain in strain_list
                }
                min_rep = min(distances, key=distances.get)
                min_val = distances[min_rep]
                min_dist, max_dist = np.min(strain_dist[min_rep]), np.max(strain_dist[min_rep])
                
                final_clusters[min_rep] = cluster_id
                clusters[cluster_id][strain_path[min_rep]] = ''
                for strain in strain_list:
                    strain_to_rep[strain] = min_rep

                out.write(f"{cluster_id}\t{count}\t{min_rep}\t{min_dist},{max_dist},{min_val}\n")

    # Process non-representative strains
    other_outfile = Path(out_dir) / f"Other_Strain_CN_{cutoff}.txt"
    cls_strains = defaultdict(dict)  # Cluster ID -> All strains in this cluster

    with open(other_outfile, 'w') as out:
        for strain in strain_index:
            if strain in final_clusters:
                cls_strains[final_clusters[strain]][strain] = ''
                continue
            
            # Calculate the closest representative strain
            closest_rep = min(
                final_clusters, key=lambda rep: strain_dist[strain][strain_index[rep]]
            )
            closest_rep_dist = strain_dist[strain][strain_index[closest_rep]]
            strain_rep_dist = strain_dist[strain][strain_index[strain_to_rep[strain]]]

            if closest_rep == strain_to_rep[strain]:
                cls_strains[final_clusters[strain_to_rep[strain]]][strain] = ''
                continue
            
            cls_strains[final_clusters[closest_rep]][strain] = ''
            # Strain_Name    Closest_Rep,Distance_To_Rep    Closest_Cluster_Rep,Distance_To_Cluster_Rep
            out.write(f"{strain}\t{strain_to_rep[strain]},{strain_rep_dist}\t{closest_rep},{closest_rep_dist}\n")

    # Regenerate cluster information
    recls_outfile = Path(out_dir) / f"hclsMap_{cutoff}_recls.txt"
    clusters_lvl2 = defaultdict(dict)

    with open(recls_outfile, 'w') as out:
        for cls_id in cls_order:
            out.write(f"{cls_id}\t{len(cls_strains[cls_id])}\t{','.join(cls_strains[cls_id])}\n")
            for strain in cls_strains[cls_id]:
                clusters_lvl2[cls_id][strain_path[strain]] = ''

    return dict(clusters), dict(clusters_lvl2)

def main():
    parser = argparse.ArgumentParser(prog="python hcls_select_rep.py", usage=usage)
    parser.add_argument("-m", "--matrix", dest="matrix", type=str, help="ANI matrix (only applicable to FastANI matrix for now).")
    parser.add_argument("-n", dest="genomes_num", type=str, help="Max genomes number used to build pangenome every species(default:10)")
    parser.add_argument("-f", "--genomes_info", dest="genomes_info", type=str, help="Genomes information file.")
    parser.add_argument("--gtdb", dest="gtdb", action="store_true", help="GTDB genome accession.")
    parser.add_argument("-hm", "--hcls_method", dest="hcls_method", default="complete", type=str, help="Hcls method. (single/complete)")
    parser.add_argument("-e", "--cutoff", dest="cutoff", default=99, type=float, help="Hcls cutoff. default: 99)")
    parser.add_argument("-mf", "--matrix_filter", dest="matrix_filter", action="store_true", help="Filter matrix if genomes in matrix not all in genomes information file.")
    parser.add_argument("-o", "--output_dir", dest="output_dir", default="hcls_res", type=str, help="Output directory.")
    parser.add_argument("-t", "--threads", dest="threads", default=12, type=int, help="Threads. default:12")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    global gtdb
    if args.gtdb:
        gtdb = True
    else:
        gtdb = False
    
    if args.matrix == "None": args.matrix = None
    if args.genomes_num == "None": args.genomes_num = None

    genomes_info_file = check_file_avail(args.genomes_info)
    out_dir = check_dir_avail(args.output_dir)
    if args.matrix:
        matrix_file = check_file_avail(args.matrix)
        if args.matrix_filter:
            matrix_file = filter_matrix(genomes_info_file, matrix_file, out_dir)
    else:
        subprocess.run(f"awk -F'\t' 'NR > 1 {{print $5}}' {genomes_info_file} > {out_dir}/genomes.txt", shell=True)
        subprocess.run(f"fastANI --rl {out_dir}/genomes.txt --ql {out_dir}/genomes.txt -o {out_dir}/ani_res --matrix -t {args.threads} >/dev/null 2>&1", shell=True)
        matrix_file = f"{out_dir}/ani_res.matrix"
    rebuild_matrix_file = rebuild_matrix(matrix_file, out_dir)
    cls_file = hcls(rebuild_matrix_file, args.hcls_method, str(100 - args.cutoff), out_dir)
    clusters, clusters_lvl2 = pick_rep(rebuild_matrix_file, cls_file, args.cutoff, out_dir)
    genomes = []
    for cluster_number in clusters:
        this_cluster_genomes = list(clusters[cluster_number])
        genomes.extend(this_cluster_genomes)
    if args.genomes_num:
        genomes = genomes[:int(args.genomes_num)]
    genomes_info = pd.read_csv(genomes_info_file, sep="\t")
    filtered_genomes_info = genomes_info[genomes_info["id"].isin(genomes)]
    filtered_genomes_info.to_csv(f"{out_dir}/hcls_filtered_genomes_info.txt", sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())
