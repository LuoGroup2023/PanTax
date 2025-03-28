#!/usr/bin/env python3
import sys, gzip
import pandas as pd
import concurrent.futures

def main():
    genomes_info_file = sys.argv[1]
    query_result_file = sys.argv[2]
    ani_threshold = sys.argv[3]
    out_dir = sys.argv[4]
    genomes_info = pd.read_csv(genomes_info_file, sep="\t")
    genome2id = dict(zip(genomes_info["genome_ID"], genomes_info["id"]))
    genome2seqid = get_genome_to_seqid(genome2id)
    sylph_result = pd.read_csv(query_result_file, sep="\t")
    filter_genomes(sylph_result, genome2seqid, genomes_info, ani_threshold, out_dir)


def filter_genomes(sylph_result, genome2seqid, genomes_info, ani_threshold, out_dir):
    sylph_result["seq_id"] = sylph_result["Contig_name"].str.split(" ").str[0]
    sylph_result = sylph_result[sylph_result["Adjusted_ANI"] >= int(ani_threshold)]
    merge_df = pd.merge(sylph_result, genome2seqid, on="seq_id", how="left")
    sylph_genomes = merge_df[["genome_ID"]]    
    filter_genomes_info = pd.merge(genomes_info, sylph_genomes, on="genome_ID")
    filter_genomes_info.to_csv(f"{out_dir}/filter_genomes_info.txt", sep="\t", index=False)


def get_genome_to_seqid(genome2id):
    genome2seqid_list = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_genome, genome, id): genome for genome, id in genome2id.items()}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            genome2seqid_list.extend(result)

    # with open(f"{out_dir}/genome2seqid.txt", "w") as f:
    #     f.write("\n".join(genome2seqid_list) + "\n")
    genome2seqid = pd.DataFrame(genome2seqid_list)
    genome2seqid.columns = ["genome_ID", "seq_id"]
    return genome2seqid

def open_file(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")

def process_genome(genome, id):
    sequence_ids = []
    try:
        with open_file(id) as f:
            for line in f:
                if line.startswith(">"):
                    seq_id = line[1:].strip().split(" ", 1)[0]
                    sequence_ids.append(seq_id)
    except FileNotFoundError:
        print(f"File not found: {id}")
        return []
    
    # return [f"{genome}\t{seq_id}" for seq_id in sequence_ids]
    return [(genome, seq_id) for seq_id in sequence_ids]

if __name__ == "__main__":
    sys.exit(main())
