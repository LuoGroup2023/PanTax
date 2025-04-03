#!/usr/bin/env python3
import sys, os, argparse, subprocess, shutil, gzip
from toolkits import Logger
import concurrent.futures

usage = "Extract complete genomes from Refseq database and remove plasmids"

def main():
    parser = argparse.ArgumentParser(prog="python extract_complete_genome.py", description=usage)
    parser.add_argument("-r", "--refseq_database", dest="refseq_database", type=str, help="Refseq database path")
    parser.add_argument("-o", "--output_database", dest="output_database", default="complete_genome_without_plasmid", type=str, help="Output complete genomes database path")
    parser.add_argument("-s", "--summary_file", dest="summary_file_path", default="assembly_summary_bacteria.txt", type=str, help="Assembly summary file path")
    parser.add_argument("-g", "--output_genomes_info", dest="output_genomes_info", type=str, help="Output genomes information file")
    # parser.add_argument("-l", "--genome_assembly_lvl", dest="genome_assembly_lvl", default="complete", type=str, help="Genome assembly level.(all/complete).")
    # parser.add_argument("-p", "--cluster", dest="species_cluster", default="all", type=str, help="Species cluster.")
    parser.add_argument("-c", "--custom", dest="custom_complete_genomes", default="custom_complete_genomes.txt", type=str, help="Specify custom complete genomes database(use without -r -s). Format: genomeID\tstrain_taxid\tspecies_taxid\ttorganism_name\tid.")
    parser.add_argument("--remove", dest="remove", default=True, type=bool, help="Wether to remove plasmid")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()
    log = Logger()
    global remove
    remove = args.remove
    if not os.path.exists(args.output_database):
        os.mkdir(args.output_database)
    args.output_database = os.path.abspath(args.output_database)
    if not os.path.exists(args.summary_file_path) and not os.path.exists(args.custom_complete_genomes):
        download_assembly_summary(args.output_database)
        args.summary_file_path = os.path.join(args.output_database, "assembly_summary_bacteria.txt")
    if not args.output_genomes_info:
        args.output_genomes_info = os.path.join(args.output_database, "genomes_info_provided_origin.txt")
    if not os.path.exists(args.custom_complete_genomes):
        complete_genomes_path = read_complete_genome(args.summary_file_path, args.refseq_database)
        get_genomes_info(args.summary_file_path, args.output_database, complete_genomes_path, args.output_genomes_info)
    else:
        complete_genomes_path = []
        map_list = []
        with open(args.custom_complete_genomes, "r") as f:
            next(f)
            for line in f:
                origin_info = line.strip().split("\t")
                complete_genomes_path.append(origin_info[4])
                origin_info[4] = os.path.join(args.output_database, os.path.basename(origin_info[4]))
                map_list.append("\t".join(origin_info))
        with open(args.output_genomes_info, "w") as f:
            f.write("genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid\n")
            f.write("\n".join(map_list) + "\n")
    complete_genomes = "complete_genomes.txt"
    with open(complete_genomes, "w") as f:
        f.write("\n".join(complete_genomes_path) + "\n")

    parallel_extract(complete_genomes, args.output_database)
    log.logger.info(f"Extract complete genomes successfully, the complete genome database in {os.path.abspath(args.output_database)}")
    return

def download_assembly_summary(db_dir, division="bacteria", q = 0):
    path = os.path.join(db_dir, "assembly_summary_" + division + ".txt")
    if os.path.exists(path):
        print(f"assembly_summary_{division}.txt exists")
    else:
        assembly_summary_dir = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/" + division + "/assembly_summary.txt"
        if q == 0:
            process = subprocess.Popen(["wget", assembly_summary_dir], stdout=subprocess.PIPE, universal_newlines=True)
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print(output.strip())
        else:
            process = subprocess.Popen(["wget", "-q", assembly_summary_dir])
        process.wait()
        process = subprocess.Popen(["mv", "assembly_summary.txt", path])
        process.wait()
        if os.path.exists(path):
            print(f"assembly_summary_{division}.txt download successfully")
        else:
            print(f"assembly_summary_{division}.txt download fail")

def read_complete_genome(summary_file_path, Refseq_database_path):
    """get complete genomes path from local Refseq database."""
    complete_genomes_path = []
    with open(summary_file_path, "r") as f:
        f.readline()
        f.readline()
        for line in f:
            line = line.strip()
            tokens = line.split("\t")
            if tokens[11] == "Complete Genome":
                complete_genome = os.path.basename(tokens[19]) + "_genomic.fna"
                absolute_genome_path = os.path.join(Refseq_database_path, complete_genome)
                if os.path.exists(absolute_genome_path):
                    complete_genomes_path.append(absolute_genome_path)
    if not complete_genomes_path:
        print("Genome selected error. Your genomes are all not in summary file. Please set custom option.")
        sys.exit(1)
    complete_genomes_path = list(set(complete_genomes_path))
    return complete_genomes_path

def open_file(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")

def extract_complete_sequence(genome_info):
    """remove plasmids"""
    genome_file = genome_info[0]
    new_genome_file_path = genome_info[1]
    new_genome_file_path = os.path.join(new_genome_file_path, os.path.basename(genome_file))
    if not os.path.exists(new_genome_file_path):
        if remove:
            all_sequence = {}
            sequence_name = None
            sequence = []
            with open_file(genome_file) as file:
                for line in file:
                    line = line.strip()
                    if line.startswith(">"):
                        if sequence_name:
                            all_sequence[sequence_name] = "\n".join(sequence)     
                        sequence_name = line
                        sequence = [] 
                    else:
                        sequence.append(line)
                if sequence_name:
                    all_sequence[sequence_name] = "\n".join(sequence)
            genome_sequence_names = [key for key in all_sequence.keys() if "plasmid" not in key]
            chromosome_genome_data = []
            if len(genome_sequence_names) >= 1:
                for genome_sequence_name in genome_sequence_names:
                    chromosome_genome_data.append(genome_sequence_name)
                    genome_sequence = all_sequence[genome_sequence_name]
                    chromosome_genome_data.append(genome_sequence)
                with open(new_genome_file_path, "w") as f:
                    f.write("\n".join(chromosome_genome_data) + "\n")
            else:
                print(genome_file)
        else:
            shutil.copyfile(genome_file, new_genome_file_path)

def parallel_extract(complete_genomes, complete_genomes_database_path):
    """parallel remove plasmids"""
    with open(complete_genomes, "r") as f:
        genome_file_set = [line.strip() for line in f]
    genome_info = list(zip(genome_file_set, [complete_genomes_database_path]*len(genome_file_set)))
    # print(genome_info[0:2])
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(extract_complete_sequence, genome_info)

def get_genomes_info(summary_file_path, output_database, complete_genomes_path, output_genomes_info):
    genomes = [os.path.basename(genome).replace("_genomic.fna", "") for genome in complete_genomes_path]
    map_list = []
    with open(summary_file_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split("\t")
            genome = os.path.basename(tokens[19])
            if genome in genomes:
                map_list.append(f"{genome}\t{tokens[5]}\t{tokens[6]}\t{tokens[7]}\t{output_database}/{genome}_genomic.fna")
                genomes.remove(genome)
    with open(output_genomes_info, "w") as f:
        f.write("genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid\n")
        for map in map_list:
            f.write(map + "\n")

if __name__ == "__main__":
    sys.exit(main())