#!/usr/bin/env python3
import sys, os, argparse, shutil, gzip, subprocess
from toolkits import Logger
from pathlib import Path
import concurrent.futures
import multiprocessing
# file_lock = multiprocessing.Lock()
usage = "Remove plasmids and obtain genomes_info file from specified assembly_summary file and existing genomes"

def main():
    parser = argparse.ArgumentParser(prog="python genomes_process.py", description=usage)
    parser.add_argument("-i", "--input_fasta", dest="input_fasta", type=str, help="The dir of input fasta genome.")
    parser.add_argument("-o", "--output_dir", dest="output_dir", default="reference_genomes", type=str, help="Output directory.")
    parser.add_argument("-s", "--summary_file", dest="summary_file_path", default="assembly_summary_bacteria.txt", type=str, help="Assembly summary file path.")
    parser.add_argument("-l", "--genome_assembly_lvl", dest="genome_assembly_lvl", default="all", type=str, help="Genome assembly level.(all/complete).")
    parser.add_argument("-p", "--cluster", dest="species_cluster", default="all", type=str, help="Species cluster.")
    parser.add_argument("-g", "--output_genomes_info", dest="output_genomes_info", type=str, help="Output genomes information file.")
    parser.add_argument("-d", "--db", dest="db", type=str, help="Database source(rs, gb, gtdb).")
    parser.add_argument("-c", "--custom", dest="custom_genomes", default=None, type=str, help="Specify custom genomes database(use without -r -s). Format: genomeID\tstrain_taxid\tspecies_taxid\ttorganism_name\tid.")
    parser.add_argument("-f", "--gtdb", dest="gtdb_metadata", type=str, help="GTDB metadata file. Specify this file will use GTDB taxonomy.")
    parser.add_argument("-t", "--threads", dest="threads", default=64, type=int, help="Number of threads.")
    parser.add_argument("--remove", dest="remove", action="store_true", help="Whether to remove plasmid.")
    parser.add_argument("-rl", "--remove_scaffold_len", dest="remove_scaffold_len", default=1, type=int, help="Whether to remove remove <length>Mbp scaffold.")
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
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    if args.custom_genomes == "None": args.custom_genomes = None
    if args.species_cluster == "None": args.species_cluster = None
    if args.db.lower() == "rs" and not os.path.exists(args.summary_file_path) and not args.custom_genomes:
        download_assembly_summary(args.output_dir)
        args.summary_file_path = os.path.join(args.output_dir, "assembly_summary_bacteria.txt")
    elif args.db == "gb" and not os.path.exists(args.summary_file_path):
        log.logger.error("Please provide GB assembly_summary file.")

    if not args.output_genomes_info:
        args.output_genomes_info = os.path.join(args.output_dir, "genomes_info_provided_origin.txt")
    if not args.custom_genomes:
        genomes_path = obtain_existing_genomes_path(args.summary_file_path, args.gtdb_metadata, args.input_fasta, args.genome_assembly_lvl, args.species_cluster)
        # genomes_path = [str(f) for f in Path(args.input_fasta).iterdir() if f.is_file()]
        genomes_path = get_genomes_info(args.summary_file_path, args.gtdb_metadata, args.output_dir, genomes_path, args.output_genomes_info, args.genome_assembly_lvl)
    elif args.custom_genomes and os.path.exists(args.custom_genomes):
        genomes_path = []
        map_list = []
        with open(args.custom_genomes, "r") as f:
            next(f)
            for line in f:
                origin_info = line.strip().split("\t")
                genomes_path.append(origin_info[4])
                origin_info[4] = os.path.join(args.output_dir, os.path.basename(origin_info[4]))
                map_list.append("\t".join(origin_info))
        with open(args.output_genomes_info, "w") as f:
            f.write("genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid\n")
            f.write("\n".join(map_list) + "\n")
    else:
        log.logger.error("No enough information.")
    genomes = "genomes.txt"
    with open(genomes, "w") as f:
        f.write("\n".join(genomes_path) + "\n")

    parallel_extract(genomes, args.output_dir, args.threads, args.remove_scaffold_len)
    log.logger.info(f"Extract genomes successfully, the genome database in {os.path.abspath(args.output_dir)}")
    return

def download_assembly_summary(db_dir, division="bacteria", q = 0):
    path = os.path.join(db_dir, "assembly_summary_" + division + ".txt")
    if os.path.exists(path):
        print(f"assembly_summary_{division}.txt exists")
    else:
        assembly_summary_dir = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/" + division + "/assembly_summary.txt"
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


def obtain_existing_genomes_path(summary_file_path, gtdb_metadata, genomes_dir, genome_assembly_lvl, species_clusters):
    """get genomes path from local database directory."""
    assert isinstance(genome_assembly_lvl, str)
    if "complete" in genome_assembly_lvl.lower():
        genome_assembly_lvl = "Complete Genome"
    
    if species_clusters:
        species_clusters = species_clusters.strip().split(",")

    genomes_path = []
    if gtdb_metadata and os.path.exists(gtdb_metadata):
        if species_clusters:
            species_clusters = [species_cluster.replace(" ", "_") for species_cluster in species_clusters]
        with open(gtdb_metadata, "r") as f:
            next(f)
            for line in f:
                tokens = line.strip().split("\t")
                if genome_assembly_lvl == "all" or tokens[48] == genome_assembly_lvl:      # tokens[48]: ncbi_assembly_level
                    if species_clusters:
                        gtdb_species = tokens[19].strip().split(";")[-1].split("__")[1].replace(" ", "_")
                        if gtdb_species in species_clusters:
                            accession = tokens[0].strip().split("_", 1)[1]
                            ncbi_assembly_name = tokens[49]
                            genome = accession + "_" + ncbi_assembly_name + "_genomic.fna"
                            absolute_genome_path = os.path.join(genomes_dir, genome)
                            gzip_absolute_genome_path = absolute_genome_path + ".gz"
                            if os.path.exists(absolute_genome_path):
                                genomes_path.append(absolute_genome_path)
                            elif os.path.exists(gzip_absolute_genome_path):
                                genomes_path.append(gzip_absolute_genome_path)
                            # else:
                            #     print(f"{absolute_genome_path} or {gzip_absolute_genome_path} does not exist.")
                    else:
                        accession = tokens[0].strip().split("_", 1)[1]
                        ncbi_assembly_name = tokens[49]
                        genome = accession + "_" + ncbi_assembly_name + "_genomic.fna"
                        absolute_genome_path = os.path.join(genomes_dir, genome)
                        gzip_absolute_genome_path = absolute_genome_path + ".gz"
                        if os.path.exists(absolute_genome_path):
                            genomes_path.append(absolute_genome_path)
                        elif os.path.exists(gzip_absolute_genome_path):
                            genomes_path.append(gzip_absolute_genome_path)
                        # else:
                        #     print(f"{absolute_genome_path} or {gzip_absolute_genome_path} does not exist.")
    elif summary_file_path and os.path.exists(summary_file_path):
        with open(summary_file_path, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                line = line.strip()
                tokens = line.split("\t")
                if genome_assembly_lvl == "all" or tokens[11] == genome_assembly_lvl:
                    if species_clusters:
                        if tokens[6] in species_clusters:
                            genome = os.path.basename(tokens[19]) + "_genomic.fna"
                            absolute_genome_path = os.path.join(genomes_dir, genome)
                            gzip_absolute_genome_path = absolute_genome_path + ".gz"
                            if os.path.exists(absolute_genome_path):
                                genomes_path.append(absolute_genome_path)
                            elif os.path.exists(gzip_absolute_genome_path):
                                genomes_path.append(gzip_absolute_genome_path)
                            else:
                                print(f"{absolute_genome_path} or {gzip_absolute_genome_path} does not exist.")
                    else:
                        genome = os.path.basename(tokens[19]) + "_genomic.fna"
                        absolute_genome_path = os.path.join(genomes_dir, genome)
                        gzip_absolute_genome_path = absolute_genome_path + ".gz"
                        if os.path.exists(absolute_genome_path):
                            genomes_path.append(absolute_genome_path)
                        elif os.path.exists(gzip_absolute_genome_path):
                            genomes_path.append(gzip_absolute_genome_path)
                        else:
                            print(f"{absolute_genome_path} or {gzip_absolute_genome_path} does not exist.")
    else:
        print(f"{summary_file_path} or {gtdb_metadata} does not exist. You must provide assembly_summary file or GTDB_metadata file.")        
        
    if not genomes_path:
        print("Genome selected error. Your genomes are all not in summary file. Please set custom option.")
        sys.exit(1)
    genomes_path = list(set(genomes_path))
    return genomes_path

def open_file(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")

def extract_complete_sequence(genome_info, remove_scaffold_len, progress_counter):
    """remove plasmids"""
    remove_scaffold_len = remove_scaffold_len*1000000
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
            genome_sequence_names = [key for key, value in all_sequence.items() if len(value) >= remove_scaffold_len]
            chromosome_genome_data = []
            if len(genome_sequence_names) >= 1:
                for genome_sequence_name in genome_sequence_names:
                    chromosome_genome_data.append(genome_sequence_name)
                    genome_sequence = all_sequence[genome_sequence_name]
                    chromosome_genome_data.append(genome_sequence)
                if genome_file.endswith("gz"):
                    with gzip.open(new_genome_file_path, "wb") as f:
                        f.write(("\n".join(chromosome_genome_data) + "\n").encode("utf-8"))
                else:
                    with open(new_genome_file_path, "w") as f:
                        f.write("\n".join(chromosome_genome_data) + "\n")
            else:
                print(f"{genome_file} all scaffolds are smaller than 1Mbp.")
            # Locking to safely write to the output file
            # with file_lock:
            #     with open(progress_counter["output_file"], "a") as output_f:
            #         output_f.write(new_genome_file_path + "\n")

        else:
            shutil.copyfile(genome_file, new_genome_file_path)

        progress_counter["completed"] += 1
        print(f"Completed {new_genome_file_path}: {progress_counter['completed']}/{progress_counter['total']}")

def parallel_extract(genomes, genomes_database_path, threads, remove_scaffold_len):
    """parallel remove plasmids"""
    with open(genomes, "r") as f:
        genome_file_set = [line.strip() for line in f]
    genome_info = list(zip(genome_file_set, [genomes_database_path]*len(genome_file_set)))
    # Shared file lock and progress counter
    
    progress_counter = multiprocessing.Manager().dict()
    progress_counter["completed"] = 0
    progress_counter["total"] = len(genome_info)
    progress_counter["output_file"] = "finished_add_genomes.txt"
    
    if threads > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(extract_complete_sequence, genome, remove_scaffold_len, progress_counter)
                for genome in genome_info
            ]
            
            for future in futures:
                future.result()
    elif threads == 1:
        for _genome_info in genome_info:
            extract_complete_sequence(_genome_info, remove_scaffold_len, progress_counter)

def get_genomes_info(summary_file_path, gtdb_metadata, output_dir, genomes_path, output_genomes_info, genome_assembly_lvl):
    assert isinstance(genome_assembly_lvl, str)
    if "complete" in genome_assembly_lvl.lower():
        genome_assembly_lvl = "Complete Genome"
    if genomes_path[0].endswith(".gz"): 
        file_gzip = True
    else:
        file_gzip = False
    if file_gzip:
        genomes = [os.path.basename(genome).replace("_genomic.fna.gz", "") for genome in genomes_path]
    else:
        genomes = [os.path.basename(genome).replace("_genomic.fna", "") for genome in genomes_path]
    
    map_list = []
    existing_genomes_path = []
    if gtdb_metadata and os.path.exists(gtdb_metadata):
        genomes_gtdb_taxonomy = {}
        with open(gtdb_metadata, "r") as f:
            next(f)
            for line in f:
                tokens = line.strip().split("\t")
                if genome_assembly_lvl == "all" or tokens[48] == genome_assembly_lvl:
                    accession = tokens[0].strip().split("_", 1)[1]
                    gtdb_taxonomy = tokens[19]
                    gtdb_species = gtdb_taxonomy.strip().split(";")[-1]
                    assert gtdb_species.startswith("s")
                    gtdb_species = gtdb_species.split("__")[1]
                    genomes_gtdb_taxonomy[accession] = [gtdb_species, tokens[65]]
        genomes = list(map(lambda x: "_".join(x.split("_")[:2]), genomes))
        count = 5000000
        for i, genome in enumerate(genomes):
            if genome in genomes_gtdb_taxonomy:
                count += 1
                map_list.append(f"{genome}\t{str(count)}\t{genomes_gtdb_taxonomy[genome][0]}\t{genomes_gtdb_taxonomy[genome][1]}\t{output_dir}/{Path(genomes_path[i]).name}")                
                existing_genomes_path.append(genomes_path[i])             
            else:
                print(f"{genome} does not exist GTDB metadata file.")
    elif summary_file_path and os.path.exists(summary_file_path):   
        taxonomy = {}
        with open(summary_file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                tokens = line.strip().split("\t")
                genome = os.path.basename(tokens[19])
                if genome_assembly_lvl == "all" or tokens[11] == genome_assembly_lvl:
                    taxonomy[genome] = [tokens[5], tokens[6], tokens[7]]
        for i, genome in enumerate(genomes):
            if genome in taxonomy:
                map_list.append(f"{genome}\t{taxonomy[genome][0]}\t{taxonomy[genome][1]}\t{taxonomy[genome][2]}\t{output_dir}/{Path(genomes_path[i]).name}")
                existing_genomes_path.append(genomes_path[i])
    else:
        print(f"{summary_file_path} or {gtdb_metadata} does not exist. You must provide assembly_summary file or GTDB_metadata file.")                
    with open(output_genomes_info, "w") as f:
        f.write("genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid\n")
        for mapping in map_list:
            f.write(mapping + "\n")
    return existing_genomes_path

if __name__ == "__main__":
    sys.exit(main())
