#!/usr/bin/env python3
# for every species only has one strain, merge all strains to a big gfa file and node length chopped by ourselves.
import os, argparse, sys
from typing import Dict, List
from staticsData import read_data
import concurrent.futures

usage = "Build variation graph for the species only one genomes"

class BuildGraph:
    def __init__(self, species_eq1_genomes_info: str, wd: str):
        self.species_eq1_genomes_info = species_eq1_genomes_info
        self.wd = wd
        self.gfa_build_done = []

    def read_fasta(self, genome_path: str) -> Dict[str, str]:
        """
        Read a genome file(.fa, .fna, .fasta).
        """
        genome = read_data(genome_path)
        genome_seqids = list(genome.keys())
        genome_name = os.path.basename(genome_path).split("_")[0] + "_" + os.path.basename(genome_path).split("_")[1]
        new_genome_seqs: Dict[str, str] = {} # [key] = seqid(i.e. NZ...), [value] = seq(fast object)
        for genome_seqid in genome_seqids:
            new_genome_seqid = genome_seqid.split(" ")[0]
            new_genome_seqid = f"{genome_name}#1#{new_genome_seqid}"
            new_genome_seqs[new_genome_seqid] = genome[genome_seqid]
        return new_genome_seqs


    def single_genome_gfa(self, species_taxid: str, genome_id: str, chopped: int = 1024) -> None:
        """
        Cut the genome into blocks of specified length, which serve as nodes. This GFA graph is like a chain.
        Node 1 -> node2 -> ... -> node n
        """
        if not os.path.exists(f"{self.wd}/{species_taxid}.gfa"):        
            # genome_path = "/home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/complete_genome/complete_genome_without_plasmid/GCF_000009485.1_ASM948v1_genomic.fna"
            # genome_path = "/home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/complete_genome/strain-level2/249567/GCF_019880305.1_ASM1988030v1_genomic.fna"
            # genome_path = "249567/GCF_019880305.1_ASM1988030v1_genomic.fna"
            # chopped = 2
            genome_seqs = self.read_fasta(genome_id)
            flag = 0
            S: List[str] = []
            L: List[str] = []
            W: List[str] = []
            for genome_seqid, genome_seq in genome_seqs.items():
                split_genome_seq = [genome_seq[i:i+chopped] for i in range(0, len(genome_seq), chopped)]
                s_validate_genome_len = 0
                for i in range(len(split_genome_seq)):
                    S.append(f"S\t{i+flag+1}\t{split_genome_seq[i]}\n")
                    s_validate_genome_len = s_validate_genome_len + len(split_genome_seq[i])
                if s_validate_genome_len != len(genome_seq):
                    print("S lines errors")
                    break
                for i in range(len(split_genome_seq)-1):
                    i = i + flag
                    L.append(f"L\t{i+1}\t+\t{i+2}\t+\t0M\n")
                
                link = ">" + ">".join(str(i+flag+1) for i in range(len(split_genome_seq)))
                tokens = genome_seqid.split("#")
                W.append(f"W\t{tokens[0]}\t{tokens[1]}\t{tokens[2]}\t0\t{len(genome_seq)}\t{link}\n")
                flag = flag + len(split_genome_seq)
            if len(S) != len(L) + len(genome_seq) and len(W) != len(genome_seqs):
                print("L lines or W lines error")
            
            with open(f"{self.wd}/{species_taxid}.gfa", "w") as f:
                f.write("H\tVN:Z:1.1\n")
                f.write("".join(S))
                f.write("".join(L))
                f.write("".join(W))
        self.gfa_build_done.append(f"{self.wd}/{species_taxid}.gfa")

    def parallel_build_gfa(self) -> None:
        with open(self.species_eq1_genomes_info, "r") as f:
            species_to_genome: Dict[str, str] = {} # [key] = species_taxid, [value] = genomeID(i.e. GCF...)
            for line in f:
                tokens = line.strip().split("\t")
                species_taxid = tokens[0]
                genome_id = tokens[1]
                species_to_genome[species_taxid] = genome_id

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {key: executor.submit(self.single_genome_gfa, key, value) for key, value in species_to_genome.items()}
            concurrent.futures.wait(futures.values())

    def write_gfa_building_done(self):
        with open(f"{self.wd}/gfa_build_done.txt", "w") as f:
            f.write("\n".join(self.gfa_build_done) + "\n")

def main():
    parser = argparse.ArgumentParser(prog="python build_graph.py", description=usage)
    parser.add_argument("species_eq1_genomes_info", type=str, help="Species eq1 genomes information file. Format: species_taxid\\tID\\n")
    parser.add_argument("wd", type=str, help="Work directory")
    args = parser.parse_args()
    build_graph=BuildGraph(args.species_eq1_genomes_info, args.wd)
    build_graph.parallel_build_gfa()
    build_graph.write_gfa_building_done()

if __name__=="__main__":
    sys.exit(main())