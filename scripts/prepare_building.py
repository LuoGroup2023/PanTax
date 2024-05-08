import sys
import os
import argparse
import pandas as pd

usage = "Prepare genomes information for every species to build pangenome with pggb"

def main():
    parser = argparse.ArgumentParser(prog="python prepare_building.py", description=usage)
    parser.add_argument("genomes_info", type=str, help="Genomes information file.Format: genomeID\tstrain_taxid\tspecies_taxid\ttorganism_name\tid.")
    args = parser.parse_args()
    if not os.path.exists("pggb_species"):
        prepare_species_genomes(args.genomes_info)

def prepare_species_genomes(genomes_info):
    genomes_info_df = pd.read_csv(genomes_info, sep="\t")
    species_ge2 = pd.DataFrame(genomes_info_df[genomes_info_df.duplicated("species_taxid", keep=False)])
    grouped = species_ge2.groupby("species_taxid")
    if len(grouped) != 0:
        os.mkdir("pggb_species")
        for species_taxid, group_df in grouped:
            genomes = group_df["id"].tolist()
            with open(f"pggb_species/{species_taxid}.txt", "w") as f:
                f.write("\n".join(genomes) + "\n")
    species_eq1 = pd.DataFrame(genomes_info_df[~genomes_info_df.duplicated("species_taxid", keep=False)])
    species_taxid = species_eq1["species_taxid"].tolist()
    genomes = species_eq1["id"].tolist()
    with open("species_eq1_genomes.txt", "w") as f:
        for i in range(len(species_taxid)):
            f.write(f"{species_taxid[i]}\t{genomes[i]}\n")

if __name__ == "__main__":
    sys.exit(main())
            

    
