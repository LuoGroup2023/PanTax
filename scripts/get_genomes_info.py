#!/usr/bin/env python3
import pandas as pd
import os, sys, argparse, re

usage = "Get the filtered genomes information in every species cluster. These genomes are used to construct pangenome"

def main():
    parser = argparse.ArgumentParser(prog="python get_genomes_info.py", description=usage) 
    parser.add_argument("-f", dest="genomes_info_provided_origin", default="genomes_info_provided_origin.txt", type=str, help="File(genomes_info_provided_origin.txt) produced in the cluster step")
    parser.add_argument("-c", "--custom_genomes", dest="genomes_info_provided_custom", type=str, help="If use custom genomes information, this step will filter unused genomes information")
    parser.add_argument("-o", "--output_dir", dest="output_dir", type=str, help="Output dir")
    parser.add_argument("--gtdb", dest="gtdb", action="store_true", help="GTDB genome accession.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    if args.genomes_info_provided_custom:
        try:
            used_genomes = pd.read_csv("used_genomes.txt", header=None).iloc[:,0].tolist()
        except:
            used_genomes = []
            raise IOError("No genomes was selected to cluster.")
        if not args.gtdb:
            used_genomes = [re.sub(r"_genomic\.fna(\.gz)?$", "", os.path.basename(genome)) for genome in used_genomes]
        else:
            used_genomes = ["_".join(os.path.basename(genome).split("_")[:2]) for genome in used_genomes]
        used_genomes_df = pd.DataFrame({'genome_ID': used_genomes})
        genomes_info_provided_custom = pd.read_csv(args.genomes_info_provided_custom, sep="\t")
        used_genomes_info = pd.merge(used_genomes_df, genomes_info_provided_custom, on="genome_ID", how="left")
        if args.gtdb:
            used_genomes_info["species_taxid"] = used_genomes_info["species_taxid"].str.replace(" ", "_")
            # print(used_genomes_info["species_taxid"].tolist()[:2])
        process_map(used_genomes_info, args.output_dir)
    elif args.genomes_info_provided_origin:
        genomes_info = pd.read_csv(args.genomes_info_provided_origin, sep="\t")
        if args.gtdb:
            genomes_info["species_taxid"] = genomes_info["species_taxid"].str.replace(" ", "_")
        process_map(genomes_info, args.output_dir)
    else:
        print("No any provided file.\n")

def process_map(genomes_info, output_dir):    
    group_counts = genomes_info["strain_taxid"].value_counts()
    genomes_info["group_index"] = genomes_info.groupby("strain_taxid").cumcount() + 1
    # mask = (genomes_info["strain_taxid"] == genomes_info["species_taxid"]) & (genomes_info["strain_taxid"].map(group_counts) > 1)
    mask = genomes_info["strain_taxid"].map(group_counts) > 1
    genomes_info["group_index"] = genomes_info["group_index"].astype(str).where(mask, "")
    genomes_info["new_strain_taxid"] = genomes_info["strain_taxid"].astype(str) + "." + genomes_info["group_index"].astype(str)
    genomes_info["new_strain_taxid"] = genomes_info["strain_taxid"].astype(str).where(genomes_info["group_index"] == "", genomes_info["new_strain_taxid"])
    genomes_info["strain_taxid"] = genomes_info["new_strain_taxid"]
    genomes_info = genomes_info.drop(["new_strain_taxid", "group_index"], axis=1)
    genomes_info.to_csv(f"{output_dir}/genomes_info.txt", index=False, sep="\t")


if __name__ == "__main__":
    sys.exit(main())