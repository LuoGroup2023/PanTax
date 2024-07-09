#!/usr/bin/env python3
import sys, os, subprocess, argparse, re
import pandas as pd
from typing import Dict, List
import concurrent.futures

usage = "Find OTU node range in reference pangenome graph"

class OtuRangeProcessor:

    def __init__(self, input_file: str, genomes_info: str, extract_hap_paths_file: str):
        self.input_file = input_file
        self.genomes_info = genomes_info
        self.extract_hap_paths_file = extract_hap_paths_file
        self.flag = self.hap_line_start_in_gfa()
        self.taxid_to_index = None

    def hap_line_start_in_gfa(self):
        """
        Obtain GAF version.
        """
        with open(self.input_file, "r") as f:
            line = next(f)
            version = line.strip().split(":")[-1].strip()
            if version == "1.0": # P line
                return 1
            elif version == "1.1": # W line
                return 0
            else:
                print(f"VersionError: GFA specification version {version} not supported")
                sys.exit(1)
                
    def get_species_to_genomes_info(self) -> Dict[str, List[str]]:
        """
        Obtain mapping relation of species taxid and genomes.
        """
        species_to_genomes: Dict[str, List[str]] = {}
        genomes_info= pd.read_csv(self.genomes_info, sep="\t")
        genomes_info["species_taxid"] = genomes_info["species_taxid"].astype(str)
        grouped_genomes_info = genomes_info.groupby("species_taxid")
        for species_taxid, group in grouped_genomes_info:
            species_to_genomes[species_taxid] = group["genome_ID"].tolist()
        return species_to_genomes

    def get_strain_to_genomes_info(self):
        genome_ID = []
        strain_taxid = []

        with open(self.genomes_info, "r") as f:
            next(f)
            for line in f:
                tokens = line.strip().split("\t")
                genome_ID.append([tokens[0]])
                strain_taxid.append(tokens[1])
        strain_to_genomes = dict(zip(strain_taxid, genome_ID))
        return strain_to_genomes

    def extract_hap_paths(self):
        """
        Extract haplotype path from GAF file.
        """
        if not os.path.exists(self.extract_hap_paths_file):
            if self.flag:
                awk_command = f"awk '$1 == \"P\" {{print}}' {self.input_file}"
            else:
                awk_command = f"awk '$1 == \"W\" {{print}}' {self.input_file}"
            subprocess.run(awk_command, shell=True, check=True, stdout=open(self.extract_hap_paths_file, "w"))

    def OTU_to_hap(self, otu_to_genomes: Dict[str, List[str]]):
        """
        get mapping relation of taxid and index in extract_hap_paths_file.
        GCF1 index:1 species_taxid:24
        GCF2 index:2 species_taxid:34
        GCF3 index:3 species_taxid:24
        taxid_to_index = [24:[1,3], 34:[2]]
        """
        df1 = pd.read_csv(self.extract_hap_paths_file, sep="\t", header=None, usecols=[1], names=["genomes"])
        if self.flag:
            df1["genomes"] = df1["genomes"].apply(lambda x: x.split("#")[0])
        df1 = df1.reset_index()
        data_list = [(taxid, "GCF_" + genome.split("_")[1]) for taxid, genomes in otu_to_genomes.items() for genome in genomes]
        df2 = pd.DataFrame(data_list, columns=["taxid", "genomes"])
        merge_df = pd.merge(df1, df2, on="genomes", how="left")
        taxid_set = merge_df["taxid"].unique().tolist()
        self.taxid_to_index: Dict[str, List[int]] = {}
        for taxid in taxid_set:
            index = merge_df.loc[merge_df["taxid"] == taxid, "index"].tolist()
            self.taxid_to_index[taxid] = index
            # if len(taxid_to_index) >= 10:
            #     break
        return self.taxid_to_index

    def get_otu_nodeID_range(self, taxid: str, index: List[int]):
        """
        get mapping relation of taxid and haplotype path, and then process. Some species have more than one path.
        """
        paths: List[int] = []
        if self.flag:
            extract_row = 2
        else:
            extract_row = 6
        with open(self.extract_hap_paths_file, "r") as file:
            for i, line in enumerate(file):
                if i in index:
                    row_data = line.strip().split('\t')
                    selected_data = row_data[extract_row]
                    paths.append(selected_data)
        min_values: List[int] = []
        max_values: List[int] = []
        for path in paths:
            if self.flag:
                path = [int(match.group()) for match in re.finditer(r'\d+', path)]                             
            else:
                path = [int(match.group()) for match in re.finditer(r'-?\d+', path)]
            min_values.append(min(path))
            max_values.append(max(path))
        min_value = min(min_values)
        max_value = max(max_values)     
        return taxid, [min_value, max_value]

    def sort_dict(self, result_dict: Dict[str, List[int]], flag: str) -> Dict[str, List[int]]:
        """
        sort and validate
        """
        # result_dict = {9: [1, 10696], 34: [190499, 1718342], 24: [10697, 190498]}
        sorted_data = dict(sorted(result_dict.items(), key=lambda item: item[1][0]))
        if flag == "species":
            sorted_keys = list(sorted_data.keys())
            count = 1
            for i in range(1, len(sorted_keys)):
                current_key = sorted_keys[i]
                prev_key = sorted_keys[i - 1]
                if sorted_data[current_key][0] == sorted_data[prev_key][1] + 1:
                    count += 1
                else:
                    print("The order is error")
                    sys.exit(1)
            assert count == len(sorted_keys)
            # if count == len(sorted_keys):
                # print("the order is right")
        return sorted_data

    def parallel(self, flag: str) -> Dict[str, List[int]]:
        result_dict = {}
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {executor.submit(self.get_otu_nodeID_range, key, value): key for key, value in self.taxid_to_index.items()}
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    taxid, range_values = result
                    result_dict[taxid] = range_values
                except Exception as e:
                    print(f"An error occurred: {e}")
        result_dict = self.sort_dict(result_dict, flag)
        # print(result_dict)
        return result_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="python otu_range.py", description=usage)
    parser.add_argument("reference_pangenome_gfa", type=str, help="Reference pangenome GFA file")
    parser.add_argument("genomes_info", type=str, help="Genomes information file")
    parser.add_argument("--species-level", dest="species_flag", action="store_true", help="OTU(species) level range. on(1)/off(0)")
    parser.add_argument("--strain-level", dest="strain_flag", action="store_true", help="OTU(species) level range. on(1)/off(0)")
    args = parser.parse_args()
    print("\nProgram settings:\n")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))
    print()
    extract_hap_paths_file = "extract_hap_paths.gfa"    
    otu_range_processor = OtuRangeProcessor(args.reference_pangenome_gfa, args.genomes_info, extract_hap_paths_file)
    otu_range_processor.extract_hap_paths()
    if args.species_flag:
        species_to_genomes = otu_range_processor.get_species_to_genomes_info()
        otu_range_processor.OTU_to_hap(species_to_genomes)
        result_dict = otu_range_processor.parallel(flag="species")
        with open("species_range.txt", "w") as f:
            for key, value in result_dict.items():
                f.write(f"{key}\t{value[0]}\t{value[1]}\n")
    if args.strain_flag:
        strain_to_genomes = otu_range_processor.get_strain_to_genomes_info()
        otu_range_processor.OTU_to_hap(strain_to_genomes)
        result_dict = otu_range_processor.parallel(flag="strain")
        with open("strain_range.txt", "w") as f:
            for key, value in result_dict.items():
                f.write(f"{key}\t{value[0]}\t{value[1]}\n")    