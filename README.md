# PanTax

## Installation
```
git clone https://github.com/LuoGroup2023/PanTax.git
cd PanTax
conda env create -f environment.yaml
sh install.sh
```

## Genome preprocessing

We recommend removing plasmids and redundancy from the genome first with --remove option and --cluster option, respectively. Eventually you will get a file containing genomic information in /path/to/database/library.
If genomes are all in NCBI refseq database, you only need to use -r option to specify the directory containing these genomes.
```
sh data_preprocessing.sh -r ref --remove --cluster
```
Otherwise, you need to provide a file containing information about the custom genomes.
```
sh data_preprocessing.sh -c genomes_info.txt --remove --cluster
```
The `genomes_info.txt` file gives a list of reference genomes in fasta format, which constitute PaxTax's original database, alongwith NCBI's taxonomic information. The input lines in the file should contain at least 5 tab-delimited fields; from left to right, they are Genome IDs, Strain taxonomic IDs, Species taxonomic IDs, Organism names, Genome absolute path.
Here is an example format of `genomes_info.txt` file:
```
genome_ID	strain_taxid	species_taxid	organism_name	id
GCF_000218545.1_ASM21854v1	593907	11	Cellulomonas gilvus ATCC 13127	/path/to/GCF_000218545.1_ASM21854v1_genomic.fna
GCF_025402875.1_ASM2540287v1	24.1	24	Shewanella putrefaciens	/path/to/GCF_025402875.1_ASM2540287v1_genomic.fna
```

## Running
* create database only 
```
sh pantax -f $genome_info --create
```
* Query with specified database
```
# long read
sh pantax -f $genome_info -l -r $fq -db $db --species-level --strain-level
# short read
sh pantax -f $genome_info -s -p -r $fq -db $db --species-level --strain-level
```

## Examples

* short read
```
cd PanTax/example/sample_hifi_data
sh ../../scripts/pantax -f ../sample_genome_info.txt -l -r long_reads.fq.gz --species-level --strain-level
```
* long read
```
cd PanTax/example/sample_ngs_data
sh ../../scripts/pantax -f ../sample_genome_info.txt -s -p -r short_reads.fq.gz --species-level --strain-level
```