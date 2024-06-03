# PanTax

PanTax is a pangenome graph-based taxonomic classification tool designed to overcome the limitations of traditional single-reference genome approaches. It excels in providing accurate taxonomic classification at the strain level, handling both short and long reads, and supporting single or multiple species. By leveraging pangenome graphs, PanTax captures the genetic diversity across multiple related genomes, thereby eliminating the biases introduced by using a single linear representative genome. 

## Installation
The dependencies of PaxTax can generally be installed through conda. There are two dependencies that require special attention. 

Firstly, [pggb](https://github.com/pangenome/pggb.git) depends on vg version 1.40, which is incompatible with the version of vg we use for read alignment. We recommend using a more recent version of [vg](https://github.com/vgteam/vg.git), specifically vg>=1.52. We have included the vg 1.52 executable file (default path) in the tools directory. If this is not available, we suggest creating a new environment with conda and installing vg. You can then place it in the tools directory via a symbolic link or specify its path directly using the `--vg` option.

Secondly, PanTax relies on the [Gurobi Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer/) Python package. While the package download is included in the environment.yaml file, a license must be obtained from the official website to use it for large model optimizations.


```
git clone https://github.com/LuoGroup2023/PanTax.git
cd PanTax
conda env create -f environment.yaml
sh install.sh
conda activate Pantax

# If vg is not available, install with conda.
conda create -n vg python=3.10
conda activate vg
conda install vg=1.52 -c bioconda
ln -fs /path/to/miniconda3/envs/vg/bin/vg PanTax/tools
```

## Genome preprocessing

We recommend removing plasmids and redundancy from the genome first with `--remove` option and `--cluster` option, respectively. Eventually you will get a file containing genomic information in /path/to/database/library.

If genomes are all in NCBI refseq database, you only need to use `-r` option to specify the directory containing these genomes.

```
/path/to/PanTax/scripts/data_preprocessing -r ref --remove --cluster
```
Otherwise, you need to provide a file containing information about the custom genomes.
```
/path/to/PanTax/scripts/data_preprocessing -c genomes_info.txt --remove --cluster
```
The `genomes_info.txt` file gives a list of reference genomes in fasta format, which constitute PaxTax's original database, alongwith NCBI's taxonomic information. The input lines in the file should contain at least 5 tab-delimited fields; from left to right, they are Genome IDs, Strain taxonomic IDs, Species taxonomic IDs, Organism names, Genome absolute path.
Here is an example format of `genomes_info.txt` file:
```
genome_ID	strain_taxid	species_taxid	organism_name	id
GCF_000218545.1_ASM21854v1	593907	11	Cellulomonas gilvus ATCC 13127	/path/to/GCF_000218545.1_ASM21854v1_genomic.fna
GCF_025402875.1_ASM2540287v1	24.1	24	Shewanella putrefaciens	/path/to/GCF_025402875.1_ASM2540287v1_genomic.fna
```

## Running
* **Create database only** 
```
/path/to/PanTax/scripts/pantax -f $genome_info --create
```
You'll need to run /path/to/PanTax/scripts/pantax -f $genome_info --create. This will generate reference_pangenome.gfa and other files in your database directory.

Due to the large size of the reference pangenome we used for testing, we provide the `genomes_info.txt` used here. You need to download these genomes from NCBI RefSeq and update the actual paths in `genomes_info.txt`. Please note that NCBI RefSeq periodically updates their database, so we cannot guarantee that all the listed genomes will be available. Building the reference pangenome takes approximately one week with this `genomes_info.txt`. 

* **Query with specified database**
```
# long read
/path/to/PanTax/scripts/pantax -f $genome_info -l -r $fq -db $db --species-level --strain-level
# short read(pair-end)
/path/to/PanTax/scripts/pantax -f $genome_info -s -p -r $fq -db $db --species-level --strain-level
```

## options
```
Usage: /path/to/PanTax/scripts/pantax -f genomes_info -s/-l -r read.fq [option]
       paired-end: /path/to/PanTax/scripts/pantax -f genomes_info -s -p -r read.fq --species-level

Strain-level taxonomic classification of metagenomic data using pangenome graphs
    General options:
        --create                          Create database only.
        --vg FILE                         Path to vg executable file.
        -v, --verbose                     Detailed database build log.
        -t, --threads INT                 Number of processes to run in parallel (default: 64).
        --help, -h                        Print this help message.
        --version                         Print the version info.
    Database creation:
        -f, --genomesInformation FILE:    A list of reference genomes in specified format (Mandatory).
        -db NAME                          Name for pantax DB (default: pantax_db).
    Index construction:
        -s, --short-read                  Short read alignment.
        --best                            Best autoindex(only used with -s).
        -l, --long-read                   Long read alignment.
        --fast                            Fast index(only used with -l).
    Read classification:
        -r, --fastq-in FILE               Read and align FASTQ-format reads from FILE (two are allowed with -p).
        -p, --paired                      For paired-end alignment.
    Abundacnce calculation:
        --species-level                   Species abundance calulation.
        --strain-level                    Strain abundance calulation.
        -a float                          Species with more than abundance threshold used for strain abundance calulation.(default: 0.0001).
        -fr float                         Unique trio nodes fraction(default: 0.3).
        -fc float                         Unique trio nodes mean count fraction(default: 0.45).
        --filter                          MAPQ-based filter.
        --min_cov int                     Minimum coverage required per strain(default: 0).
        --min_depth int                   Graph nodes with sequence depth less than <min_depth>(default: 0).
        -gt int                           Gurobi threads(default: 1).
        -g, --save                        Save species graph information.
        -S, --classified-out FILENAME     File for alignment output(suffix).
        -R, --report FILENAME             File for read classification output(suffix).
        -o, --ouput FILENAME              File for abundance output(suffix).
```

## PanTax output
* **Alignment output**

You can specify the `-S` option to obtain this file, which is in GAF format, from [vg](https://github.com/vgteam/vg.git) or [Graphaligner](https://github.com/maickrau/GraphAligner). This file indicates the pangenome positions to which each read aligns. For details on the GAF format, please refer to the [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).

* **Classification output**(default output: pantax_report.tsv)

You can specify `-R` option to obtain this file. The following example shows classification assignments for a read. The assignment output has 4 columns.
```
read_id mapq    species_taxid   read_length
S0R1806	60	34	10654
The first column is the read ID from a raw sequencing read.
The second column is the mapping quality(mapq) for the read alignment, used to measure the quality of the alignment. (60 is the best)
The third column is the taxonomic ID of the read.
The fourth column is the length of the read.
```

* **Abundance profiling**
1. species level

The following example shows a classification summary at species level. 
```
species_taxid	predicted_abundance	predicted_coverage
34	0.5005489240249426	6.723225501680235

The first column is species taxonomic ID.
The second column is the abundance of this species normalized by its genomic length.
The third column is the average coverage of this species.
```

2. strain level

The following example shows a classification summary at strain level. 
```
species_taxid	strain_taxid	genome_ID	predicted_coverage	predicted_abundance
34	34.4	GCF_006401215.1_ASM640121v1	5.0	0.3945153945153945

The first column is species taxonomic ID.
The second column is strain taxonomic ID.
The third column is the name of a genome.
The second column is the abundance of this strain normalized by its genomic length.
The third column is the average coverage of this strain.
```

## Examples

* long read
```
cd PanTax/example/sample_hifi_data
# species level
sh ../../scripts/pantax -f ../sample_genome_info.txt -l -r long_reads.fq.gz --species-level
# strain level
sh ../../scripts/pantax -f ../sample_genome_info.txt -l -r long_reads.fq.gz --species-level --strain-level
```
* short read
```
cd PanTax/example/sample_ngs_data
# species level
sh ../../scripts/pantax -f ../sample_genome_info.txt -s -p -r short_reads.fq.gz --species-level
# strain level
sh ../../scripts/pantax -f ../sample_genome_info.txt -s -p -r short_reads.fq.gz --species-level --strain-level
```