# PanTax: Strain-level metagenomic profiling using pangenome graphs

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/pantax.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/pantax)
[![License](https://img.shields.io/github/license/LuoGroup2023/PanTax)](https://www.gnu.org/)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/LuoGroup2023/PanTax)](https://github.com/LuoGroup2023/PanTax/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/LuoGroup2023/PanTax/total.svg?style=social&logo=github&label=Download)](https://github.com/LuoGroup2023/PanTax/releases)
[![GitHub stars](https://img.shields.io/github/stars/LuoGroup2023/PanTax?color=yellow)](https://github.com/LuoGroup2023/PanTax/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/LuoGroup2023/PanTax?color=lightblue&label=fork)](https://github.com/LuoGroup2023/PanTax/network/members) 
<!-- [![visitors](https://visitor-badge.laobi.icu/badge?page_id=LuoGroup2023.PanTax)](https://github.com/LuoGroup2023/PanTax) -->


## Table of Contents

+ [Overview](#overview)
+ [Installation](#installation)
+ [Installation (v2.0.0)](#installation-version-200)
+ [Genome preprocessing](#genome-preprocessing)
+ [Running](#running)
+ [Options](#options)
+ [PanTax output](#pantax-output)
+ [Examples](#examples)
+ [Possible issues during installation](#possible-issues-during-installation)
+ [Change](#change)
+ [TODO](#todo)
+ [Citation](#citation)

> [!IMPORTANT]
> **Pantax v2.0.0 dev is released in the rust_dev branch.**
> 
> The species and strain profiling module, the conversion of a single strain into a GFA module, the graph node and path information serialization and compression module, and the long read GAF filtering module are all rewritten using Rust. They are all integrated into the subcommand of pantaxr. It runs more than 25x faster in profiling.

## Overview

PanTax is a pangenome graph-based taxonomic profiling tool designed for accurate strain-level classification of metagenomic sequencing data. Unlike traditional methods that rely on multiple linear reference genomes, PanTax leverages pangenome graphs to better represent genetic variation and relationships across related genomes. It supports both short and long reads, works across single or multiple species, and delivers superior precision or recall at the strain level. PanTax provides a robust, scalable solution to taxonomic classification for strain resolution, overcoming key limitations of existing tools.

## Installation
The dependencies of PaxTax can generally be installed through conda. There are two dependencies that require special attention. 

Firstly, PanTax relies on the [Gurobi Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer/) Python package. While the package download is included in the environment.yaml file, a license must be obtained from the official website to use it for large model optimizations.

Secondly, the [pggb](https://github.com/pangenome/pggb.git) on Bioconda depends on a specific version of [vg](https://github.com/vgteam/vg.git) 1.59, preventing timely updates to the latest vg version. This implies that installing PanTax from Bioconda would face the same issue. Therefore, we recommend installing PanTax from the source. You need create a new environment with conda and installing latest vg, and then place it in the tools directory via a symbolic link or specify its path directly using the `--vg` option.

* **From bioconda**

```
conda install -c bioconda -c conda-forge -c gurobi pantax gurobi 

## Run pantax.
pantax -h
```

* **From source**
```
git clone https://github.com/LuoGroup2023/PanTax.git --recursive
conda create -n pantax python=3.10
conda activate pantax
conda install -c bioconda -c conda-forge -c gurobi -c defaults \
    r-base=4.2 \
    samtools=1.19.2 \
    bcftools=1.19 \
    htslib=1.19.1 \
    pggb \
    vg \
    graphaligner \
    sylph \
    h5py \
    pandas \
    tqdm \
    numpy \
    networkx \
    pyarrow \
    gurobi
cd PanTax
bash install.sh

# If vg is not available, install with conda.(optional)
conda create -n vg python=3.10
conda activate vg
conda install vg -c bioconda
cd tools
ln -fs $(which vg) ./

# Run pantax
cd ../scripts
./pantax -h
```

* **From source (Version for conducting experiments)**

> The version of the main tools we used in the paper's experiments are `pggb` 0.6.0, `vg` 1.52, `sylph` 0.6.1, `gurobipy` 11.0.2.

```
git clone https://github.com/LuoGroup2023/PanTax.git
conda env create -f environment.yaml
conda activate pantax_v12
cd PanTax
bash install.sh

# If vg is not available, install with conda.(optional)
conda create -n vg python=3.10
conda activate vg
conda install vg=1.52 -c bioconda
cd tools
ln -fs $(which vg) ./

# Run pantax
cd ../scripts
chmod +x pantax pantax_utils data_preprocessing 
./pantax -h
```

## Installation (Version 2.0.0)
* **From bioconda**
Bioconda installation is not supported for now.

* **From source**
```
git clone https://github.com/LuoGroup2023/PanTax.git -b rust_dev
conda create -n pantax
conda activate pantax
conda install -c bioconda -c conda-forge -c gurobi -c defaults \
    python=3.10 \
    r-base=4.2 \
    pggb=0.6.0 \
    vg \
    graphaligner=1.0.17 \
    sylph=0.6.1 \
    fastani=1.33 \
    pandas \
    tqdm \
    numpy \
    networkx \
    pyarrow \
    gurobi=11 \
    clang \
    rust=1.82
cd PanTax
bash install.sh v2

# If vg is not available, install with conda.(optional)
conda create -n vg python=3.10
conda activate vg
conda install vg=1.59 -c bioconda
cd tools
ln -fs $(which vg) ./

# Run pantax
cd ../scripts
./pantax -h
```

You may also choose not to specify the version of the tool, but the impact of using the latest version has not yet been tested.

## Genome preprocessing

We recommend removing plasmids and redundancy from the genome first with `--remove`, `--compute`, `--cluster` option. Eventually you will get final file containing genomic information in /path/to/database.

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
pantax -f $genome_info --create
```
You'll need to run `/path/to/PanTax/scripts/pantax -f $genome_info --create`. This will generate reference_pangenome.gfa and other files in your database directory.

Due to the large size of the reference pangenome we used for testing, we provide the `genomes_info.txt` used here. You need to download these genomes from NCBI RefSeq and update the actual paths in `genomes_info.txt`. Please note that NCBI RefSeq periodically updates their database, so we cannot guarantee that all the listed genomes will be available. Building the reference pangenome takes approximately one week with this `genomes_info.txt`. 

* **Query with specified database**
1. species level and strain level 
```
# long read
pantax -f $genome_info -l -r $fq -db $db --species-level --strain-level
# short read(pair-end)
pantax -f $genome_info -s -p -r $fq -db $db --species-level --strain-level
```
2. species level and then strain level
```
# long read
# species level 
pantax -f $genome_info -l -r $fq -db $db --species-level -n
# strain level
pantax -f $genome_info -l -r $fq -db $db --strain-level -n
```

* **Fast query with specified database**
```
# long read
pantax -f $genome_info -l -r $fq -db $db --species-level --strain-level --fast
# short read(pair-end)
pantax -f $genome_info -s -p -r $fq -db $db --species-level --strain-level --fast
```


## Options
```
Usage: /path/to/PanTax/scripts/pantax -f genomes_info -s/-l -r read.fq [option]
       paired-end: /path/to/PanTax/scripts/pantax -f genomes_info -s -p -r read.fq --species-level --strain-level

Strain-level metagenomic profiling using pangenome graphs with PanTax
    General options:
        -f, --genomesInformation file:    A list of reference genomes in specified format. (Mandatory)
        -db dir                           Name for pantax DB (default: pantax_db).
        -T dir                            Temporary directory (default: pantax_db_tmp).
        -n, --next                        Keep the temporary folder for later use at the strain level (resume).
        -t, --threads int                 Number of processes to run in parallel. (default: all available)
        -v, --verbose                     Detailed database build log.
        --vg file                         Path to vg executable file.
        --debug                           Keep the temporary folder for any situation.
        --help, -h                        Print this help message.
        --version                         Print the version info.
    Database creation:
        --create                          Create the database only.
        --fast                            Create the database using genomes filtered by sylph query instead of all genomes.
        -g, --save                        Save species graph information.
            --lz                          Serialized zip graph file saved with lz4 format (for save option).
            --zstd                        Serialized zip graph file saved with zstd format (for save option).
        --force                           Force to rebuild pangenome.
        -e file                           Path to pangenome building tool executable file. (default: pggb)
        -A, --ani float                   ANI threshold for sylph query result filter. (default: 99)
    Index construction(for vg giraffe):
        --index                           Create the index only.
        --best                            Best autoindex, which corresponds to vg autoindex. (only used with -s)
        --fast-aln                        Long read fast alignment with vg instead of Graphaligner. (only used with -l)
    Read alignment:
        -r, --fastq-in file               Read and align FASTQ-format reads from FILE (two are allowed with -p).
        -s, --short-read                  Short read alignment.
        -p, --paired                      For paired-end alignment.
        -l, --long-read                   Long read alignment.
        -lt, --long-read-type str         Long read type (hifi, clr, ontr9, ontr10). Set precise-clipping based on read type.
        --precise-clipping float          clip the alignment ends with arg as the identity cutoff between correct / wrong alignments. (default: 0.66)
    Abundacnce calculation:
        --species-level                   Species abundance calulation.
        --strain-level                    Strain abundance calulation.
        -a float                          Species with relative abundance above the threshold are used for strain abundance estimation. (default: 0.0001)
        -fr float                         The fraction of strain-specific triplet nodes covered by reads for one strain. The larger, the stricter. (default: multi species 0.3/ single species 0.43)
        -fc float                         The divergence between first rough and second refined strain abundance estimates. The smaller, the stricter. (default: 0.46)
        -sh, --shift bool                 Shifting fraction of strain-specific triplet nodes. (multi-species: off, single-species: on)
        -sd float                         Coverage depth difference between the strain and its species, with only one strain. (default: 0.2)
        -sr float                         Rescued strain retention score. (default: 0.85)
        --min_cov int                     Minimum coverage depth required per strain. (default: 0)
        --min_depth int                   Graph nodes with sequence coverage depth more than <min_depth>. (default: 0)
        -gt int                           Gurobi threads. (default: 1)
        -S, --classified-out str          File for alignment output(prefix).
        -R, --report str                  File for read classification(binning) output(prefix).
        -o, --ouput str                   File for abundance output(prefix).
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
species_taxid	strain_taxid	genome_ID	predicted_coverage	predicted_abundance	path_base_cov	unique_trio_fraction	uniq_trio_cov_mean	first_sol	strain_cov_diff	total_cov_diff
34	34.4	GCF_006401215.1_ASM640121v1	16.0	0.39983790355261384	0.9967217217217217	1.0	15.54	16.0	0.01	0.0010005002501250622

The first column is species taxonomic ID.
The second column is strain taxonomic ID.
The third column is the name of a genome (assembly accession).
The fourth column is the average coverage depth of this strain.
The fifth column is the relative abundance of this strain in the sample, which normalized by its genomic length.
The sixth column is the coverage of the strain.
The seventh column is the coverage of strain-specific triplet nodes.
The eighth column is the average abundance of all strain-specific triplet nodes for this strain.
The ninth column is the solution of the first path abundance optimization, which represents the average coverage depth of the strain in the first solution.
The tenth column reflects the divergence between the values in the ninth and tenth columns.
The eleventh column represents the difference between the sum of the average coverage depths of all strains of the species and the average coverage depth of the species.
```

## Examples

* long read
```
cd PanTax/example/hifi
# species level
pantax -f ../example_genomes_info.txt -l -r long_reads.fq.gz --species-level -n
# strain level
pantax -f ../example_genomes_info.txt -l -r long_reads.fq.gz --strain-level -n
```
* short read
```
cd PanTax/example/ngs
# species level
pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --species-level -n
# strain level
pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --strain-level -n
```

## Possible issues during installation
* `gsl` incompatibility

During the "Building reference pangenome" step, there were no errors reported by `pantax` and an interrupt occurred without obtaining any results. You should run the command `wfmash -h` to check whether `wfmash` can work. If you encounter a error `symbol lookup error: ../lib/libgsl.so.25: undefined symbol: cblas_ctrmv` after running it, this indicates a GSL incompatibility issue. The `gsl` should be installed from the `conda-forge` channel `https://anaconda.org/conda-forge/gsl/files` instead of the `anaconda` channel `https://anaconda.org/anaconda/gsl/files`. This issue is caused by the channel settings of conda, so you should specify the channel with the `-c` flag during installation.
```
conda install pantax -c bioconda -c conda-forge -c gurobi
```

## Change

### Version: V2.0.0 dev (update at 2024-06-14)

<details>
<summary>Click here to check the log of all updates</summary>

#### *__[Update - 2024 - 04 - 08]__* :  <BR/>

*V0.0.1: Initial commit! <BR/>*

#### *__[Update - 2024 - 09 - 09]__* :  <BR/>

*V1.0.1 <BR/>*
* *Fix a bug caused by the presence of the same path or a total of less than 3 nodes per path in GFA <BR/>*
* *database creation, species level and strain level work step by step now <BR/>*
* *Optimization of strain level calculation <BR/>*
* *Sampling nodes(500) are used for small model testing (set for codeocean) <BR/>*
* *Improve the logger for strain level calculation <BR/>*

#### *__[Update - 2024 - 10 - 09]__* :  <BR/>

*V1.0.2 <BR/>*
* *Fix a bug in numpy matrix caused by the graph has the same path or the total number of nodes per path less than 3 <BR/>*
* *Kill the long-term vg alignment and then continue the analysis <BR/>*

#### *__[Update - 2024 - 10 - 12]__* :  <BR/>

* *PanTax can be downloaded from bioconda now <BR/>*

#### *__[Update - 2024 - 11 - 09]__* :  <BR/>

*V1.1.0 <BR/>*
* *Estimate abundance with GAF file, no longer json file.<BR/>*
* *Fix the bug when using GAF file to estimate abundance. In the GAF file, the read_end columns may be * instead of a number.<BR/>*

#### *__[Update - 2025 - 04 - 03]__* :  <BR/>

*V1.2.0 <BR/>*
* *PanTax can accept interleaved and un-interleaved paired-end reads now.<BR/>*
* *Implementing pangenome construction task scheduling to accelerate construction.<BR/>*
* *Adding hierarchical clustering methods for single species redundancy reduction and code simplification.<BR/>*
* *Add the --fast option to execute the fast version of PanTax. This version utilizes the ultrafast species-level metagenomic profiler sylph for rapid       species/strain filtering, significantly reducing the number of species-level pan-genome constructions for metagenomic data.<BR/>*
* *PanTax adds a strain rescue step.<BR/>*
* *For single-species strain-level tasks, the -fr option is set as a dynamic metric, becoming stricter as the average coverage depth of the strain increases.<BR/>*
* *The strain profiling result file is now more detailed. For specifics, please refer to the description in the new README file.<BR/>*
* *PanTax is now compatible with strain profiling using the GTDB database.<BR/>*
* *data_preprocessing is applicable for deduplication filtering using the GTDB database with --db gtdb.<BR/>*
* *data_preprocessing now can use hierarchical clustering methods for single species redundancy reduction，and optimize the data_preprocessing code.<BR/>*
* *The genome information file's id column can now include genomes in gzip-compressed (.gz) file formats.<BR/>*

</details>

#### *__[Update - 2025 - 06 - 14]__* :  <BR/>
*V2.0.0 dev <BR/>*
* *Fix bug: filter a few reads. The mapping of very few reads in GAF shows end > start, and only one node of these reads maps to the graph. <BR/>*
* *The species and strain profiling module, the conversion of a single strain into a GFA module, the graph node and path information serialization and compression module, and the long read GAF filtering module are all rewritten using Rust. They are all integrated into the subcommand of pantaxr. It runs more than 25x faster in profiling. In general, the results will change slightly, with higher precision and lower recall. <BR/>*



## TODO
+ Performance comparison between `vg giraffe` long read alignment and `Graphaligner` alignment.

## Citation
```
@article{luo2025strain,
  title={Strain-level metagenomic profiling using pangenome graphs with PanTax},
  author={Luo, Xiao and Zhang, Wenhai and Liu, Yuansheng and Li, Guangyi and Xu, Jialu and Chen, Enlian and Schonhuth, Alexander},
  journal={bioRxiv},
  pages={2025--04},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```