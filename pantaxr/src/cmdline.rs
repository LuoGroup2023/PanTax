
use clap::{Parser, Args, Subcommand};
use std::path::PathBuf;


#[derive(Parser)]
#[clap(author, version, arg_required_else_help = true, about = "pantaxr")]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}


#[derive(Subcommand)]
pub enum Mode {
    #[clap(arg_required_else_help = true, display_order = 1)]
    Profile(ProfileArgs),
    #[clap(arg_required_else_help = true, display_order = 2)]
    Fastixe(FastixeArgs),
    #[clap(arg_required_else_help = true, display_order = 3)]
    Rcls(RclsArgs),
    #[clap(arg_required_else_help = true, display_order = 4)]
    Stat(StatArgs),
}

#[derive(Args, Default, Debug)]
pub struct ProfileArgs {
    #[clap(short = 'd', long = "db", help_heading = "INPUT", help = "PanTax database directory.")]
    pub db: PathBuf, 

    #[clap(short = 'w', long = "wd", default_value = "./", help_heading = "INPUT", help = "PanTax work directory.")]
    pub wd: PathBuf, 

    #[clap(short = 'g', long = "metadata", help_heading = "INPUT", help = "Genomes information(metadata) file. (Optional, default included in db)")]
    pub genomes_metadata: Option<PathBuf>,    

    #[clap(short = 's', long = "range", help_heading = "INPUT", help = "Species range file. (Optional, default included in db)")]
    pub range_file: Option<PathBuf>,

    #[clap(short = 'm', long = "gaf", help_heading = "SPECIES PROFILE INPUT FILE", help = "Mapping gaf file.")]
    pub input_aln_file: Option<PathBuf>,

    #[clap(short = 'l', long = "len", help_heading = "SPECIES PROFILE INPUT FILE", help = "Average genome length of species file. (Optional, default included in db)")]
    pub species_len_file: Option<PathBuf>,

    #[clap(short = 'b', long = "bin", help_heading = "STRAIN PROFILE INPUT FILE", help = "Read binning file from species profiling.")]
    pub reads_binning_file: Option<PathBuf>,    

    #[clap(short = 'a', long = "abund", help_heading = "STRAIN PROFILE INPUT FILE", help = "Species abundance file from species profiling.")]
    pub species_abund_file: Option<PathBuf>,  

    #[clap(long = "binning_out", help_heading = "SPECIES PROFILE PARAS OPTIONS", help = "Reads binning file output.")]
    pub out_binning_file: Option<PathBuf>,      

    #[clap(long = "min_species_abundance", default_value_t = 1e-04, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Minimum species abundance(filter low abundance species).")]
    pub min_species_abundance: f64,    

    #[clap(short = 'f', long = "unique_trio_nodes_fraction", default_value_t = 0.3, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Unique trio nodes fraction.")]
    pub unique_trio_nodes_fraction: f64,

    #[clap(short = 'r', long = "unique_trio_nodes_mean_count_f", default_value_t = 0.46, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Unique trio nodes mean count fraction.")]
    pub unique_trio_nodes_mean_count_f: f64,

    #[clap(long = "scr", default_value_t = 0.85, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Single species strain level coverage ratio (Identify strains with high confidence).")]
    pub single_cov_ratio: f64,    

    #[clap(long = "scd", default_value_t = 0.2, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Coverage difference for the species with single strain.")]
    pub single_cov_diff: f64, 

    #[clap(long = "mmcov", default_value_t = 0.0, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Strain min coverage.")]
    pub minimization_min_cov: f64,

    #[clap(short = 'c', long = "min_cov", default_value_t = 0, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Minimum coverage required per strain.")]
    pub min_cov: i64,    

    #[clap(long = "min_depth", default_value_t = 0, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Graph nodes with sequence coverage depth more than <min_depth>.")]
    pub min_depth: i64,

    #[clap(long="species", help_heading = "GENERAL OPTIONS", help = "Perform species level profiling.")]
    pub species: bool,

    #[clap(long="strain", help_heading = "GENERAL OPTIONS", help = "Perform strain level profiling.")]
    pub strain: bool,

    #[clap(short = 'o', long = "out", default_value = "./", help_heading = "GENERAL OPTIONS", help = "Output directory.")]
    pub output_dir: PathBuf,

    #[clap(short = 't', long = "threads", default_value_t = 1, help_heading = "GENERAL OPTIONS", help = "Number of threads.")]
    pub threads: usize, 

    #[clap(long = "gthreads", default_value_t = 1, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Number of gurobi threads.")]
    pub gurobi_threads: i32,

    #[clap(long="sample_test", help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Sampling nodes are used for small model testing.")]
    pub sample_test: bool,    

    #[clap(long="sample", default_value_t = 500000, help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Sampling nodes are used for the graph with too many nodes.")]
    pub sample_nodes: usize,

    #[clap(long="ds", help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Designated species.")]
    pub designated_species: Option<String>,

    #[clap(long="mode", help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Select species to perform strain profiling [0: the species with one strains, 1: the species with more strains].")]
    pub mode: Option<i32>,

    #[clap(long="full", help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Retain full float precision (default: 2 decimals).")]
    pub full: bool,

    #[clap(long="shift", help_heading = "STRAIN PROFILE PARAS OPTIONS", help = "Unique_trio_nodes_fraction varies with strain coverage depth when estimating single species.")]
    pub shift: bool,    

    #[clap(long="filtered", help_heading = "SPECIES PROFILE PARAS OPTIONS", help = "MAPQ-based filtered.")]
    pub filtered: bool,
    
    #[clap(long="force", help = "Force rerun.")]
    pub force: bool,

    #[clap(long="trace", help = "Trace output (caution: very verbose).")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output.")]
    pub debug: bool,
}

#[derive(Args, Default, Debug)]
pub struct RclsArgs {

    #[clap(short = 'm', long = "gaf", help_heading = "INPUT FILE", help = "Mapping gaf file.")]
    pub input_aln_file: PathBuf,

    #[clap(short = 's', long = "range", help_heading = "INPUT FILE", help = "Species range file.")]
    pub range_file: PathBuf,

    #[clap(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
    pub threads: usize,

    #[clap(short = 'o', long = "out", help_heading = "OUTPUT", help = "Output file name.")]
    pub output_name: PathBuf,

    #[clap(long="lazy", help = "Lazy reading.")]
    pub lazy_read: bool,

    #[clap(long="struct", help = "Struct reading.")]
    pub struct_read: bool
    
}

#[derive(Args, Default, Debug)]
pub struct FastixeArgs {

    #[clap(short = 's', long = "input-files", help_heading = "INPUT FILE", num_args = 1.., help = "Multiple input files.")]
    pub input_files: Option<Vec<PathBuf>>,

    #[clap(short = 'l', long = "input-genome-list", help_heading = "INPUT FILE", help = "Input genome list.")]
    pub input_list: Option<PathBuf>,

    #[clap(short = 'd', long = "input-dir", help_heading = "INPUT FILE", help = "Input directory containing FASTA files.")]
    pub input_directory: Option<PathBuf>,

    #[clap(short = 'o', long = "out-dir", default_value = "genomes", help_heading = "OUTPUT", help = "Output directory.")]
    pub out_directory: PathBuf,

    #[clap(short = 'p', long = "prefix", help_heading = "RENAME", help = "Prefix to add to headers.")]
    pub prefix: Option<String>,

    #[clap(short = 'r', long = "regex", default_value_t = String::from("[^_]+_[^_]+"), help_heading = "RENAME", help = "File name regex")]
    pub reg: String,

    #[clap(short, long="up", help_heading = "Sequence", help = "All bases are converted to uppercase letters.")]
    pub uppercase: bool,

    #[clap(short = 'e', long="output-file-name", default_value_t = String::from("merged.fa"), help_heading = "MERGE OUTPUT", help = "Merge output file path.")]
    pub merge_output_file_path: String,

    #[clap(short = 'b', long="bgz", help_heading = "MERGE OUTPUT", help = "Merge bgzip output.")]
    pub merge_bgzip_output: bool,    

    #[clap(long = "level", help_heading = "COMPRESSION LEVEL", help = "Compression (0-9).")]
    pub compression_level: Option<u32>,

    #[clap(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
    pub threads: usize,

    #[clap(long="trace", help = "Trace output (caution: very verbose).")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output.")]
    pub debug: bool,
}

#[derive(Args, Default, Debug)]
pub struct StatArgs {
    #[clap(short = 'i', long = "input_file", help_heading = "INPUT FILE", help = "Input genome file.")]
    pub input_file: Option<PathBuf>,

    #[clap(short = 'l', long = "input_genome_list", help_heading = "INPUT FILE", help = "Input genome list.")]
    pub input_list: Option<PathBuf>,  

    #[clap(short = 'g', long = "genome_info", help_heading = "INPUT FILE", help = "Input genome metadata.")]
    pub genome_metadata_file: Option<PathBuf>,

    #[clap(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
    pub threads: usize,

    #[clap(short = 'o', long = "out", default_value = "genome_statics.txt", help_heading = "OUTPUT", help = "Output file name.")]
    pub output_name: PathBuf,

}