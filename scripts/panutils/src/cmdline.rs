
use clap::{Parser, Args, Subcommand};
use std::path::PathBuf;


#[derive(Parser)]
#[clap(author, version, arg_required_else_help = true, about = "panutils")]
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

}

#[derive(Args, Default, Debug)]
pub struct ProfileArgs {
    #[clap(short = 'm', long = "gaf", help_heading = "INPUT FILE", help = "Mapping gaf file.")]
    pub input_aln_file: PathBuf,

    #[clap(short = 's', long = "range", help_heading = "INPUT FILE", help = "Species range file.")]
    pub range_file: PathBuf,   

    #[clap(short = 'l', long = "len", help_heading = "SPECIES PROFILE INPUT FILE", help = "Average genome length of species file.")]
    pub species_len_file: PathBuf,     

    #[clap(short = 'd', long = "db", help_heading = "INPUT FILE", help = "PanTax database directory.")]
    pub db: PathBuf,      

    #[clap(short = 't', long = "threads", default_value_t = 1, help = "Number of threads [default: 1].")]
    pub threads: usize, 


    #[clap(short = 'o', long = "out", help_heading = "OUTPUT", help = "Output file name.")]
    pub output_name: PathBuf,

    #[clap(long="ds", help = "Designated species .")]
    pub designated_species: Option<String>,

    #[clap(long="filtered", help = "MAPQ-based filtered.")]
    pub filtered: bool,

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

    #[clap(short = 't', long = "threads", default_value_t = 1, help = "Number of threads [default: 1].")]
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

    #[clap(short = 't', long = "threads", default_value_t = 1, help = "Number of threads [default: 1].")]
    pub threads: usize,

    #[clap(long="trace", help = "Trace output (caution: very verbose).")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output.")]
    pub debug: bool,
}