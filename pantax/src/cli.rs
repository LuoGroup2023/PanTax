
use crate::constants::{CLI_HEADINGS, DEFAULT_SOLVER};
use clap::{Parser, ValueEnum};
use std::path::PathBuf;


#[derive(Parser)]
#[command(
    author = "Wenhai Zhang", 
    version = "2.0.1", 
    about = "Strain-level metagenomic profiling using pangenome graphs with PanTax",
    arg_required_else_help = true
)]


#[derive(Default, Debug)]
pub struct Cli {

    // #[arg(num_args = 1..3, required = true, value_name = "FASTQ(.gz)")]
    // pub read_files: Vec<PathBuf>,

    /// Read and align FASTQ-format reads from FILE (two are allowed with -p)
    #[arg(short='r', long="reads", num_args = 1..3, help_heading = CLI_HEADINGS[0])]
    pub read_files: Option<Vec<PathBuf>>,    

    /// A list of reference genomes in specified format.
    #[arg(short = 'f', long = "genomesInformation", help_heading = CLI_HEADINGS[0])]
    pub genomes_info: Option<PathBuf>,

    /// Name for pantax database.
    #[arg(short = 'd', long, default_value = "pantax_db", help_heading = CLI_HEADINGS[0])]
    pub db: PathBuf, 

    /// Short read alignment.
    #[arg(short = 's', long, help_heading = CLI_HEADINGS[0])]
    pub short_read: bool,

    /// For paired-end alignment.
    #[arg(short = 'p', long, help_heading = CLI_HEADINGS[0])]
    pub paired: bool,

    /// Long read alignment.
    #[arg(short = 'l', long, help_heading = CLI_HEADINGS[0])]
    pub long_read: bool,

    /// Species abundance estimation.
    #[arg(long, visible_alias = "species-level", help_heading = CLI_HEADINGS[0])]
    pub species: bool,

    /// Strain abundance estimation.
    #[arg(long, visible_alias = "strain-level", help_heading = CLI_HEADINGS[0])]
    pub strain: bool,

    /// Number of processes to run in parallel.
    #[arg(short, long, default_value_t = default_threads(), help_heading = CLI_HEADINGS[0])]
    pub threads: usize,

    /// Create the database only.
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub create: bool,

    /// Create the database using genomes filtered by sylph query instead of all genomes.
    #[arg(long = "fast", help_heading = CLI_HEADINGS[1])]
    pub fast_query: bool,

    /// old paras, --mode 1 same as --fast. 
    #[arg(long, hide = true)]
    pub mode: Option<i32>,

    /// Specify Sylph syldb files. Sketching first can reduce runtime when processing multiple samples.
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub syldb: Option<String>,

    /// ANI threshold for sylph query result filter.
    #[arg(short = 'A', long, default_value_t = 99., help_heading = CLI_HEADINGS[1])]
    pub ani: f64,

    /// Save species graph information.
    #[arg(short = 'g', long, help_heading = CLI_HEADINGS[1])]
    pub save: bool,

    /// Serialized zip graph file saved with lz4 format (for save option).
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub lz: bool,    

    /// Serialized zip graph file saved with zstd format (for save option).
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub zstd: bool,      

    /// Create the index only. Make sure the pangenome construction is already complete and successful.
    #[arg(long, help_heading = CLI_HEADINGS[2])]
    pub index: bool,

    /// Vg autoindex for read alignment.
    #[arg(long, help_heading = CLI_HEADINGS[2])]
    pub auto: bool,

    /// Long read aligner (GraphAligner, vg >= 1.71). Use vg need to set long read type with --lt, vg only support hifi and r10.
    #[arg(long, default_value = "GraphAligner", help_heading = CLI_HEADINGS[3])]
    pub lr_aligner: PathBuf,

    /// Long read type (hifi, clr, ontr9, ontr10). Set precise clipping based on read type, and some empirical ANI for fast query, default is None.
    #[arg(long, visible_alias = "lt", help_heading = CLI_HEADINGS[3])]
    pub long_read_type: Option<String>,

    /// clip the alignment ends with arg as the identity cutoff between correct / wrong alignments.
    #[arg(long, default_value_t = 0.66, help_heading = CLI_HEADINGS[3])]
    pub precise_clipping: f64,

    /// fstrain. The fraction of strain-specific triplet nodes covered by reads for one strain. The larger, the stricter. (default: short 0.3/ long 0.5)
    #[arg(long = "fr", help_heading = CLI_HEADINGS[4])]
    pub unique_trio_nodes_fraction: Option<f64>,

    /// dstrain. The divergence between first rough and second refined strain abundance estimates. The smaller, the stricter. (default: 0.46)
    #[arg(long = "fc", help_heading = CLI_HEADINGS[4])]
    pub unique_trio_nodes_count: Option<f64>,

    /// Species with relative abundance above the threshold are used for strain abundance estimation.
    #[arg(short = 'a', default_value_t = 1e-04, help_heading = CLI_HEADINGS[4])]
    pub min_species_abundance: f64, 

    /// Rescued strain retention score.
    #[arg(long = "sr", default_value_t = 0.85, help_heading = CLI_HEADINGS[4])]
    pub single_cov_ratio: f64,

    /// Coverage depth difference between the strain and its species, with only one strain.
    #[arg(long = "sd", default_value_t = 0.2, help_heading = CLI_HEADINGS[4])]
    pub single_cov_diff: f64,   

    /// Shifting fraction of strain-specific triplet nodes. (multi-species: off, single-species: on) [true, false]
    #[arg(long, visible_alias = "sh", help_heading = CLI_HEADINGS[4])] 
    pub shift: Option<String>,

    /// Minimum coverage depth required per strain. 
    #[arg(long = "min_cov", default_value_t = 0, help_heading = CLI_HEADINGS[4])]
    pub min_cov: i64,

    /// Graph nodes with sequence coverage depth more than <min_depth>.
    #[arg(long = "min_depth", default_value_t = 0, help_heading = CLI_HEADINGS[4])]
    pub min_depth: i64,

    /// File for read classification(binning) output(prefix).
    #[arg(short = 'R', long = "report", help_heading = CLI_HEADINGS[4])]
    pub pantax_report: Option<PathBuf>,

    /// File for alignment output(prefix).
    #[arg(short = 'S', long = "classified-out", help_heading = CLI_HEADINGS[4])]
    pub read_aln: Option<String>,

    ///  File for abundance output(prefix).
    #[arg(short = 'o', long = "ouput", help_heading = CLI_HEADINGS[4])]
    pub pantax_output: Option<String>,

    /// MLP solver. (Gurobi, cplex, Cbc, highs, glpk)
    #[clap(long = "solver", default_value = DEFAULT_SOLVER, help_heading = CLI_HEADINGS[4])]
    pub solver: String,

    /// Solver threads.
    #[clap(long = "gthreads", default_value_t = 1, help_heading = CLI_HEADINGS[4])]
    pub gurobi_threads: i32,

    /// Temporary directory.
    #[arg(short = 'T', default_value = "pantax_db_tmp", help_heading = CLI_HEADINGS[5])]
    pub tmp_dir: PathBuf,

    /// Keep the temporary folder for later use at the strain level (resume).
    #[arg(short, long = "next", help_heading = CLI_HEADINGS[5])]
    pub next_for_strain: bool,

    /// Debug.
    #[arg(long, help_heading = CLI_HEADINGS[5])]
    pub debug: bool,

    /// Verbose.
    #[arg(short, long, help_heading = CLI_HEADINGS[5])]
    pub verbose: bool,

    /// Specifies a folder for the log files. 
    #[arg(long, help_heading = CLI_HEADINGS[6], hide = true)]
    pub log_dir: Option<PathBuf>,

    /// Add a marker String is added to the log file name. (e.g. sample name) 
    #[arg(long, help_heading = CLI_HEADINGS[6], hide = true)]
    pub log_m: Option<String>,        

    /// log level
    #[arg(long = "log", value_enum, default_value = "debug", help_heading = CLI_HEADINGS[6], hide = true)]
    pub log_level: LogLevel,

    /// only query and filter.
    #[arg(long = "qt", hide = true)]
    pub query_and_filter: bool,

    /// turn off multiple species parallel building.
    #[arg(long, hide = true)]
    pub no_parallel: bool,

    /// Force to rebuild pangenome.
    #[arg(long, hide = true)]
    pub force: bool,

    /// Reference genome for mc.
    #[arg(long, hide = true)]
    pub reference: Option<PathBuf>,

    /// Path to vg executable file.
    #[arg(long, default_value = "vg", hide = true)]
    pub vg: PathBuf,

    /// vg executable file version.
    #[arg(long, default_value = "v1", hide = true)]
    pub vg_lvl: VgLevel,

    /// Path to pangenome building tool executable file.
    #[arg(long, default_value = "pggb", hide = true)]
    pub pangenome_building_exe: PathBuf,   

    /// MAPQ-based filter (default: True/on).
    #[arg(long, hide = true)]
    pub no_filter: bool,    

    /// Sampling nodes(500000) are used for the graph with too many nodes.
    #[arg(long="sample", default_value_t = 500000, hide = true)]
    pub sample_nodes: usize,   

    /// Sampling nodes(500) are used for small model testing (set for codeocean).
    #[arg(long, hide = true)]
    pub sample_test: bool,  

    /// Only return the strain result of designated species. 
    #[arg(long="ds", hide = true)]
    pub designated_species: Option<String>,

    /// 0 (use all genomes), 1 (use sylph filter genomes, dataset specificity).
    #[arg(long, hide = true)]
    pub smode: Option<i32>,

    /// test, for save ori_strain_abundance.txt
    #[arg(long, hide = true)]
    pub test: bool,    
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum VgLevel {
    V1,
    V2,
}

impl Default for VgLevel {
    fn default() -> Self {
        VgLevel::V1
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Preset{
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl Default for LogLevel {
    fn default() -> Self {
        LogLevel::Debug
    }
}


impl Cli {
    pub fn log_level_filter(&self) -> log::LevelFilter {
        match self.log_level {
            LogLevel::Error => log::LevelFilter::Error,
            LogLevel::Warn => log::LevelFilter::Warn,
            LogLevel::Info => log::LevelFilter::Info,
            LogLevel::Debug => log::LevelFilter::Debug,
            LogLevel::Trace => log::LevelFilter::Trace,
        }
    }

    pub fn to_string(&self) -> String {
        format!("{:?}", self)
    }

    pub fn is_save(&self) -> bool {
        if self.save || self.lz || self.zstd {
            true
        } else {
            false
        }
    }
}

fn default_threads() -> usize {
    let t = (num_cpus::get() + 1) / 2;
    if t == 0 { 1 } else { t }
}