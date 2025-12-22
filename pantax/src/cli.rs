
use crate::constants::CLI_HEADINGS;
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
    #[arg(short, long, default_value_t = num_cpus::get(), help_heading = CLI_HEADINGS[0])]
    pub threads: usize,

    /// Create the database only.
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub create: bool,

    /// Create the database using genomes filtered by sylph query instead of all genomes.
    #[arg(long = "fast", visible_alias = "mode", help_heading = CLI_HEADINGS[1])]
    pub fast_query: bool,

    /// ANI threshold for sylph query result filter.
    #[arg(short = 'A', long, default_value_t = 99., help_heading = CLI_HEADINGS[1])]
    pub ani: f64,

    /// Create the index only.
    #[arg(long, help_heading = CLI_HEADINGS[2])]
    pub index: bool,

    /// Long read type (hifi, clr, ontr9, ontr10). Set precise clipping based on read type, and some empirical ANI for fast query.
    #[arg(long, visible_alias = "lt", help_heading = CLI_HEADINGS[3])]
    pub long_read_type: Option<String>,

    /// Temporary directory.
    #[arg(short = 'T', default_value = "pantax_db_tmp", help_heading = CLI_HEADINGS[4])]
    pub tmp_dir: PathBuf,

    /// Keep the temporary folder for later use at the strain level (resume).
    #[arg(short, long = "next", help_heading = CLI_HEADINGS[4])]
    pub next_for_strain: bool,

    /// Debug.
    #[arg(long, help_heading = CLI_HEADINGS[4])]
    pub debug: bool,

    /// Specifies a folder for the log files. 
    #[arg(long, help_heading = CLI_HEADINGS[5])]
    pub log_dir: Option<PathBuf>,

    /// Add a marker String is added to the log file name. (e.g. sample name) 
    #[arg(long, help_heading = CLI_HEADINGS[5])]
    pub log_m: Option<String>,        

    /// log level
    #[arg(long = "log", value_enum, default_value = "debug", help_heading = CLI_HEADINGS[5])]
    pub log_level: LogLevel,

    /// only query and filter
    #[arg(long = "qt", hide = true)]
    pub query_and_filter: bool,
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
}


