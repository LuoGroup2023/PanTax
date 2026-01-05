
use serde::{Serialize, Deserialize};
use std::{collections::BTreeMap, path::PathBuf};
#[derive(Default, Debug)]
pub struct CheckPoints {
    pub reconstruction: bool,
    pub fast_query: bool,
    pub fast_query_and_filter: bool,
    pub need_index: bool,
    pub alignemnt: bool,
    // pub species_profiling: bool,
    // pub strain_profiling: bool,
    pub profiling: bool,
}

#[derive(Clone)]
pub struct GenomesInfo {
    pub genome_id: String,
    pub strain_taxid: String,
    pub species_taxid: String,
    pub organism_name: String,
    pub path_id: String,
}

#[derive(Debug, PartialEq)]
pub enum DataType {
    ShortReadSingle,
    ShortReadPaired,
    ShortReadPairedInter,
    LongReadSingle,
}

impl DataType {
    pub fn needs_index(&self) -> bool {
        match self {
            DataType::LongReadSingle => false,
            _ => true,
        }
    }    
}

// #[derive(Debug, Encode, Decode)]
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph {
    pub nodes_len: Vec<i64>,
    pub paths: BTreeMap<String, Vec<usize>>,
}

pub struct  ProfilingConfig {
    pub db: PathBuf,
    pub wd: PathBuf,
    pub genomes_metadata: Option<PathBuf>,
    pub range_file: Option<PathBuf>,
    pub input_aln_file: Option<PathBuf>,
    pub species_len_file: Option<PathBuf>,
    pub reads_binning_file: Option<PathBuf>,
    pub species_abund_file: Option<PathBuf>,  
    pub out_binning_file: Option<PathBuf>,
    pub min_species_abundance: f64, 
    pub unique_trio_nodes_fraction: f64,
    pub unique_trio_nodes_mean_count_f: f64,
    pub single_cov_ratio: f64,
    pub single_cov_diff: f64,
    pub minimization_min_cov: f64,
    pub min_cov: i64,
    pub min_depth: i64,
    pub species: bool,
    pub strain: bool,
    pub output_dir: PathBuf,
    pub solver: String,
    pub gurobi_threads: i32,
    pub sample_test: bool,
    pub sample_nodes: usize,
    pub designated_species: Option<String>,
    pub mode: Option<i32>,
    pub full: bool,
    pub shift: bool, 
    pub filtered: bool,
    pub zip: Option<String>,
    pub force: bool,
    pub trace: bool,
    pub debug: bool,
}