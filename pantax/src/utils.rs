
use crate::types::*;
use crate::cli::Cli;
use crate::constants::GENOMES_INFO;
use std::process::Command;
use std::io::{Write, BufWriter};
use std::fs::File;
use std::path::PathBuf;
use anyhow::{Result, Context};
use std::sync::OnceLock;

/// Some custom directory path and file path
pub struct GlobalConfig {
    pub genomes_metadata: PathBuf,
    pub query_res: PathBuf,
    pub query_done: PathBuf,
    pub pangeno_dir: PathBuf,
    pub single_species_file: PathBuf,
    pub gfa_build_wd: PathBuf,
    pub gfa_done_path: PathBuf,
    pub zip_wd: PathBuf,
    pub range_file: PathBuf,
    pub reference_pangenome_vg: PathBuf,
    pub reference_pangenome_gfa: PathBuf,
    pub species_graph_info: PathBuf,
    pub species_gfa: PathBuf,
    pub reference_pangenome_gfa_db: PathBuf,
    pub range_file_db: PathBuf,
    pub species_genomes_stats: PathBuf,
    pub index_files: Vec<PathBuf>,
    pub gfa_mapped: PathBuf,
    pub pantax_report: Option<PathBuf>,
    pub species_abund_file: PathBuf,
    pub strain_abund_file: PathBuf,
    pub ori_strain_abund_file: PathBuf,
}

static GLOBAL_CONFIG: OnceLock<GlobalConfig> = OnceLock::new();

pub fn init_global_config(args: &Cli) {
    let query_res = args.tmp_dir.join("sylph_db/query_result");
    let query_done = args.tmp_dir.join("sylph_db/sylph_done");
    let pangeno_dir = args.tmp_dir.join("species_pangenome");
    let single_species_file = args.tmp_dir.join("species_eq1_genomes.txt");
    let gfa_build_wd = args.tmp_dir.join("gfa_build");
    let gfa_done_path = gfa_build_wd.join("gfa_build_done.txt");
    let zip_wd = PathBuf::from(format!("{}2", gfa_build_wd.display()));
    let range_file = args.tmp_dir.join("species_range.txt");
    let reference_pangenome_vg = args.tmp_dir.join("reference_pangenome.vg");
    let reference_pangenome_gfa = args.tmp_dir.join("reference_pangenome.gfa");
    let gfa_mapped = args.tmp_dir.join("gfa_mapped.gaf");
    let pantax_report = if args.pantax_report.is_none() && (args.next_for_strain || args.debug) {
        Some(args.tmp_dir.join("reads_classification.tsv"))
    } else if args.pantax_report.is_some() {
        Some(args.pantax_report.as_ref().unwrap().with_extension("tsv"))
    } else {
        None
    };
    let species_abund_file = args.tmp_dir.join("species_abundance.txt");
    let strain_abund_file = args.tmp_dir.join("strain_abundance.txt");
    let ori_strain_abund_file = args.tmp_dir.join("ori_strain_abundance.txt");

    let genomes_metadata = args.db.join(GENOMES_INFO);
    let species_graph_info = args.db.join("species_graph_info");
    let species_gfa = args.db.join("species_gfa");
    let reference_pangenome_gfa_db = args.db.join("reference_pangenome.gfa");
    let range_file_db = args.db.join("species_range.txt");
    let species_genomes_stats = args.db.join("species_genomes_stats.txt");
    let reference_pangenome_min = args.db.join("reference_pangenome.min");
    let reference_pangenome_giraffe_gbz = args.db.join("reference_pangenome.giraffe.gbz");
    let reference_pangenome_dist = args.db.join("reference_pangenome.dist");
    
    let config = GlobalConfig {
        genomes_metadata,
        query_res,
        query_done,
        pangeno_dir,
        single_species_file,
        gfa_build_wd,
        gfa_done_path,
        zip_wd,
        range_file,
        reference_pangenome_vg,
        reference_pangenome_gfa,
        species_graph_info,
        species_gfa,
        reference_pangenome_gfa_db,
        range_file_db,
        species_genomes_stats,
        index_files: vec![reference_pangenome_giraffe_gbz, reference_pangenome_dist, reference_pangenome_min],
        gfa_mapped,
        pantax_report,
        species_abund_file,
        strain_abund_file,
        ori_strain_abund_file,
    };
    let _res = GLOBAL_CONFIG.set(config);
}

pub fn get_config() -> &'static GlobalConfig {
    GLOBAL_CONFIG.get().expect("GlobalConfig not initialized. Call init_global_config first.")
}


pub fn output_genomes_info(genomes_info: &Vec<GenomesInfo>, genomes_info_path: &PathBuf) -> Result<()> {
    let file = File::create(genomes_info_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid")?;
    for genome_info in genomes_info.iter() {
        writeln!(writer, "{}\t{}\t{}\t{}\t{}", genome_info.genome_id, genome_info.strain_taxid, genome_info.species_taxid, genome_info.organism_name, genome_info.path_id)?;
    }
    Ok(())
}

pub fn run_shell_command(command: &str) -> Result<()> {        
    let output = Command::new("bash")
        .arg("-c")
        .arg(command)
        .output()
        .with_context(|| format!("Failed to execute command: {}", command))?;
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        log::error!("Command failed with status: {:?}", output.status);
        log::error!("Command: {}", command);
        if !stderr.trim().is_empty() {
            log::error!("Stderr: {}", stderr);
        }
        if !stdout.trim().is_empty() {
            log::error!("Stdout: {}", stdout);
        }
        
        return Err(anyhow::anyhow!(
            "Command failed with exit code: {}\nCommand: {}",
            output.status.code().unwrap_or(-1),
            command
        ));
    }

    let stderr = String::from_utf8_lossy(&output.stderr);
    // let _stdout = String::from_utf8_lossy(&output.stdout);
    
    if !stderr.trim().is_empty() {
        log::warn!("Command succeeded with stderr output:\n{}", stderr);
    }
       
    Ok(())
}