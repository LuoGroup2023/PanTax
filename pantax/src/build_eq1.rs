
use crate::cli::Cli;
use crate::utils::get_config;
use crate::utils::run_shell_command;
use crate::zip::{CompressType, zip};
use fxhash::FxHashMap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::fs::{File, create_dir_all, remove_file};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use anyhow::{anyhow, Result};

fn read_fasta(path: &str) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let mut records = Vec::new();
    let mut reader = parse_fastx_file(path)?;
    while let Some(record) = reader.next() {
        let seqrec = record?;
        let id = String::from_utf8_lossy(seqrec.id()).to_string();
        let seq = String::from_utf8_lossy(&seqrec.seq()).to_string();
        records.push((id, seq));
    }
    Ok(records)
}

fn chunk_sequence(seq: &str, chunk_size: usize) -> Vec<&str> {
    let mut chunks = Vec::new();
    let mut start = 0;
    let len = seq.len();
    while start < len {
        let end = (start + chunk_size).min(len);
        chunks.push(&seq[start..end]);
        start = end;
    }
    chunks
}

fn build_gfa_for_single_genome(species_taxid: &str, genome_path: &str, wd: &str, chopped: usize) -> Option<String> {
    let out_path = format!("{}/{}.gfa", wd, species_taxid);
    if Path::new(&out_path).exists() {
        return Some(out_path);
    }

    let genome_records = match read_fasta(genome_path) {
        Ok(records) => records,
        Err(err) => {
            eprintln!("Failed to read {}: {}", genome_path, err);
            return None;
        }
    };

    let genome_file_name = Path::new(genome_path).file_name()?.to_str()?;
    let mut parts = genome_file_name.split('_');
    let genome_name = format!("{}_{}", parts.next()?, parts.next()?);

    let mut s_lines = Vec::new();
    let mut l_lines = Vec::new();
    let mut w_lines = Vec::new();

    let mut node_id = 1;

    for (seqid, sequence) in genome_records {
        let seqid = seqid.split_whitespace().next().unwrap();
        let new_seqid = format!("{}#1#{}", genome_name, seqid);
        let chunks = chunk_sequence(&sequence, chopped);
        let seq_len = sequence.len();

        // Build S lines
        for chunk in &chunks {
            s_lines.push(format!("S\t{}\t{}\n", node_id, chunk));
            node_id += 1;
        }

        // Build L lines
        let first_id = node_id - chunks.len();
        for i in 0..(chunks.len() - 1) {
            l_lines.push(format!(
                "L\t{}\t+\t{}\t+\t0M\n",
                first_id + i,
                first_id + i + 1
            ));
        }

        // Build W line
        let walk = (first_id..node_id)
            .map(|i| format!(">{}", i))
            .collect::<Vec<_>>()
            .join("");
        let tokens: Vec<&str> = new_seqid.split('#').collect();
        w_lines.push(format!(
            "W\t{}\t{}\t{}\t0\t{}\t{}\n",
            tokens[0], tokens[1], tokens[2], seq_len, walk
        ));
    }

    if s_lines.len() != node_id - 1 {
        eprintln!("S lines mismatch: expected {}, got {}", node_id - 1, s_lines.len());
        return None;
    }

    if w_lines.len() != 1 {
        eprintln!("W lines error for taxid {}", species_taxid);
        return None;
    }

    let mut writer = BufWriter::new(File::create(&out_path).expect("Cannot write GFA file"));
    writeln!(writer, "H\tVN:Z:1.1").unwrap();
    for s in s_lines {
        writer.write_all(s.as_bytes()).unwrap();
    }
    for l in l_lines {
        writer.write_all(l.as_bytes()).unwrap();
    }
    for w in w_lines {
        writer.write_all(w.as_bytes()).unwrap();
    }

    Some(out_path)
}

pub fn parallel_build_gfa(single_species2genome: &FxHashMap<String, PathBuf>, wd: &PathBuf, chunk_size: usize) -> Vec<String> {
    log::info!("Pangenome building second step: building pangenome for species with one genomes...");
    if !wd.is_dir() {
        create_dir_all(wd).expect("Failed to create working directory");
    }
    let wd = wd.to_string_lossy();
    let gfa_paths = single_species2genome
        .par_iter()
        .filter_map(|(taxid, genome_path)| {
            let path_str = match genome_path.to_str() {
                Some(s) => s,
                None => {
                    eprintln!("Invalid path for taxid: {}", taxid);
                    return None;
                }
            };
            
            build_gfa_for_single_genome(taxid, path_str, &wd, chunk_size)
        })
        .collect();
    gfa_paths
}

pub fn convert(gfa_path: &PathBuf, vg: &PathBuf) -> Result<()> {
    let vg_path = gfa_path.with_extension("vg");
    if vg_path.exists() { return Ok(());}
    let threads = 1;
    let cmd = format!(
        "{} convert -g {} -p -t {} > {}",
        vg.display(),
        gfa_path.display(),
        threads,
        vg_path.display(),
    );
    
    run_shell_command(&cmd)?;
    if !vg_path.exists() {
        return Err(anyhow!("VG file not created: {:?}", vg_path));
    }
    Ok(())
}

pub fn single_zip(gfa_path: &PathBuf, output_dir: &PathBuf, range_file: &PathBuf, compress_type: &CompressType) -> Result<()> {
    let threads = 1;
    zip(&gfa_path, output_dir, range_file, threads, &compress_type)?;
    remove_file(gfa_path)?;
    Ok(())
}

pub fn parallel_convert_and_zip(gfa_paths: &Vec<String>, args: &Cli) -> Result<()> {
    let gloal_config = get_config();
    if args.is_save() {
        let output_dir = &gloal_config.zip_wd;
        let range_file = &gloal_config.range_file;
        let compress_type = CompressType::from_args(&args)?;
        gfa_paths
            .par_iter()
            .try_for_each(|gfa_path_s| -> anyhow::Result<()> {
                let gfa_path = PathBuf::from(gfa_path_s);
                convert(&gfa_path, &args.vg)?;
                single_zip(&gfa_path, output_dir, range_file, &compress_type)?;
                Ok(())
            })?;
    } else {
        gfa_paths
            .par_iter()
            .try_for_each(|gfa_path_s| -> anyhow::Result<()> {
                let gfa_path = PathBuf::from(gfa_path_s);
                convert(&gfa_path, &args.vg)?;
                Ok(())
            })?;        
    }
    log::info!("Pangenome building second step completed successfully.");
    Ok(())
}
