
use crate::cmdline::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

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

fn parallel_build_gfa(info_path: &PathBuf, wd: &PathBuf, chunk_size: usize) -> Vec<String> {
    let file = File::open(info_path).expect("Cannot open species_eq1_genomes_info");
    let reader = BufReader::new(file);
    let tasks: Vec<(String, String)> = reader
        .lines()
        .filter_map(|line| {
            line.ok().and_then(|l| {
                let toks: Vec<&str> = l.trim().split('\t').collect();
                if toks.len() == 2 {
                    Some((toks[0].to_string(), toks[1].to_string()))
                } else {
                    None
                }
            })
        })
        .collect();
    
    if !wd.is_dir() {
        create_dir_all(wd).expect("Failed to create working directory");
    }
    let wd = wd.to_string_lossy();
    tasks
        .into_par_iter()
        .filter_map(|(taxid, genome_path)| {
            build_gfa_for_single_genome(&taxid, &genome_path, &wd, chunk_size)
        })
        .collect()
}




pub fn build(args: BuildArgs) -> Result<(), Box<dyn std::error::Error>> {

    let info_path = &args.species_eq1_genomes_info;
    let wd = &args.wd;
    let chunk_size = 1024;

    let gfa_paths = parallel_build_gfa(info_path, wd, chunk_size);

    let done_path = wd.join("gfa_build_done.txt");
    let mut done_writer = BufWriter::new(File::create(&done_path).expect("Cannot write done file"));
    for path in gfa_paths {
        writeln!(done_writer, "{}", path).unwrap();
    }

    Ok(())
}