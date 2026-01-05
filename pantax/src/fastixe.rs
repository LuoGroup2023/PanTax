
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufWriter, Write};
use regex::Regex;
use rayon::prelude::*;
use needletail::parse_fastx_file;
use crossbeam_channel::unbounded;
use rust_htslib::bgzf::{Writer as BGZFWriter, CompressionLevel};
use rust_htslib::tpool::ThreadPool;
use anyhow::Result;

pub fn process_all_fasta_and_merge(files: &[PathBuf], output_file_path: &Path) -> Result<()> {
    let reg = String::from("[^_]+_[^_]+");
    let uppercase = true;
    let (sender, receiver) = unbounded();
    files.par_iter().for_each_with(sender.clone(), |s, file_path| {
        if let Ok(results) = process_fasta_needle(file_path.as_ref(), &reg, uppercase) {
            for line in results {
                s.send(line).expect("Failed to send!");
            }
        }

    });
    drop(sender);

    // bgzip output
    let merge_bgzip_output = true;
    // compress default
    let compression_level = None;
    let threads = 1;
    let mut writer = create_all_fasta_and_merge_writer(output_file_path, merge_bgzip_output, compression_level, threads)?;

    for line in receiver {
        writer.write_all(line.as_bytes())?;
    }

    Ok(())
} 

fn process_fasta_needle(file_path: &Path, regex: &str, uppercase: bool) -> Result<Vec<String>> {
    let mut results = vec![];
    
    let prefix = extract_prefix_from_path(file_path, regex)?;

    if let Ok(mut reader) = parse_fastx_file(file_path) {
        while let Some(record) = reader.next() {
            if let Ok(seqrec) = record {
                if let Some(first_record_id) = std::str::from_utf8(seqrec.id()).unwrap().split_whitespace().next() {
                    if let Ok(seq) = std::str::from_utf8(seqrec.raw_seq()) {
                        let seq_formatted = if uppercase {
                            seq.to_ascii_uppercase()
                        } else {
                            seq.to_string()
                        };
                        results.push(format!(">{}{}\n{}\n", prefix, first_record_id, seq_formatted));
                    } else {
                        eprintln!("Invalid UTF-8 sequence in file: {:?}", file_path);
                    }
                } else {
                    eprintln!("Missing recored id in file: {:?}", file_path);
                }
                
            }
        }
    }
    Ok(results)
}

fn extract_prefix_from_path(file_path: &Path, regex: &str) -> Result<String, std::io::Error> {
    let input_file_name = file_path.file_stem()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidInput, "Invalid file path"))?
        .to_string_lossy();

    let re = Regex::new(regex)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Invalid regex: {}", e)))?;

    match re.captures(&input_file_name) {
        Some(caps) => {
            let matched = caps.get(0).map(|m| m.as_str()).unwrap_or("");
            Ok(format!("{}#0#", matched))
        }
        // None => Err(std::io::Error::new(
        //     std::io::ErrorKind::InvalidInput,
        //     format!(
        //         "Filename '{}' doesn't match regex '{}'",
        //         input_file_name, regex
        //     ),
        // )),
        None => {
            Ok(format!("{}#0#", input_file_name))
        }
    }
}

fn create_all_fasta_and_merge_writer(output_file_path: &Path, bgzip_output: bool, compression_level: Option<u32>, threads: usize) -> Result<Box<dyn Write>> {
    if bgzip_output {
        let compression = compression_level
            .map(|l| CompressionLevel::Level(l as i8))
            .unwrap_or(CompressionLevel::Default);
        let mut writer = BGZFWriter::from_path_with_level(output_file_path, compression)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        let tpool = ThreadPool::new(threads as u32)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        writer
            .set_thread_pool(&tpool)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        Ok(Box::new(writer))

    } else {
        let output_file = File::create(output_file_path)?;
        Ok(Box::new(BufWriter::new(output_file)))
    }
}
