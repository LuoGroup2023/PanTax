
use std::path::Path;
use std::fs::{read_dir, File, create_dir_all};
use std::io::{BufReader, BufRead, BufWriter, Write};
use rayon::prelude::*;
use crossbeam_channel::unbounded;
use needletail::parse_fastx_file;
use regex::Regex;
use crate::cmdline::*;
use log::*;
use rust_htslib::bgzf::{Writer as BGZFWriter, CompressionLevel};
use rust_htslib::tpool::ThreadPool;
use rust_htslib::faidx::build;


fn check_args_valid(args: &FastixeArgs) {
    let level: LevelFilter;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
    } else {
        level = log::LevelFilter::Info
    }

    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    if args.input_list.is_none()
        && args.input_directory.is_none()
        && args.input_files.is_none()
    {
        error!("No genome found! Exiting.");
        std::process::exit(1);    
    }

    // if args.prefix.is_none() {
    //     warn!("No prefix provided; use default regex.");
    // }
}

fn parse_line_file(path: &Path, vec: &mut Vec<String>) -> std::io::Result<()> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let path: &Path = line.as_ref();
        if path.exists() {
            vec.push(line);
        } else {
            eprintln!("{:?} does not exist!", path);
        }
    }
    Ok(())
}

fn is_fasta<P: AsRef<Path>>(file: P) -> bool {
    let path = file.as_ref();

    if let Some(file_name) = path.to_str() {
        file_name.ends_with(".fa") ||
            file_name.ends_with(".fna") ||
            file_name.ends_with(".fasta") ||
            file_name.ends_with(".fa.gz") ||
            file_name.ends_with(".fna.gz") ||
            file_name.ends_with(".fasta.gz")
    } else {
        eprintln!("{:?} can not convert to utf-8 string", path);
        false
    }

}

fn parse_files(args: &FastixeArgs, input_genomes: &mut Vec<String>) {
    let mut all_files = vec![];

    if let Some(ref input_files) = args.input_files {
        for input_file in input_files {
            if input_file.is_file() && is_fasta(&input_file) {
                all_files.push(input_file.to_string_lossy().to_string());
            }
        }
    }

    if let Some(ref input_list) = args.input_list {
        if input_list.exists() {
            parse_line_file(&input_list, &mut all_files).unwrap();
        }
    }

    if let Some(ref input_directory) = args.input_directory {
        if input_directory.exists() && input_directory.is_dir() {
            for entry in read_dir(&input_directory).unwrap() {
                let entry = entry.unwrap();
                let path = entry.path();
                if path.is_file() && is_fasta(&path) {
                    all_files.push(path.to_string_lossy().to_string());
                }
            }
        }
    }

    input_genomes.extend(all_files);

}


fn extract_prefix_from_path(file_path: &Path, regex: &str) -> Result<String, std::io::Error> {
    let input_file_name = file_path.file_name()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidInput, "Invalid file path"))?
        .to_string_lossy();

    let re = Regex::new(regex)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Invalid regex: {}", e)))?;

    match re.captures(&input_file_name) {
        Some(caps) => {
            let matched = caps.get(0).map(|m| m.as_str()).unwrap_or("");
            Ok(format!("{}#0#", matched))
        }
        None => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "Filename '{}' doesn't match regex '{}'",
                input_file_name, regex
            ),
        )),
    }
}

fn process_fasta_needle(file_path: &Path, regex: &str, uppercase: bool) -> Result<Vec<String>, std::io::Error> {
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


fn create_all_fasta_and_merge_writer(output_file_path: &Path, bgzip_output: bool, compression_level: Option<u32>, threads: usize) -> std::io::Result<Box<dyn Write>> {
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

fn process_all_fasta_and_merge(args: &FastixeArgs, files: &[String], output_file_path: &Path) -> std::io::Result<()> {
    let (sender, receiver) = unbounded();
    files.par_iter().for_each_with(sender.clone(), |s, file_path| {
        if let Ok(results) = process_fasta_needle(file_path.as_ref(), &args.reg, args.uppercase) {
            for line in results {
                s.send(line).expect("Failed to send!");
            }
        }

    });
    drop(sender);

    let mut writer = create_all_fasta_and_merge_writer(output_file_path, args.merge_bgzip_output, args.compression_level, args.threads)?;

    for line in receiver {
        writer.write_all(line.as_bytes())?;
    }

    Ok(())
} 

pub fn fastixe(args: FastixeArgs) -> Result<(), Box<dyn std::error::Error>> {
    let mut input_genomes = vec![];

    check_args_valid(&args);
    parse_files(&args, &mut input_genomes);
    create_dir_all(&args.out_directory)?;
    let mut merged_path = Path::new(&args.out_directory).join(&args.merge_output_file_path);
    if args.merge_bgzip_output {
        merged_path.set_extension("gz");
    };
    process_all_fasta_and_merge(&args,&input_genomes, &merged_path)?;
    build(&merged_path).expect("Failed to build FASTA index"); 

    Ok(())
}
