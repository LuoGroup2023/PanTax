
use crate::cli::Cli;
use crate::task_scheduling::TaskConfig;
use crate::types::Graph;
use std::path::{Path, PathBuf};
use std::process::exit;
use bincode;
use lz4_flex::frame::{FrameEncoder, FrameDecoder};
use zstd::stream::{Encoder, Decoder};
use std::io::{BufRead, BufReader};
use std::io::{BufWriter, Write};
use std::collections::BTreeMap;
use regex::Regex;
use std::fs::{OpenOptions, File};
use fs2::FileExt;
use std::time::Duration;
use std::thread::sleep;
use anyhow::Result;
#[cfg(feature = "h5")]
use hdf5_metno::{File as H5File, Result};


#[derive(Debug, Hash, PartialEq, Eq)]
pub enum CompressType {
    Lz4,
    Zstd,
    Serialized,
    #[cfg(feature = "h5")]
    Hdf5,
}

impl CompressType {
    pub fn from_config(args: &TaskConfig) -> Result<CompressType> {
        if args.lz4 {
            Ok(CompressType::Lz4)
        } else if args.zstd {
            Ok(CompressType::Zstd)
        } else if args.save {
            Ok(CompressType::Serialized)
        } else {
            eprintln!("Error: no compression format specified!");
            exit(1)
        }
    }

    pub fn from_args(args: &Cli) -> Result<CompressType> {
        if args.lz {
            Ok(CompressType::Lz4)
        } else if args.zstd {
            Ok(CompressType::Zstd)
        } else if args.save {
            Ok(CompressType::Serialized)
        } else {
            eprintln!("Error: no compression format specified!");
            exit(1)
        }
    }

    fn _from_file(file: &PathBuf) -> Result<CompressType> {
        match file.extension().and_then(|s| s.to_str()) {
            Some("lz4") => Ok(CompressType::Lz4),
            Some("zst") => Ok(CompressType::Zstd),
            Some("bin") => Ok(CompressType::Serialized),
            Some(_) => {
                eprintln!("Error: it does not support this format!");
                exit(1)
            }
            None => {
                eprintln!("Error: file '{}' has no extension; cannot determine compression type. It will be considered uncompressed.", file.display());
                exit(1)
            }
        }
    }
}



fn read_and_zip_gfa(gfa_file: &Path, previous: usize, output_dir: &Path, compressed: &CompressType, compress_lvl: i32, threads: usize) -> Result<(String, usize, usize, usize)> {
    let file = File::open(gfa_file)?;
    let reader = BufReader::new(file);

    let mut nodes_len = Vec::new();
    let mut paths: BTreeMap<String, Vec<usize>> = BTreeMap::new();

    // extract nodes from paths
    let re_w = Regex::new(r"-?\d+").unwrap();  
    let re_p = Regex::new(r"\d+").unwrap(); 

    let mut node_index = 0usize;
    let mut min_values: Vec<usize> = Vec::new();
    let mut max_values: Vec<usize> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('S') {
            // parse line starts with 'S', that is nodes
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() < 3 {
                continue;
            }
            let node_id: usize = parts[1].parse()?;
            let node_id_adj = node_id - 1 - previous;
            assert_eq!(node_index, node_id_adj, "Node ID out of order or mismatch");
            node_index += 1;

            let seq = parts[2];
            let node_len = seq.len() as i64;
            assert!(node_len > 0, "Error: Node length 0 appears in the GFA!");
            nodes_len.push(node_len);
        } else if line.starts_with('W') || line.starts_with('P') {
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.is_empty() {
                continue;
            }

            let reverse_flag ;
            let haplotype_id: String;
            let mut path: Vec<usize>;

            if parts[0] == "W" {
                // W line
                haplotype_id = parts[1].to_string();
                let last_field = parts.last().unwrap_or(&"");
                reverse_flag = last_field.starts_with('<');

                path = re_w.find_iter(last_field)
                    .map(|mat| {
                        let v: isize = mat.as_str().parse().unwrap();
                        (v as isize - 1 - previous as isize) as usize
                    })
                    .collect();
                // println!("path: {:?}", path);
            } else {
                // P line
                haplotype_id = parts[1].split('#').next().unwrap_or("").to_string();
                let path_field = parts.get(2).unwrap_or(&"");
                reverse_flag = path_field.split(',').next().unwrap_or("").ends_with('-');

                path = re_p.find_iter(path_field)
                    .map(|mat| {
                        let v: usize = mat.as_str().parse().unwrap();
                        v - 1 - previous
                    })
                    .collect();
            }

            if reverse_flag {
                path.reverse();
            }

            min_values.push(*path.iter().min().unwrap());
            max_values.push(*path.iter().max().unwrap());

            // merge same haplotype_id
            // # Multiple chromosomes from the same genome merge into one path. 
            // # Although this would introduce two non-existent trio nodes for each additional chromosome, 
            // # it is not worth mentioning that only two nodes are relative to all nodes in the entire graph.
            paths.entry(haplotype_id).and_modify(|p| p.extend(&path)).or_insert(path);
        }
    }

    let min_value = *min_values.iter().min().unwrap();
    let max_value = *max_values.iter().max().unwrap();

    let is_pan = if paths.len() > 1 {
        1
    } else {
        0
    };

    let graph = Graph { nodes_len, paths };
    // println!("graph: {:?}", graph);

    let stem = gfa_file.file_stem().unwrap().to_str().unwrap();
    let mut zip_file = PathBuf::from(output_dir);
    match compressed {
        CompressType::Serialized => {
            zip_file.push(format!("{stem}.bin"));
            if !zip_file.exists() {
                let mut serialized_gfa_file = BufWriter::new(
                    File::create(&zip_file)
                        .expect(&format!("{:?} path not valid; exiting ", zip_file)),
                );        
                // bincode::encode_into_std_write(&graph, &mut serialized_gfa_file, bincode::config::standard()).unwrap();
                bincode::serialize_into(&mut serialized_gfa_file, &graph).unwrap();

            } else {
                log::debug!("- {:?} exists. Skipping.", zip_file)
            }
        },
        CompressType::Lz4 => {
            // let serialized_stream = bincode::encode_to_vec(&graph, bincode::config::standard()).unwrap();
            // let serialized_zip_gfa_file = BufWriter::new(
            //     File::create("ser.bin.lz4")
            //         .expect(&format!("{} path not valid; exiting ", "ser.bin.lz4")),
            // ); 
            // let mut compressor = frame::FrameEncoder::new(serialized_zip_gfa_file);
            // compressor.write_all(&serialized_stream)?;
            // compressor.finish()?;
            zip_file.push(format!("{stem}.bin.lz4"));
            if !zip_file.exists() {
                let file = File::create(zip_file)?;
                let mut compressed_writer = FrameEncoder::new(BufWriter::new(file));
                // bincode::encode_into_std_write(&graph, &mut compressed_writer, bincode::config::standard())?;
                bincode::serialize_into(&mut compressed_writer, &graph).unwrap();
                compressed_writer.finish()?;
            } else {
                log::debug!("- {:?} exists. Skipping.", zip_file)
            }

        },
        // there is pure rust crate ruzstd, but it does not yet reach the speed, ratio or configurability of the original zstd library.
        CompressType::Zstd => {
            zip_file.push(format!("{stem}.bin.zst"));
            if !zip_file.exists() {
                let file = File::create(zip_file)?;
                let mut compressed_writer = Encoder::new(file, compress_lvl)?;
                compressed_writer.multithread(threads as u32)?;
                // bincode::encode_into_std_write(&graph, &mut compressed_writer, bincode::config::standard())?;
                bincode::serialize_into(&mut compressed_writer, &graph).unwrap();
                compressed_writer.finish()?;
            } else {
                log::debug!("- {:?} exists. Skipping.", zip_file)
            }
        },
        #[cfg(feature = "h5")]
        CompressType::Hdf5 => {
            eprintln!("Don't support hdf5 zip!");
            exit(1)
        }
    }  

    Ok((stem.into(), min_value + 1, max_value + 1, is_pan))
}

pub fn load_from_zip_graph(
    file_path: &Path,
    compress_type: CompressType,
) -> Result<Graph, Box<dyn std::error::Error>> {
    match compress_type {
        CompressType::Serialized => {
            let file = File::open(file_path)?;
            let mut reader = BufReader::new(file);
            // let graph: Graph = bincode::decode_from_std_read(&mut reader, bincode::config::standard())?;
            let graph: Graph = bincode::deserialize_from(&mut reader)?;
            Ok(graph)
        }
        CompressType::Lz4 => {
            let file = File::open(file_path)?;
            let decoder = FrameDecoder::new(BufReader::new(file));
            let mut reader = BufReader::new(decoder);
            // let graph: Graph = bincode::decode_from_std_read(&mut reader, bincode::config::standard())?;
            let graph: Graph = bincode::deserialize_from(&mut reader)?;
            Ok(graph)
        }
        CompressType::Zstd => {
            let file = File::open(file_path)?;
            let decoder = Decoder::new(BufReader::new(file))?;
            let mut reader = BufReader::new(decoder);
            // let graph: Graph = bincode::decode_from_std_read(&mut reader, bincode::config::standard())?;
            let graph: Graph = bincode::deserialize_from(&mut reader)?;
            Ok(graph)
        }
        #[cfg(feature = "h5")]
        CompressType::Hdf5 => {
            let file = H5File::open(file_path)?;
            let paths_group = file.group("paths")?;
            let mut paths = BTreeMap::new();
            for obj in paths_group.member_names()? {
                let ds = paths_group.dataset(&obj)?;
                let data: Vec<usize> = ds.read_1d::<usize>()?.to_vec();
                paths.insert(obj, data);
            }
        
            let node_len_group = file.group("node_len")?;
            let node_len_ds = node_len_group.dataset("node_len")?;
            let nodes_len: Vec<i64> = node_len_ds.read_1d::<i64>()?.to_vec();
        
            Ok(Graph { nodes_len, paths })            
        }
    }

}


fn write_locked_line<P: AsRef<Path>>(
    path: P,
    taxid: &str,
    min_value: usize,
    max_value: usize,
    is_pan: usize,
) -> std::io::Result<()> {
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .read(true) 
        .open(&path)?;

    loop {
        match file.try_lock_exclusive() {
            Ok(_) => break,
            Err(_) => sleep(Duration::from_millis(100)),
        }
    }

    {
        let mut writer = BufWriter::new(&file);
        writeln!(writer, "{}\t{}\t{}\t{}", taxid, min_value, max_value, is_pan)?;
        writer.flush()?;
    }

    fs2::FileExt::unlock(&file)?;
    Ok(())
}

pub fn zip(gfa_file: &PathBuf, output_dir: &PathBuf, output_range_file: &PathBuf, threads: usize, compress_type: &CompressType) -> Result<()> {
    let zstd_threads = if threads > 8 {
        8
    } else {
        threads
    };    
    let compress_lvl = 3;
    let (taxid, min_value, max_value, is_pan) = read_and_zip_gfa(&gfa_file, 0, &output_dir, &compress_type, compress_lvl, zstd_threads)?;

    write_locked_line(output_range_file, &taxid, min_value, max_value, is_pan)?;
    Ok(())
}

