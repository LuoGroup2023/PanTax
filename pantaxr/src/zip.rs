
use crate::cmdline::ZipArgs;
use crate::profile::Graph;
use crate::clog::init_logger;
use std::path::{Path, PathBuf};
use std::fs::{OpenOptions, File, create_dir, read_to_string};
use std::io::{BufRead, BufReader};
use std::io::{BufWriter, Write};
use std::process::exit;
use std::collections::BTreeMap;
use regex::Regex;
use bincode;
use lz4_flex::frame::{FrameEncoder, FrameDecoder};
use zstd::stream::{Encoder, Decoder};
use rayon::prelude::*;

use log::*;
use fs2::FileExt;
use std::time::Duration;
use std::thread::sleep;

#[derive(Debug, Hash, PartialEq, Eq)]
pub enum CompressType {
    Lz4,
    Zstd,
    Serialized,
}

impl CompressType {
    fn from_args(args: &ZipArgs) -> Result<Self, &'static str> {
        if args.lz {
            Ok(CompressType::Lz4)
        } else if args.zstd {
            Ok(CompressType::Zstd)
        } else if args.serialized {
            Ok(CompressType::Serialized)
        } else {
            eprintln!("Error: no compression format specified!");
            exit(1)
        }
    }

    fn from_file(file: &PathBuf) -> Result<Self, &'static str> {
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

// use std::io::{Read, Write};
// use zstd::stream::read::Decoder;
// use zstd::stream::write::Encoder;
// use zstd::zstd_safe::{CCtx, DCtx};
// use dashmap::DashMap;
// use std::cell::RefCell;
// use zstd::zstd_safe::CParameter;

// thread_local! {
//     static CTX_MAP: RefCell<DashMap<CompressType, CCtx<'static>>> = RefCell::new(DashMap::new());
//     static DCTX_MAP: RefCell<DashMap<CompressType, DCtx<'static>>> = RefCell::new(DashMap::new());
// }

// pub struct DefaultCompressor {}

// impl DefaultCompressor {
//     pub fn new() -> Self {
//         DefaultCompressor {}
//     }

//     pub async fn compress_raw(
//         &self,
//         data: &[u8],
//         compress_type: CompressType,
//     ) -> Result<Vec<u8>, Error> {
//         match compress_type {
//             CompressType::Zstd => {
//                 let ret = CTX_MAP.with(|map_cell| {
//                     let map = map_cell.borrow();
//                     let mut ctx_entry = map.entry(compress_type).or_default();
//                     let writer = Vec::new();
//                     let mut o = Encoder::with_context(writer, ctx_entry.value_mut());
//                     o.write_all(data)?;
//                     o.finish()
//                 });
//                 Ok(ret?)
//             }
//             _ => Ok(data.to_vec()),
//         }
//     }

//     pub async fn decompress_raw(
//         &self,
//         data: &[u8],
//         compress_type: CompressType,
//     ) -> Result<Vec<u8>, Error> {
//         match compress_type {
//             CompressType::Zstd => DCTX_MAP.with(|map_cell| {
//                 let map = map_cell.borrow();
//                 let mut ctx_entry = map.entry(compress_type).or_default();
//                 let mut decoder = Decoder::with_context(data, ctx_entry.value_mut());
//                 let mut output = Vec::new();
//                 decoder.read_to_end(&mut output)?;
//                 Ok(output)
//             }),
//             _ => Ok(data.to_vec()),
//         }
//     }
// }



fn check_args_valid(args: &ZipArgs) -> std::io::Result<()> {

    let level = if args.trace {
        "trace"
    } else if args.debug {
        "debug"
    } else {
        "info"
    };

    init_logger(level);

    if let Some(input_gfa_file) = &args.input_gfa {
        if !input_gfa_file.is_file() {
            eprintln!("Error: {:?} does not exists", input_gfa_file);
            exit(1)
        }
    } else if let Some(input_gfa_list) = &args.input_gfa_list {
        if !input_gfa_list.is_file() {
            eprintln!("Error: {:?} does not exists", input_gfa_list);
            exit(1)
        }        
    } else {
        eprintln!("Error: please provide gfa file");
        exit(1)
    }

    if !args.output_dir.exists() {
        create_dir(&args.output_dir)?;
    }

    Ok(())
}

fn read_and_zip_gfa(gfa_file: &Path, previous: usize, output_dir: &Path, compressed: &CompressType, compress_lvl: i32, threads: usize) -> Result<(String, usize, usize, usize), Box<dyn std::error::Error>> {
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
                debug!("- {:?} exists. Skipping.", zip_file)
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
                debug!("- {:?} exists. Skipping.", zip_file)
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
                debug!("- {:?} exists. Skipping.", zip_file)
            }

        },

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

pub fn zip(args: ZipArgs) -> Result<(), Box<dyn std::error::Error>> {
    check_args_valid(&args)?;
    let zstd_threads = if args.threads > 8 {
        8
    } else {
        args.threads
    };

    if args.decompress {
        let compress_type = CompressType::from_file(&args.input_gfa.as_ref().unwrap())?;
        let _graph = load_from_zip_graph(args.input_gfa.as_ref().unwrap(), compress_type)?;
        // println!("graph: {:?}", _graph);
    } else {
        let compress_type = CompressType::from_args(&args)?;
        if let Some(input_gfa) = &args.input_gfa {
            let (taxid, min_value, max_value, is_pan) = read_and_zip_gfa(input_gfa, 0, &args.output_dir, &compress_type, args.compress_lvl, zstd_threads)?;
            // println!("taxid: {}, min: {}, max: {}, is_pan: {}", taxid, min_value, max_value, is_pan);
            // println!("{}\t{}\t{}\t{}", taxid, min_value, max_value, is_pan);
            write_locked_line(args.output_range_file, &taxid, min_value, max_value, is_pan)?;
        } else if let Some(input_gfa_list) = &args.input_gfa_list {
            rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
            let paths: Vec<PathBuf> = read_to_string(input_gfa_list).unwrap().lines().map(|s| s.into()).collect();
            let mut results: Vec<(usize, String, usize, usize, usize)> = paths
                .par_iter()
                .enumerate()
                .map(|(i, path)| {
                    let (taxid, min_value, max_value, is_pan) = read_and_zip_gfa(
                        path,
                        0,
                        &args.output_dir,
                        &compress_type,
                        args.compress_lvl,
                        zstd_threads,
                    ).expect("read_and_zip_gfa failed");
                    (i, taxid, min_value, max_value, is_pan)
                })
                .collect();
            
            results.sort_by_key(|x| x.0);
            
            let mut adjusted_results: Vec<(String, usize, usize, usize)> = Vec::with_capacity(results.len());
            
            let mut offset = 0;
            for (_i, taxid, mut start, mut end, is_pan) in results {
                start += offset;
                end += offset;
                adjusted_results.push((taxid, start, end, is_pan));
                offset = end;
            }
            
            let mut writer = BufWriter::new(File::create(args.output_range_file)?);
            for (taxid, start, end, is_pan) in adjusted_results.iter() {
                writeln!(writer, "{}\t{}\t{}\t{}", taxid, start, end, is_pan)?;
            }
        }
        
                
    }

    Ok(())
}