
use crate::cli::Cli;
use crate::types::*;
use crate::constants::*;
use crate::utils::output_genomes_info;
use sylph::api::contain_new;
use std::{collections::HashSet, path::Path};
use std::fs::{create_dir_all, File, remove_dir_all};
use anyhow::Result;
// use rayon::prelude::*;
use needletail::parse_fastx_file;
use fxhash::FxHashMap;

pub fn construct(genomes_metadata: &mut Vec<GenomesInfo>, data_type: &DataType, check_points: &CheckPoints, args: &Cli) -> Result<()> {
    if check_points.fast_query {
        log::info!("Sylph query and filter with ANI threshold {}...", args.ani);
        create_dir_all(args.tmp_dir.join("sylph_db"))?;
        if !args.tmp_dir.join(SYLPH_QUERY_DONE).exists() {
            let reference_genomes_paths: Vec<String> = genomes_metadata.iter().map(|x| x.path_id.clone()).collect();
            let read_files: Vec<String> = args.read_files.as_ref().unwrap().iter().map(|x| x.to_string_lossy().into_owned()).collect();
            let is_interleaved = match data_type {
                DataType::ShortReadPairedInter => true,
                _ => false
            };
            // println!("read_files: {:?}", read_files);
            // println!("is_interleaved: {}", is_interleaved);
            let out_prefix = Some(args.tmp_dir.join("sylph_db/query_result").to_string_lossy().into_owned());
            let mut query_res = contain_new(&reference_genomes_paths, &read_files, &out_prefix, is_interleaved);
            query_res.retain(|x| x.2 >= args.ani);
            // for (_, s, _) in query_res.iter_mut() {
            //     if let Some(first) = s.split_whitespace().next() {
            //         *s = first.to_string();
            //     }
            // }
            // let seq2genome: FxHashMap<String, String> = genomes_metadata
            //     .par_iter()
            //     .flat_map(|genomes_info| process_genome(&genomes_info.genome_id, &genomes_info.path_id))
            //     .collect();
            let retained_genomes_path: Vec<String> = query_res.iter().map(|x| x.0.clone()).collect::<HashSet<String>>().into_iter().collect();
            genomes_metadata.retain(|x| retained_genomes_path.contains(&x.path_id));
            output_genomes_info(&genomes_metadata, &args.tmp_dir.join("sylph_db/filter_genomes_info.txt"))?;
            output_genomes_info(&genomes_metadata, &args.db.join("genomes_info.txt"))?;
            File::create(args.tmp_dir.join(SYLPH_QUERY_DONE))?;
            log::info!("Sylph query and filter with {}...done", args.ani);
        } else {
            log::info!("Sylph query and filter with {} has been done. Skipping.", args.ani);
        }
        if check_points.fast_query_and_filter && !args.debug {
            remove_dir_all(&args.tmp_dir)?;
        }
    }

    Ok(())
}

fn _process_genome(genome: &String, file_path: &String) -> FxHashMap<String, String> {
    let mut seq_ids = FxHashMap::default();

    let reader = parse_fastx_file(Path::new(file_path));
    if let Err(e) = reader {
        eprintln!("Failed to open file {}: {}", file_path, e);
        return seq_ids;
    }

    let mut reader = reader.unwrap();
    while let Some(record) = reader.next() {
        match record {
            Ok(rec) => {
                let seq_id = String::from_utf8_lossy(rec.id()).split_whitespace().next().unwrap().to_string();
                seq_ids.insert(seq_id, genome.to_string());
            }
            Err(e) => eprintln!("Error reading record in {}: {}", file_path, e),
        }
    }

    seq_ids
}