
use log::*;
use chrono::Local;
use std::fs::{File, metadata};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use regex::Regex;
use rayon::prelude::*;
use polars::frame::DataFrame;
use polars::prelude::*;
use nalgebra::DMatrix;
use crate::{cmdline::ProfileArgs, rcls};
use rcls::save_output_to_file;

fn check_args_valid(args: &ProfileArgs) {
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

    for (name, path) in &[
        ("GAF mapping file", &args.input_aln_file),
        ("Species range file", &args.range_file),
        ("Species length file", &args.species_len_file),
    ] {
        if !path.exists() {
            panic!("Error: {} not found at path {:?}", name, path);
        }
        if !path.is_file() {
            panic!("Error: {} path {:?} is not a file", name, path);
        }

        let metadata = metadata(path).expect("Failed to get metadata");
        if metadata.len() == 0 {
            panic!("Error: {} at {:?} is empty", name, path);
        }
    }

}

fn equal_length_read_cls(df: &DataFrame, read_len: i64, isfilter: bool) -> PolarsResult<DataFrame> {

    let lazy_df = df.clone().lazy();
    let grouped_df = if !isfilter {
        lazy_df
            .group_by_stable([col("species")])
            .agg([(len()*lit(read_len)).alias("base_count")])
            .collect()?
    } else {

        // all reads count，without mapq filtered
        let read_count_df = lazy_df.clone()
            .group_by([col("species")])
            .agg([len().alias("read_count")]);
        
        // count less_multi (reads with mapq 3~60) and uniq_count (reads with mapq 60)
        let filtered = lazy_df
            .filter(col("mapq").gt_eq(lit(3)).and(col("mapq").lt_eq(lit(60))));
        
        let filtered_agg_df = filtered
            .group_by([col("species")])
            .agg([
                len().alias("less_multi"),
                col("mapq").eq(lit(60)).sum().alias("uniq_count"),
            ]);
        
        // join read_count + less_multi + uniq_count
        let joined = filtered_agg_df
            .join(read_count_df, [col("species")], [col("species")], JoinType::Inner.into());
        
        // filter
        let result = joined
            .filter(
                col("uniq_count").gt(lit(0))
                    .and(col("less_multi").gt(
                        col("read_count").cast(DataType::Float64) / lit(10.0),
                    ))
            )
            .select([col("species"), (col("read_count") * lit(read_len)).alias("base_count")])
            .collect()?;
        result
    };

    Ok(grouped_df)
}

fn non_equal_length_read_cls(df: &DataFrame, isfilter: bool) -> PolarsResult<DataFrame> {
    let lazy_df = df.clone().lazy();
    let grouped_df = if !isfilter {
        lazy_df
            .group_by_stable([col("species")])
            .agg([(col("read_len").sum()).alias("base_count")])
            .collect()?
    } else {

        // all reads count，without mapq filtered
        let read_count_df = lazy_df.clone()
            .group_by([col("species")])
            .agg([len().alias("read_count"), (col("read_len").sum()).alias("base_count")]);
        
        // count less_multi (reads with mapq 3~60) and uniq_count (reads with mapq 60)
        let filtered = lazy_df
            .filter(col("mapq").gt_eq(lit(3)).and(col("mapq").lt_eq(lit(60))));
        
        let filtered_agg_df = filtered
            .group_by([col("species")])
            .agg([
                len().alias("less_multi"),
                col("mapq").eq(lit(60)).sum().alias("uniq_count"),
            ]);
        
        // join read_count + less_multi + uniq_count
        let joined = filtered_agg_df
            .join(read_count_df, [col("species")], [col("species")], JoinType::Inner.into());
        
        // filter
        let result = joined
            .filter(
                col("uniq_count").gt(lit(0))
                    .and(col("less_multi").gt(
                        col("read_count").cast(DataType::Float64) / lit(10.0),
                    ))
            )
            .select([col("species"), col("base_count")])
            .collect()?;
        result
    };

    Ok(grouped_df)
}

fn species_profile(args: &ProfileArgs, rcls_df: &DataFrame) -> PolarsResult<DataFrame> {
    // read species length file
    let schema = Schema::from_iter(vec![
        Field::new("species".into(), DataType::String),
        Field::new("len".into(), DataType::Float64),
    ]);
    let species_len_df = LazyCsvReader::new(&args.species_len_file)
        .with_has_header(false)
        .with_separator(b'\t')
        .with_schema(Some(Arc::new(schema)))
        .finish()?;
    
    // Determine if the read length is of equal length
    let read_len_unique = rcls_df
        .clone()
        .lazy()
        .select([col("read_len")])
        .limit(1000)
        .unique(None, UniqueKeepStrategy::First)
        .collect()?;
    let unique_len = read_len_unique.height();
    let grouped_df: LazyFrame;
    if unique_len == 1 {
        let read_len= read_len_unique.column("read_len")?.i64()?.get(0).unwrap();
        grouped_df = equal_length_read_cls(&rcls_df, read_len, args.filtered)?.lazy();
    } else {
        grouped_df = non_equal_length_read_cls(&rcls_df, args.filtered)?.lazy();
    }

    
    // let species_grouped_df = grouped_df
    //     .join(&species_len_df, ["species"], ["species"], JoinArgs::new(JoinType::Left), None)?;

    let species_grouped_df = grouped_df
        .join(species_len_df, [col("species")], [col("species")], JoinType::Left.into())
        .select([
            col("species"),
            (col("base_count") / col("len")).alias("absolute_abund")
        ]);
    let species_profile = species_grouped_df
        .select([
            col("species").alias("species_taxid"),
            (col("absolute_abund") / col("absolute_abund").sum()).alias("predicted_abundance"),
            col("absolute_abund").alias("predicted_coverage")
        ])
        .sort(["predicted_abundance"], SortMultipleOptions::new().with_order_descending(true))
        .collect()?;
    let mut write_species_profile = species_profile.clone();
    save_output_to_file(&mut write_species_profile, &args.output_name.to_string_lossy(), true)?;
    Ok(species_profile)
}

pub struct Graph {
    pub nodes_len: Vec<usize>,
    pub paths: HashMap<String, Vec<usize>>,
}

pub fn read_gfa(gfa_file: &Path, previous: usize) -> Result<Graph, Box<dyn std::error::Error>> {
    let file = File::open(gfa_file)?;
    let reader = BufReader::new(file);

    let mut nodes_len = Vec::new();
    let mut paths: HashMap<String, Vec<usize>> = HashMap::new();

    // extract nodes from paths
    let re_w = Regex::new(r"-?\d+").unwrap();  
    let re_p = Regex::new(r"\d+").unwrap(); 

    let mut node_index = 0usize;

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
            let node_len = seq.len();
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

            // merge same haplotype_id
            // # Multiple chromosomes from the same genome merge into one path. 
            // # Although this would introduce two non-existent trio nodes for each additional chromosome, 
            // # it is not worth mentioning that only two nodes are relative to all nodes in the entire graph.
            paths.entry(haplotype_id).and_modify(|p| p.extend(&path)).or_insert(path);
        }
    }

    Ok(Graph { nodes_len, paths })
}

pub struct SpeciesRange {
    pub species: String,
    pub start: u32,
    pub end: u32,
}

fn load_species_range(args: &ProfileArgs, species_profile_df: &DataFrame) -> PolarsResult<Vec<SpeciesRange>> {
    let schema = Schema::from_iter(vec![
        Field::new("species".into(), DataType::String),
        Field::new("start".into(), DataType::Int64),
        Field::new("end".into(), DataType::Int64),
        Field::new("is_pan".into(), DataType::Int32),
    ]);
    let species_range_df = LazyCsvReader::new(&args.range_file.clone())
        .with_has_header(false)
        .with_separator(b'\t')
        .with_schema(Some(Arc::new(schema)))
        .finish()?;

    let species_profile_lazy_df = species_profile_df.clone().lazy().rename(&["species_taxid"], &["species"], false);
    
    let species_profile_to_range = species_profile_lazy_df
        .join(species_range_df, [col("species")], [col("species")], JoinArgs::new(JoinType::Left));

    let filtered_species_range = if let Some(species_str) = &args.designated_species {
        let species_vec: Vec<&str> = species_str
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
        let designated_species_series = Series::new("designated_species_series".into(), species_vec);
        species_profile_to_range
            .filter(col("species").is_in(lit(designated_species_series)))
            .select([
                col("species"),
                col("start"),
                col("end")
            ])
            .collect()?
    } else {
        species_profile_to_range
            .select([
                col("species"),
                col("start"),
                col("end")
            ])
            .collect()?
    };    



    let species = filtered_species_range.column("species")?.str()?.into_no_null_iter();
    let starts = filtered_species_range.column("start")?.i64()?.into_no_null_iter();
    let ends = filtered_species_range.column("end")?.i64()?.into_no_null_iter();
    Ok(species
        .zip(starts)
        .zip(ends)
        .map(|((s, start), end)| SpeciesRange {
            species: s.to_string(),
            start: start as u32,
            end: end as u32,
        })
        .collect())
    
}

fn trio_nodes_info(graph: &Graph) -> (HashMap<(usize, usize, usize), usize>, Vec<usize>, DMatrix<u8>) {
    let mut trio_nodes_set = HashSet::new();
    let mut hap_trio_paths = HashMap::new();
    let haps: Vec<String> = graph.paths.keys().cloned().collect();

    // construct hap → trio_path
    for (hap, path) in &graph.paths {
        let trio_path: Vec<(usize, usize, usize)> = path.windows(3)
            .map(|w| (w[0], w[1], w[2]))
            .collect();
        trio_nodes_set.extend(trio_path.iter().cloned());
        hap_trio_paths.insert(hap.clone(), trio_path);
    }

    let trio_nodes: Vec<(usize, usize, usize)> = trio_nodes_set.into_iter().collect();
    let trio_index_map: HashMap<_, _> = trio_nodes.iter().enumerate().map(|(i, t)| (*t, i)).collect();
    let mut presence_matrix = DMatrix::zeros(trio_nodes.len(), graph.paths.len());

    let mut count_per_trio = vec![0; trio_nodes.len()];
    for (hap_idx, hap) in haps.iter().enumerate() {
        if let Some(trios) = hap_trio_paths.get(hap) {
            for t in trios {
                if let Some(&idx) = trio_index_map.get(t) {
                    presence_matrix[(idx, hap_idx)] = 1;
                    count_per_trio[idx] += 1;
                }
            }
        }
    }

    // save unique trio nodes
    let mut unique_trio_nodes = HashMap::new();
    let mut unique_lengths = Vec::new();
    let mut rows = Vec::new();
    for (i, count) in count_per_trio.iter().enumerate() {
        if *count == 1 {
            let trio = trio_nodes[i];
            unique_trio_nodes.insert(trio, unique_trio_nodes.len());
            let trio_len = graph.nodes_len[trio.0] + graph.nodes_len[trio.1] + graph.nodes_len[trio.2];
            unique_lengths.push(trio_len);
            rows.push(presence_matrix.row(i).transpose());
        }
    }

    let final_matrix = if !rows.is_empty() {
        DMatrix::from_columns(&rows)
    } else {
        DMatrix::zeros(0, graph.paths.len())
    };

    (unique_trio_nodes, unique_lengths, final_matrix)
}

#[derive(Debug)]
pub struct OptimizeRes {
    pub otu: String,
    pub absolute_abund: f64,
}

fn optimize_otu(args: &ProfileArgs, otu: &str, start: u32, end: u32) -> Option<Vec<OptimizeRes>> {
    log::debug!("Reading {} graph information, start: {}, end: {}", otu, start, end);
    let species_gfa_dir = args.db.join("species_gfa");
    let graph = if species_gfa_dir.exists() {
        let gfa_file = species_gfa_dir.join(format!("{}.gfa", otu));
        read_gfa(&gfa_file, 0).map_err(|e| format!("GFA read error: {}", e)).unwrap()
    } else {
        // panic!("{:?} does not exist!", species_gfa_dir);
        eprintln!("{:?} does not exist! Skipping.", species_gfa_dir);
        return None;
    };

    let (unique_trio_nodes, unique_trio_nodes_len, hap2unique_trio_nodes_m) = trio_nodes_info(&graph);


    Some(vec![
        OptimizeRes { otu: otu.to_string(), absolute_abund: 10.0 }
    ])
}

pub fn profile(args: ProfileArgs) -> Result<(), Box<dyn std::error::Error>> {
    check_args_valid(&args);

    log::info!("{} - Read classification...", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let rcls_df = rcls::rcls_profile(&args)?;
    let filter_unmapped_rcls_df = rcls_df
        .clone()
        .lazy()
        .filter(col("species").neq(lit("U")))
        .collect()?;

    log::info!("{} - Species level profiling...", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let species_profile_df = species_profile(&args, &filter_unmapped_rcls_df)?;

    log::info!("{} - Strain level profiling...", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let species_ranges = load_species_range(&args, &species_profile_df)?;
    let total = species_ranges.len();
    let counter = Arc::new(Mutex::new(0));
    
    let results: Vec<_> = species_ranges
        .par_iter()
        .filter_map(|range| {
            let otu: &String = &range.species;
            let result: Option<Vec<OptimizeRes>> = optimize_otu(&args, otu, range.start, range.end);
    
            let mut count = counter.lock().unwrap();
            *count += 1;
            if *count % 10 == 0 {
                let percent = (*count as f64) * 100.0 / total as f64;
                println!("percentage: {:.1}%", percent);
            }
    
            result
        })
        .flatten()
        .collect();
    println!("{:?}", results);
    Ok(())
}