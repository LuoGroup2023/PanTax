
use crate::{cmdline::ProfileArgs, rcls};
use rcls::{save_output_to_file, load_gaf_file_lazy};
use crate::zip::{load_from_zip_graph, CompressType};
use crate::clog::init_logger;
use log::*;
use regex::Regex;
use rayon::prelude::*;

use polars::frame::DataFrame;
use polars::prelude::*;
use polars::functions::concat_df_horizontal;

use std::fs::{File, create_dir_all};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;
use std::io::{BufRead, BufReader};

use cpu_time::ProcessTime;
use dashmap::DashMap;

use grb::prelude::*;
use grb::expr::LinExpr;

use nalgebra::{DMatrix, RowDVector};

use rand::{SeedableRng, rngs::StdRng, seq::IndexedRandom};

// use bincode::{Encode, Decode};
use serde::{Serialize, Deserialize};
// use std::io::{BufWriter, Write}; 

#[derive(Debug, Default, Clone)]
pub struct InputFile {
    // db: Option<PathBuf>,
    pub gaf_file: Option<PathBuf>,
    pub range_file: Option<PathBuf>,
    len_file: Option<PathBuf>,
    reads_binning_file: Option<PathBuf>,
    species_abund_file: Option<PathBuf>,
    strain_abund_file: Option<PathBuf>,
    genomes_metadata_file: Option<PathBuf>,
}

fn check_args_valid(args: &ProfileArgs, input_file: &mut InputFile) {

    if !args.species && !args.strain {
        panic!("Please choose profiling level with --species or/and --strain.");
    }

    let level = if args.trace {
        "trace"
    } else if args.debug {
        "debug"
    } else {
        "info"
    };

    init_logger(level);

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    let db = &args.db;
    if !(db.exists() && db.is_dir()) {
        panic!("Error: Specified PanTax database directory {:?} is not a valid directory path", db);
    }

    let wd = &args.wd;
    if !(wd.exists() && wd.is_dir()) {
        panic!("Error: Specified PanTax work directory {:?} is not a valid directory path", wd);
    }

    fn is_existing_file(name: &str, opt_path: &Option<PathBuf>) {
        let is_exist = opt_path
            .as_ref()
            .map_or(false, |p| p.exists() && p.is_file());
        if !is_exist {
            let path_str = opt_path
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "<none>".to_string());
            panic!("Error: Specified {} '{}' is not a valid file path", name, path_str);
        }
    }

    fn choose_existing_file_from_two_files(
        name: &str,
        path1: &Option<PathBuf>,
        path2: &Option<PathBuf>,
    ) -> PathBuf {
        let exist1 = path1.as_ref().filter(|p| p.exists() && p.is_file());
        let exist2 = path2.as_ref().filter(|p| p.exists() && p.is_file());
    
        match (exist1, exist2) {
            (Some(p1), Some(_)) => p1.clone(), // if two file both exist, return first file
            (Some(p1), None) => p1.clone(),
            (None, Some(p2)) => p2.clone(),
            (None, None) => {
                let path_str1 = path1
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| "<none>".to_string());
                let path_str2 = path2
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| "<none>".to_string());
                panic!(
                    "Error: Neither {} '{}' nor '{}' is a valid file path",
                    name, path_str1, path_str2
                );
            }
        }
    }

    let species_abund_file_is_exist = if !args.force {
        let species_abund_file = args.wd.join("species_abundance.txt");
        let exists = species_abund_file.exists() && species_abund_file.is_file();
        if exists {
            input_file.species_abund_file = Some(species_abund_file);
        }
        exists
    } else {
        false
    };

    let strain_abund_file_is_exist = if !args.force {
        let strain_abund_file = args.wd.join("strain_abundance.txt");
        let exists = strain_abund_file.exists() && strain_abund_file.is_file();
        if exists {
            input_file.strain_abund_file = Some(strain_abund_file);
        }
        exists
    } else {
        false
    };

    if args.species && !species_abund_file_is_exist {
        is_existing_file("GAF mapping file", &args.input_aln_file);
        input_file.gaf_file = args.input_aln_file.clone();

        let range_file = choose_existing_file_from_two_files("species range file", &args.range_file, &Some(args.db.join("species_range.txt")));
        input_file.range_file = Some(range_file);

        let len_file = choose_existing_file_from_two_files("species length file", &args.species_len_file, &Some(args.db.join("species_genomes_stats.txt")));
        input_file.len_file = Some(len_file);

        if args.strain {
            let genomes_metadata_file = choose_existing_file_from_two_files("genomes metadata file", &args.genomes_metadata, &Some(args.db.join("genomes_info.txt")));
            input_file.genomes_metadata_file = Some(genomes_metadata_file);    
        }
    } else if args.strain && !strain_abund_file_is_exist {
        is_existing_file("GAF mapping file", &args.input_aln_file);
        input_file.gaf_file = args.input_aln_file.clone();

        let range_file = choose_existing_file_from_two_files("species range file", &args.range_file, &Some(args.db.join("species_range.txt")));
        input_file.range_file = Some(range_file);

        let default_reads_binning_file = args.wd.join("reads_classification.tsv");
        // is_existing_file("reads binning file", &Some(reads_binning_file.clone())); 
        let reads_binning_file = choose_existing_file_from_two_files("reads binning file", &args.reads_binning_file, &Some(default_reads_binning_file.clone()));
        input_file.reads_binning_file = Some(reads_binning_file);

        let species_abund_file = args.wd.join("species_abundance.txt");
        is_existing_file("species abundance file", &Some(species_abund_file.clone()));
        input_file.species_abund_file = Some(species_abund_file);

        let genomes_metadata_file = choose_existing_file_from_two_files("genomes metadata file", &args.genomes_metadata, &Some(args.db.join("genomes_info.txt")));
        input_file.genomes_metadata_file = Some(genomes_metadata_file);


    }

    if !args.output_dir.exists() {
        create_dir_all(args.output_dir.clone())
            .unwrap_or_else(|e| panic!("Failed to create directory {:?}: {}", args.output_dir.clone(), e));        
    }

}

fn file_exists(path_opt: &Option<PathBuf>) -> bool {
    match path_opt {
        Some(path) => path.exists(),
        None => false,
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

fn species_profiling(args: &ProfileArgs, input_file: &InputFile, rcls_df: &DataFrame) -> PolarsResult<DataFrame> {
    // read species length file
    let schema = Schema::from_iter(vec![
        Field::new("species".into(), DataType::String),
        Field::new("len".into(), DataType::Float64),
    ]);
    let species_len_df = LazyCsvReader::new(&input_file.len_file.as_ref().unwrap())
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
    save_output_to_file(&mut write_species_profile, &args.output_dir.join("species_abundance.txt"), true)?;
    Ok(species_profile)
}

#[derive(Debug, Clone)]
struct Record {
    read_id: String,
    path: String,
    read_path_len: i64,
    read_start: i64,
    read_end: i64,
    species: String,
}

fn dataframe_to_records_and_check_unique(df: &DataFrame) -> (Vec<Record>, bool) {
    let read_id_series = df.column("read_id").unwrap().str().unwrap();
    let path_series = df.column("path").unwrap().str().unwrap();
    let read_path_len_series = df.column("read_path_len").unwrap().i64().unwrap();
    let read_start_series = df.column("read_start").unwrap().i64().unwrap();
    let read_end_series = df.column("read_end").unwrap().i64().unwrap();
    let species_series = df.column("species").unwrap().str().unwrap();

    let mut seen = HashSet::new();
    let mut unique = true;

    let mut records = Vec::with_capacity(df.height());

    for i in 0..df.height() {
        let read_id = read_id_series.get(i).unwrap();
        if !seen.insert(read_id) {
            unique = false;
        }
        
    if let (Some(path), Some(read_path_len), Some(read_start), Some(read_end), Some(species)) = (
        path_series.get(i),
        read_path_len_series.get(i),
        read_start_series.get(i),
        read_end_series.get(i),
        species_series.get(i),
    ) {
        records.push(Record {
            read_id: read_id.to_string(),
            path: path.to_string(),
            read_path_len,
            read_start,
            read_end,
            species: species.to_string(),
        });
    } else {
        // S0R4970856/2    150     0       150     +       <104888926<104888924<104888923<104888921        82      5       *
        // some reads like this, very few 
        debug!("row {} -> path: {:?}, read_path_len: {:?}, read_start: {:?}, read_end: {:?}, species: {:?}", i, path_series.get(i), read_path_len_series.get(i), read_start_series.get(i), read_end_series.get(i), species_series.get(i));
    }

    }

    (records, unique)
}

fn process_with_duplicates(records: Vec<Record>) -> BTreeMap<String, Vec<Record>> {
    let mut grouped: HashMap<String, Vec<Record>> = HashMap::new();
    for r in records {
        grouped.entry(r.read_id.clone()).or_default().push(r);
    }

    let mut species_map: BTreeMap<String, Vec<Record>> = BTreeMap::new();

    for (_read_id, group) in grouped.into_iter() {
        let species_set: HashSet<_> = group.iter().map(|r| &r.species).collect();
        if species_set.len() == 1 {
            let species = group[0].species.clone();
            let entry = species_map.entry(species.clone()).or_default();
            for (i, r) in group.into_iter().enumerate() {
                let mut id = r.read_id.clone();
                if i > 0 {
                    id = format!("{}_{}", id, i + 1);
                }
                entry.push(Record {
                    read_id: id,
                    path: r.path.clone(),
                    read_path_len: r.read_path_len,
                    read_start: r.read_start,
                    read_end: r.read_end,
                    species: species.clone(),
                });
            }
        }
    }

    species_map
}

fn group_reads_by_species(df: &DataFrame) -> BTreeMap<String, Vec<Record>> {
    let (records, is_read_id_unique) = dataframe_to_records_and_check_unique(&df);
    if is_read_id_unique {
        // rayon
        let map: DashMap<String, Vec<Record>> = DashMap::new();
        records.par_iter()
            .for_each(|r| {
                map.entry(r.species.clone())
                    .or_default()
                    .push(r.clone());
            });

        // DashMap convert to BTreeMap
        let mut final_map: BTreeMap<String, Vec<Record>> = BTreeMap::new();
        for (key, value) in map.into_iter() {
            final_map.insert(key, value);
        }

        final_map
    } else {
        // fallback
        process_with_duplicates(records)
    }
}



// #[derive(Debug, Encode, Decode)]
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph {
    pub nodes_len: Vec<i64>,
    pub paths: BTreeMap<String, Vec<usize>>,
}

fn read_gfa(gfa_file: &Path, previous: usize) -> Result<Graph, Box<dyn std::error::Error>> {
    let file = File::open(gfa_file)?;
    let reader = BufReader::new(file);

    let mut nodes_len = Vec::new();
    let mut paths: BTreeMap<String, Vec<usize>> = BTreeMap::new();

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
            let node_len = seq.len() as i64;
            assert!(node_len > 0, "Error: Node length 0 appears in the GFA!");
            nodes_len.push(node_len);
        } else if line.starts_with('W') || line.starts_with('P') {
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.is_empty() {
                continue;
            }

            // let reverse_flag ;
            let haplotype_id: String;
            let path: Vec<usize>;

            if parts[0] == "W" {
                // W line
                haplotype_id = parts[1].to_string();
                let last_field = parts.last().unwrap_or(&"");
                // reverse_flag = last_field.starts_with('<');

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
                // reverse_flag = path_field.split(',').next().unwrap_or("").ends_with('-');

                path = re_p.find_iter(path_field)
                    .map(|mat| {
                        let v: usize = mat.as_str().parse().unwrap();
                        v - 1 - previous
                    })
                    .collect();
            }

            // if reverse_flag {
            //     path.reverse();
            // }

            // merge same haplotype_id
            // # Multiple chromosomes from the same genome merge into one path. 
            // # Although this would introduce two non-existent trio nodes for each additional chromosome, 
            // # it is not worth mentioning that only two nodes are relative to all nodes in the entire graph.
            paths.entry(haplotype_id).and_modify(|p| p.extend(&path)).or_insert(path);
        }
    }

    Ok(Graph { nodes_len, paths })
}

struct SpeciesRange {
    species: String,
    start: u32,
    end: u32,
}

fn load_species_range(args: &ProfileArgs, input_file: &InputFile, species_profile_df: &DataFrame) -> PolarsResult<Vec<SpeciesRange>> {
    let schema = Schema::from_iter(vec![
        Field::new("species".into(), DataType::String),
        Field::new("start".into(), DataType::Int64),
        Field::new("end".into(), DataType::Int64),
        Field::new("is_pan".into(), DataType::Int32),
    ]);
    let species_range_df = LazyCsvReader::new(&input_file.range_file.as_ref().unwrap())
        .with_has_header(false)
        .with_separator(b'\t')
        .with_schema(Some(Arc::new(schema)))
        .finish()?;

    let mode_filtered_species_range_df = match args.mode {
        Some(0) => {
            species_range_df.filter(col("is_pan").eq(lit(0)))
        }
        Some(1) => {
            species_range_df.filter(col("is_pan").eq(lit(1)))
        }
        _ => {
            species_range_df
        }
    };

    let ds_filtered_species_range_df = if let Some(species_str) = &args.designated_species {
        if !species_str.trim().is_empty() && species_str != "None" {
            let species_vec: Vec<&str> = species_str
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .collect();  
            let designated_species_series = Series::new("designated_species_series".into(), species_vec);
            mode_filtered_species_range_df.filter(col("species").is_in(lit(designated_species_series)))
        } else {
            mode_filtered_species_range_df
        }
    } else {
        mode_filtered_species_range_df
    };

    let ds_filtered_species_range_df_check = ds_filtered_species_range_df.clone().collect()?;
    if ds_filtered_species_range_df_check.height() == 0 {
        warn!("The filtering before strain profiling has removed all species, and the program is about to exit without performing strain profiling.");
        std::process::exit(0);
    }

    let species_profile_lazy_df = species_profile_df.clone().lazy().rename(&["species_taxid"], &["species"], false);
    
    let abund_filtered_species_profile_lazy_df = species_profile_lazy_df.filter(col("predicted_abundance").gt(lit(args.min_species_abundance)));

    let species_profile_to_range = abund_filtered_species_profile_lazy_df
        .join(ds_filtered_species_range_df, [col("species")], [col("species")], JoinArgs::new(JoinType::Inner));

    let final_species_to_range = species_profile_to_range.select([
        col("species"),
        col("start"),
        col("end")
    ]).collect()?;

    // let species_profile_to_range = species_profile_lazy_df
    //     .join(species_range_df, [col("species")], [col("species")], JoinArgs::new(JoinType::Inner));

    // let filtered_species_range = if let Some(species_str) = &args.designated_species {
    //     let species_vec: Vec<&str> = species_str
    //         .split(',')
    //         .map(|s| s.trim())
    //         .filter(|s| !s.is_empty())
    //         .collect();
    //     let designated_species_series = Series::new("designated_species_series".into(), species_vec);
    //     species_profile_to_range
    //         .filter(col("species").is_in(lit(designated_species_series)))
    //         .select([
    //             col("species"),
    //             col("start"),
    //             col("end")
    //         ])
    //         .collect()?
    // } else {
    //     species_profile_to_range
    //         .select([
    //             col("species"),
    //             col("start"),
    //             col("end")
    //         ])
    //         .collect()?
    // };    



    let species = final_species_to_range.column("species")?.str()?.into_no_null_iter();
    let starts = final_species_to_range.column("start")?.i64()?.into_no_null_iter();
    let ends = final_species_to_range.column("end")?.i64()?.into_no_null_iter();
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

fn trio_nodes_info(graph: &Graph) -> (HashMap<(usize, usize, usize), usize>, Vec<i64>, DMatrix<u8>) {
    let mut trio_nodes_set = HashSet::new();
    let mut hap_trio_paths = HashMap::new();
    let haps: Vec<String> = graph.paths.keys().cloned().collect();
    
    // println!("haps: {:?}", haps);

    // construct hap → trio_path
    for (hap, path) in &graph.paths {
        // if path.windows(3).any(|w| w == [27, 29, 30]) {
        //     println!("Found (27,29,30) in path {:?}", hap);
        // }
        let trio_path: Vec<(usize, usize, usize)> = path.windows(3)
            // .map(|w| (w[0], w[1], w[2]))
            .map(|w| {
                if w[0] > w[2] {
                    (w[2], w[1], w[0])
                } else {
                    (w[0], w[1], w[2])
                }
            })
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
                    // // for debug
                    // if *t == (293380, 293382, 293383) {
                    //     println!("count_per_trio: {}", count_per_trio[idx])
                    // }
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
            unique_lengths.push(trio_len as i64);
            rows.push(presence_matrix.row(i).into_owned());
        }
    }

    let final_matrix = if !rows.is_empty() {
        DMatrix::from_rows(&rows)
    } else {
        DMatrix::zeros(0, graph.paths.len())
    };

    // // for debug
    // use std::io::Write;
    // let file = File::create("unique_trio_nodes.txt").unwrap();
    // let mut writer = std::io::BufWriter::new(file);
    // // writeln!(writer, "i,j,k,value").unwrap();
    // // for (&(i, j, k), &v) in unique_trio_nodes.iter() {
    // //     writeln!(writer, "{},{},{},{}", i, j, k, v).unwrap();
    // // }
    // let mut sorted_entries: Vec<_> = unique_trio_nodes.iter().collect();
    // sorted_entries.sort_by_key(|(&(i, j, k), _)| (i, j, k));
    
    // for (&(i, j, k), &v) in sorted_entries {
    //     writeln!(writer, "{},{},{}\t{}", i, j, k, v).unwrap();
    // }

    (unique_trio_nodes, unique_lengths, final_matrix)
}

#[allow(unused_variables)]
fn get_node_abundances(
    otu: &String,
    nodes_len: &Vec<i64>,
    trio_nodes: &HashMap<(usize, usize, usize), usize>,
    trio_nodes_len: &Vec<i64>,
    start: usize,
    reads_cluster: &[Record],
) -> (Vec<f64>, Vec<f64>, Vec<usize>) {

    // For now we do not consider the insertions and deletions of the nodes. For the start node 
    // and end node of a read, we get matched base which is equal to the node length minus the offset.
    // For the intermediate nodes, we keep the node length as matched base. In contract, by parsing 
    // json format or gam format alignment, the intermediate nodes' matched base will consider indels.
    // In previous version, the matched base number is the node length minus the number of deletions, 
    // but don't add the number of insertions(keep node length). 
    // Actually, we can consider indels by the cs tag or cg tag in gaf file. 

    // // for debug
    // let file = File::create("nodes_len.txt").unwrap();
    // let mut writer = BufWriter::new(file);

    // for value in nodes_len {
    //     writeln!(writer, "{}", value).unwrap();
    // }

    // this regex is used to extract from W lines instead of P lines
    let re = Regex::new(r"-?\d+").unwrap();
    let trio_nodes = Arc::new(trio_nodes.clone());
    let nodes_len = Arc::new(nodes_len.clone());
    let trio_nodes_len = Arc::new(trio_nodes_len.clone());

    let bases_per_node = DashMap::new();
    let trio_nodes_bases_count = DashMap::new();
    let node_base_cov_info: DashMap<usize, (u8, usize, Vec<u8>)> = DashMap::new();
    
    for (i, &len) in nodes_len.iter().enumerate() {
        bases_per_node.insert(i, 0);
        node_base_cov_info.insert(i, (0, 0, vec![0; len as usize]));
    }

    for i in 0..trio_nodes.len() {
        trio_nodes_bases_count.insert(i, 0);
    }

    reads_cluster.par_iter().for_each(|read| {
        let read_nodes: Vec<usize> = re
            .find_iter(&read.path)
            .map(|m| m.as_str().parse::<isize>().unwrap() -1 - start as isize)
            .map(|x| x as usize)
            .collect();

        if read_nodes.is_empty() {
            return;
        }

        let start_node = read_nodes[0];
        let end_node = *read_nodes.last().unwrap();
        let mut target_len = read.read_end - read.read_start;
        let mut seen = 0;
        let mut read_nodes_len: HashMap<usize, i64> = HashMap::new();
        let mut undup_read_nodes = HashSet::new();

        // Set the initial base coverage of the node to 0    
        for &node in &read_nodes {
            read_nodes_len.insert(node, 0);
        }

        // only one node, read is shorter than the node of the graph
        if start_node == end_node && read_nodes.len() == 1 {
            // https://github.com/vgteam/vg/issues/4249
            // S0R4490169/2    150     0       150     +       <77966596       1118    1010    551 
            // the start position is beyond to end postion
            // json alignment only displays offset which is the start position
            // assert!(target_len as usize <= entry.2.len(), "read_id {} {:?} read len {} is longer than {} node len {}", read.read_id, read_nodes, target_len, start_node, entry.2.len() );
            // very few, we filter them.

            // The old Python code did not make this judgment, resulting in a smaller node abundance. 
            // Although these reads will be relatively small, sometimes the target_len will be large, especially for long reads
            if target_len < 0 {
                debug!(
                    "read: {}, read start: {}, read end: {}",
                    read.read_id, read.read_start, read.read_end
                );
                return;
            }
            *read_nodes_len.entry(start_node).or_default() += target_len;
            bases_per_node.entry(start_node).and_modify(|v| *v += target_len);

            if let Some(mut entry) = node_base_cov_info.get_mut(&start_node) {
                if read.read_start < read.read_end && read.read_end <= entry.2.len() as i64 {
                    for j in read.read_start as usize..read.read_end as usize {
                        entry.2[j] = 1;
                    }
                } else {
                    debug!(
                        "read: {} {:?} range out of bounds: [{}..{}) for entry of length {}",
                        read.read_id, read_nodes, read.read_start, read.read_end, entry.2.len()
                    );
                }
                

                entry.1 = entry.2.iter().map(|x| *x as usize).sum();
                entry.0 = if entry.1 == nodes_len[start_node] as usize { 1 } else { 0 };
            }
        } else {
            for (i, &node) in read_nodes.iter().enumerate() {
                let node_len = nodes_len[node];
                // i = 0, start node. always not all cov.
                // i = len()-1, end node. always not all cov.
                // intermediate node. all cov.
                let (node_aln_len, start_idx) = if i == 0 {
                    assert!(read.read_start <= node_len, "read start is bigger than node len: {} > {}", read.read_start, node_len);
                    // dbg!(otu, &read.read_id, i, node, node_len, read.read_start);
                    (node_len - read.read_start, read.read_start)
                } else if i == read_nodes.len() - 1 {
                    if target_len < seen { target_len = seen }
                    (target_len - seen, 0)
                } else {
                    (node_len, 0)
                };
                
                // // for debug
                // if node == 146059usize {
                //     let origin_node = node + 1 + start;
                //     println!("read id: {}, node: {}, origin_node: {}, node aln len: {}, start idx: {}", read.read_id, node, origin_node, node_aln_len, start_idx)
                // }

                if let Some(mut entry) = node_base_cov_info.get_mut(&node) {
                    for j in start_idx as usize..((start_idx + node_aln_len) as usize).min(entry.2.len()) {
                        entry.2[j] = 1;
                    }
                    entry.1 = entry.2.iter().map(|x| *x as usize).sum();
                    entry.0 = if entry.1 == nodes_len[node] as usize { 1 } else { 0 };
                }

                seen += node_aln_len;
                if undup_read_nodes.insert(node) {
                    *read_nodes_len.entry(node).or_default() += node_aln_len;
                    bases_per_node.entry(node).and_modify(|v| *v += node_aln_len);
                }
            }
        }

        if read_nodes.len() < 3 {
            return;
        }

        let read_trio_nodes: Vec<(usize, usize, usize)> = read_nodes
            .windows(3)
            .map(|w| (w[0], w[1], w[2]))
            .collect();

        for trio_node in read_trio_nodes {
            let (a, b, c) = trio_node;
            let len_sum = [a, b, c]
                .iter()
                .map(|n| read_nodes_len.get(n).copied().unwrap_or(0))
                .sum::<i64>();

            if let Some(&i) = trio_nodes
                .get(&trio_node)
                .or_else(|| trio_nodes.get(&(trio_node.2, trio_node.1, trio_node.0)))
            {
                trio_nodes_bases_count.entry(i).and_modify(|v| *v += len_sum);
            } 
            // // for debug
            // else {
            //     use std::fs::OpenOptions;
            //     let mut file = OpenOptions::new()
            //         .create(true)
            //         .append(true)
            //         .open("missing_trio_nodes.txt")
            //         .expect("Failed to open missing_trio_nodes.txt");
            //     writeln!(file, "{}\t{:?}", read.read_id, trio_node).expect("Failed to write trio_node");
            // }
        }
    });

    // // for debug
    // if otu.clone() == "Myxococcus_xanthus".to_string() || otu == "173" {
    //     use std::io::Write;
    //     let file = File::create(format!("{otu}_trio_nodes_bases_count.txt")).unwrap();
    //     let mut writer = std::io::BufWriter::new(file);
    //     for entry in trio_nodes_bases_count.iter() {
    //         if *entry.value() > 0 { writeln!(writer, "{}\t{}", entry.key(), entry.value()).unwrap() }
    //     }
    // }

    // // for debug
    // if otu == "173" || otu == "9" {
    //     use std::io::Write;
    //     let file = File::create(format!("{otu}_trio_nodes_bases_count.txt")).unwrap();
    //     let mut writer = std::io::BufWriter::new(file);
    //     let mut data: Vec<(usize, i64)> = trio_nodes_bases_count
    //         .par_iter()
    //         .filter_map(|e| {
    //             let k = *e.key();
    //             let v = *e.value();
    //             (v > 0).then_some((k, v))
    //         })
    //         .collect();
    
    //     data.par_sort_unstable_by_key(|(k, _)| *k);
    
    //     let mut buffer = String::with_capacity(data.len() * 12);
    //     for (i, val) in data {
    //         buffer.push_str(&format!("{i}\t{val}\n"));
    //     }
    
    //     writer.write_all(buffer.as_bytes()).unwrap();
    // }

    // // for debug
    // if otu == "173" || otu == "9" {
    //     use std::io::Write;
    //     let file = File::create(format!("{otu}_base_aln.txt")).unwrap();
    //     let mut writer = std::io::BufWriter::new(file);
    //     let mut data: Vec<(usize, i64)> = bases_per_node
    //         .par_iter()
    //         .filter_map(|e| {
    //             let k = *e.key();
    //             let v = *e.value();
    //             (v > 0).then_some((k, v))
    //         })
    //         .collect();
    
    //     data.par_sort_unstable_by_key(|(k, _)| *k);
    
    //     let mut buffer = String::with_capacity(data.len() * 12);
    //     for (i, val) in data {
    //         buffer.push_str(&format!("{i}\t{val}\n"));
    //     }
    
    //     writer.write_all(buffer.as_bytes()).unwrap();
    // }

    // rayon calculate abundance
    let node_abundance_vec: Vec<f64> = nodes_len
        .par_iter()
        .enumerate()
        .map(|(i, &len)| {
            if len == 0 {
                panic!("Node length 0 for node {}", i);
            }
            let base_cov = bases_per_node.get(&i).map(|v| *v).unwrap_or(0);
            base_cov as f64 / len as f64
        })
        .collect();

    // // for debug
    // if otu == "173" || otu == "9" {
    //     use std::io::Write;
    //     let file = File::create(format!("{otu}_node_depth.txt")).unwrap();
    //     let mut writer = std::io::BufWriter::new(file);
    
    //     let mut buffer = String::with_capacity(node_abundance_vec.len());
    //     for (i, val) in node_abundance_vec.iter().enumerate() {
    //         if *val > 0.0 { buffer.push_str(&format!("{i}\t{val}\n")) }
    //     }
    
    //     writer.write_all(buffer.as_bytes()).unwrap();
    // }

    let trio_node_abundance_vec: Vec<f64> = trio_nodes_len
        .par_iter()
        .enumerate()
        .map(|(i, &len)| {
            if len == 0 {
                panic!("Trio node length 0 for node {}", i);
            }
            let trio_nodes_cov = trio_nodes_bases_count.get(&i).map(|v| *v).unwrap_or(0);
            trio_nodes_cov as f64 / len as f64
        })
        .collect();

    let mut node_base_cov = vec![0; nodes_len.len()];
    for i in 0..nodes_len.len() {
        if let Some(v) = node_base_cov_info.get(&i) {
            node_base_cov[i] = v.1;
        }
    }

    (node_abundance_vec, trio_node_abundance_vec, node_base_cov)
}

fn zscore_filter(data: &[f64], threshold: f64) -> Vec<f64> {
    if data.is_empty() {
        return vec![];
    }

    if data.iter().any(|&x| x.is_nan()) {
        panic!("Input data contains NaN values.");
    }

    let mean: f64 = data.iter().sum::<f64>() / data.len() as f64;
    let std: f64 = (data.iter()
        .map(|&x| (x - mean).powi(2))
        .sum::<f64>() / data.len() as f64)
        .sqrt();

    if std == 0.0 {
        return vec![];
    }

    data.iter()
        .cloned()
        .filter(|&x| ((x - mean) / std).abs() < threshold)
        .collect()
}

#[derive(Debug, Clone)]
struct GurobiOptVar {
    otu: String,
    hap_metrics: Vec<HapMetrics>,
    possible_paths_idx: Vec<usize>,
    second_possible_paths_idx: Vec<usize>,
    orign_n_haps: usize,
    hap2trio_nodes_m_size: usize,
    same_path_flag: bool,
    second_opt: bool,
}

#[derive(Debug, Default, Clone)]
struct HapMetrics {
    otu: Option<String>,
    hap_id: Option<String>,
    unique_trio_nodes_fraction: Option<f64>,
    frequencies_mean: Option<f64>,
    path_cov_ratio: Option<f64>,
    first_sol: Option<f64>,
    divergence: Option<f64>,
    second_sol: Option<f64>,
    is_rescue: Option<bool>,
    total_cov_diff: Option<f64>,
    // genome_ID: Option<String>,
}

fn first_filter_paths(
    gurobi_opt_var: &mut GurobiOptVar,
    paths: &BTreeMap<String, Vec<usize>>,
    hap2trio_nodes_m: &DMatrix<u8>,
    trio_node_abundances: &Vec<f64>,
    node_abundance_vec: &Vec<f64>,
    args: &ProfileArgs,
) {
    for (i, hap_id) in paths.keys().enumerate() {
        gurobi_opt_var.hap_metrics[i].otu = Some(gurobi_opt_var.otu.clone());
        gurobi_opt_var.hap_metrics[i].hap_id = Some(hap_id.clone());
    }

    let orign_n_haps = paths.len();
    let hap2trio_nodes_m_size = hap2trio_nodes_m.len();
    gurobi_opt_var.orign_n_haps = orign_n_haps;
    gurobi_opt_var.hap2trio_nodes_m_size = hap2trio_nodes_m_size;
    // debug!("{}\thap2trio_nodes_m_size:{}", gurobi_opt_var.otu, hap2trio_nodes_m_size);
    if orign_n_haps != 1 && hap2trio_nodes_m_size != 0 {
        // // for debug
        // let file = File::create("hap2trio_nodes_m.csv").unwrap();
        // let mut writer = BufWriter::new(file);
        // for i in 0..hap2trio_nodes_m.nrows() {
        //     for j in 0..hap2trio_nodes_m.ncols() {
        //         write!(writer, "{}", hap2trio_nodes_m[(i, j)]).unwrap();
        //         if j < hap2trio_nodes_m.ncols() - 1 {
        //             write!(writer, ",").unwrap();
        //         }
        //     }
        //     writeln!(writer).unwrap();
        // }    

        for (hap_idx, (hap_id, _)) in paths.iter().enumerate() {
            // Retrieve the unique trio node index that this path passes through in hap2trio_nodes_m
            let trio_idxs: Vec<usize> = (0..hap2trio_nodes_m.nrows())
                .filter(|&i| hap2trio_nodes_m[(i, hap_idx)] > 0)
                .collect();
            // debug!("{}\ttrio_idxs:{}", gurobi_opt_var.otu, trio_idxs.len());
            let trio_idxs_len = trio_idxs.len();
            if trio_idxs_len == 0 { continue; }

            // Obtain the abundance of the unique trio node for this path
            let trio_idx_set: std::collections::HashSet<usize> = trio_idxs.iter().cloned().collect();
            let abundances: Vec<f64> = trio_node_abundances
                .iter()
                .enumerate()
                .map(|(i, &v)| if trio_idx_set.contains(&i) { v } else { 0.0 })
                .collect();

            let non_zero_frequencies: Vec<f64> = abundances
                .iter()
                .cloned() 
                .filter(|&x| x > 0.0)
                .collect(); 

            let unique_trio_nodes_fraction = non_zero_frequencies.len() as f64 / trio_idxs_len as f64;
            let unique_trio_nodes_fraction_rounded: f64 = (unique_trio_nodes_fraction * 100.0).round() / 100.0;
            // first metric
            gurobi_opt_var.hap_metrics[hap_idx].unique_trio_nodes_fraction = Some(unique_trio_nodes_fraction_rounded);

            if args.shift {
                // filter outliers
                let filter_non_zero_frequencies = zscore_filter(&non_zero_frequencies, 3.0);
                let frequencies_mean: f64 = if filter_non_zero_frequencies.is_empty() {
                    0.0
                } else {
                    filter_non_zero_frequencies.iter().sum::<f64>() / filter_non_zero_frequencies.len() as f64
                };
                let shift_unique_trio_nodes_fraction = if frequencies_mean >= 1.0 {
                    let _shift_unique_trio_nodes_fraction = args.unique_trio_nodes_fraction + (0.8-args.unique_trio_nodes_fraction) * frequencies_mean / 100.0;
                    if _shift_unique_trio_nodes_fraction > 0.8 {
                        0.8
                    } else {
                        _shift_unique_trio_nodes_fraction
                    }
                } else {
                    unique_trio_nodes_fraction * frequencies_mean
                };
                debug!("\t\t{} unique trio node abundance > 0 ratio: {}, shift unique trio nodes fraction: {}, frequencies mean: {}", hap_id, unique_trio_nodes_fraction, shift_unique_trio_nodes_fraction, frequencies_mean);
                
                if unique_trio_nodes_fraction < shift_unique_trio_nodes_fraction {
                    continue; // filter
                }
                // let frequencies_mean_rounded: f64 = (frequencies_mean * 100.0).round() / 100.0;
                // second metric
                gurobi_opt_var.hap_metrics[hap_idx].frequencies_mean = Some(frequencies_mean);
            } else {
                debug!("\t\t{} unique trio node abundance > 0 ratio: {}", hap_id, non_zero_frequencies.len());
                if unique_trio_nodes_fraction < args.unique_trio_nodes_fraction {
                    debug!("\t\t {}\t{}", gurobi_opt_var.otu, unique_trio_nodes_fraction);
                    continue; // filter
                }
                let filter_non_zero_frequencies = zscore_filter(&non_zero_frequencies, 3.0);
                let frequencies_mean: f64 = if filter_non_zero_frequencies.is_empty() {
                    0.0
                } else {
                    filter_non_zero_frequencies.iter().sum::<f64>() / filter_non_zero_frequencies.len() as f64
                };          
                // let frequencies_mean_rounded: f64 = (frequencies_mean * 100.0).round() / 100.0;
                // second metric
                gurobi_opt_var.hap_metrics[hap_idx].frequencies_mean = Some(frequencies_mean);  
            }

            gurobi_opt_var.possible_paths_idx.push(hap_idx);
        };
        debug!("\t\tFisrt filter #strains / #paths = {} / {}", gurobi_opt_var.possible_paths_idx.len(), orign_n_haps);
        // filtered_haps.push(hap_id.clone());
    } else if orign_n_haps != 1 && hap2trio_nodes_m_size == 0 {
        let mut paths_vec = paths.values();
        let first_path = paths_vec.next().unwrap();
        let all_same = paths_vec.all(|x| x == first_path);
        if all_same {
            gurobi_opt_var.same_path_flag = true;
            let non_zero_frequencies: Vec<f64> = node_abundance_vec
                .iter()
                .cloned() 
                .filter(|&x| x > 0.0)
                .collect(); 
            let frequencies_mean: f64 = if non_zero_frequencies.is_empty() {
                0.0
            } else {
                non_zero_frequencies.iter().sum::<f64>() / non_zero_frequencies.len() as f64
            }; 
            let frequencies_mean_rounded: f64 = (frequencies_mean * 100.0).round() / 100.0;
            gurobi_opt_var.hap_metrics[0].frequencies_mean = Some(frequencies_mean_rounded);   
            gurobi_opt_var.possible_paths_idx.push(0);
        } else {
            warn!("{} species path more than 1, but have trio node is None.", gurobi_opt_var.otu);
            gurobi_opt_var.possible_paths_idx = (0..orign_n_haps).collect();
        }

    } else if orign_n_haps == 1 {
        let non_zero_frequencies: Vec<f64> = node_abundance_vec
            .iter()
            .cloned() 
            .filter(|&x| x > 0.0)
            .collect(); 
        let frequencies_mean: f64 = if non_zero_frequencies.is_empty() {
            0.0
        } else {
            non_zero_frequencies.iter().sum::<f64>() / non_zero_frequencies.len() as f64
        };
        let frequencies_mean_rounded: f64 = (frequencies_mean * 100.0).round() / 100.0;
        gurobi_opt_var.hap_metrics[0].frequencies_mean = Some(frequencies_mean_rounded);   
        gurobi_opt_var.possible_paths_idx.push(0);         
    } 
    // println!("end {:?}", gurobi_opt_var.possible_paths_idx);
}

fn second_filter_paths(
    gurobi_opt_var: &mut GurobiOptVar,
    args: &ProfileArgs
) {
    let mut filter_possible_paths_idx = Vec::new();
    if gurobi_opt_var.orign_n_haps != 1 && gurobi_opt_var.hap2trio_nodes_m_size > 0 {
        gurobi_opt_var.second_opt = true;
        for possible_path_idx in gurobi_opt_var.possible_paths_idx.iter() {
            let frequencies_mean = gurobi_opt_var.hap_metrics[*possible_path_idx].frequencies_mean.unwrap_or(0.0);
            if frequencies_mean == 0.0 { continue; }
            let sol = gurobi_opt_var.hap_metrics[*possible_path_idx].first_sol.unwrap();
            let f = (sol - frequencies_mean).abs() / (sol + frequencies_mean);
            let f_rounded: f64 = (f * 100.0).round() / 100.0;
            gurobi_opt_var.hap_metrics[*possible_path_idx].divergence = Some(f_rounded);
            debug!("\t\thap_id:{}\tfrequencies_mean:{}\tfirst_sol:{}\tdivergence:{}", gurobi_opt_var.hap_metrics[*possible_path_idx].hap_id.as_ref().unwrap(), gurobi_opt_var.hap_metrics[*possible_path_idx].frequencies_mean.as_ref().unwrap(), gurobi_opt_var.hap_metrics[*possible_path_idx].first_sol.as_ref().unwrap(), f_rounded);
            
            // let epsilon = 0.0;
            // // debug 
            // if gurobi_opt_var.otu == "142" || gurobi_opt_var.otu == "195" {
            //     println!("{}\tf_rounded {} {} {} {}", gurobi_opt_var.otu, f_rounded, f > args.unique_trio_nodes_mean_count_f, f - args.unique_trio_nodes_mean_count_f, (f - args.unique_trio_nodes_mean_count_f).abs() < epsilon);
            // }
            
            if f_rounded > args.unique_trio_nodes_mean_count_f {
                if f_rounded <= 0.6 {
                    let this_strain_single_cov_ratio = gurobi_opt_var.hap_metrics[*possible_path_idx].unique_trio_nodes_fraction.unwrap() * gurobi_opt_var.hap_metrics[*possible_path_idx].path_cov_ratio.unwrap();
                    if this_strain_single_cov_ratio < args.single_cov_ratio || sol == 0.0 {
                        continue;
                    } else {
                        gurobi_opt_var.hap_metrics[*possible_path_idx].is_rescue = Some(true);
                        filter_possible_paths_idx.push(*possible_path_idx);
                    }
                } else {
                    continue;
                }
            } else if f_rounded <= args.unique_trio_nodes_mean_count_f && sol != 0.0 {
                filter_possible_paths_idx.push(*possible_path_idx);
            }

        }
        gurobi_opt_var.second_possible_paths_idx = filter_possible_paths_idx;
    } else if (gurobi_opt_var.orign_n_haps != 1 && gurobi_opt_var.hap2trio_nodes_m_size == 0 && gurobi_opt_var.same_path_flag) || (gurobi_opt_var.orign_n_haps == 1) {
        let frequencies_mean = gurobi_opt_var.hap_metrics[0].frequencies_mean.unwrap();
        if frequencies_mean > 0.0 {
            let sol = gurobi_opt_var.hap_metrics[0].first_sol.unwrap();
            let f = (sol - frequencies_mean).abs() / (sol + frequencies_mean);            
            let f_rounded: f64 = (f * 100.0).round() / 100.0;
            debug!("\t\thap_id:{}\tfrequencies_mean:{}\tfirst_sol:{}\tdivergence:{}", gurobi_opt_var.hap_metrics[0].hap_id.as_ref().unwrap(), gurobi_opt_var.hap_metrics[0].frequencies_mean.as_ref().unwrap(), gurobi_opt_var.hap_metrics[0].first_sol.as_ref().unwrap(), f_rounded);
            gurobi_opt_var.hap_metrics[0].divergence = Some(f_rounded);
            gurobi_opt_var.hap_metrics[0].second_sol = Some(sol);
        }
    } else if gurobi_opt_var.orign_n_haps != 1 && gurobi_opt_var.hap2trio_nodes_m_size == 0 && !gurobi_opt_var.same_path_flag  {
        for possible_path_idx in gurobi_opt_var.possible_paths_idx.iter() {
            gurobi_opt_var.hap_metrics[*possible_path_idx].second_sol = gurobi_opt_var.hap_metrics[*possible_path_idx].first_sol.clone();
        }
    }

}

fn sample_sorted(vec: &Vec<usize>, sample_size: usize, seed: u64) -> Vec<usize> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut sampled: Vec<usize> = vec
        .choose_multiple(&mut rng, sample_size)
        .cloned()
        .collect();
    sampled.sort_unstable();
    sampled
}

#[allow(unused_variables)]
fn gurobi_opt(
    gurobi_opt_var: &mut GurobiOptVar,
    nvert: usize,
    paths: &BTreeMap<String, Vec<usize>>,
    node_abundance_vec: &Vec<f64>,
    node_base_cov: &Vec<usize>,
    node_len: &Vec<i64>,
    args: &ProfileArgs,    
) -> Result<(), grb::Error> {
    // println!("start {:?}", gurobi_opt_var.possible_paths_idx);
    // if gurobi_opt_var.otu == "Escherichia_coli".to_string() {
    //     return Ok(());
    // }

    let mut m = Model::new("lp")?;
    let mut obj = LinExpr::new();
    let npaths = gurobi_opt_var.possible_paths_idx.len();
    
    let max_val = node_abundance_vec
        .iter()             
        .copied()  
        .fold(f64::NEG_INFINITY, f64::max);
    let mut var_map = HashMap::new();
    let mut x_vec = Vec::new();
    for i in 0..npaths {
        let var = add_var!(
            m,
            Continuous,
            name: &format!("x{}", i),
            bounds: 0.0..1.05 * max_val
        )?;
        var_map.insert(i, var);
        x_vec.push(var);
    }; 

    let mut coeff_matrix = DMatrix::<f32>::zeros(nvert, npaths);
    for (i, path) in paths.values().enumerate() {
        if let Some(pos_idx) = gurobi_opt_var.possible_paths_idx.iter().position(|&x| x == i) {
            for &v in path.iter() {
                // assert!(v < nvert, "Node index {} exceeds nvert {}", v, nvert);
                // assert!(i < npaths, "Path index {} exceeds npaths {}, {:?}", i, npaths, gurobi_opt_var.possible_paths_idx);
                coeff_matrix[(v, pos_idx)] = 1.0; 
            }
        }
    }

    let node_cov_matrix = RowDVector::<f32>::from_iterator(
        nvert,
        node_base_cov.iter().map(|&x| x as f32)
    );
    //1 × nvert  @  nvert × npaths  = 1 × npaths
    let path_cov = node_cov_matrix * &coeff_matrix;

    let node_len_matrix = RowDVector::<f32>::from_iterator(
        nvert,
        node_len.iter().map(|&x| x as f32)
    );
    //1 × nvert  @  nvert × npaths  = 1 × npaths
    let path_len = node_len_matrix * &coeff_matrix;
    let path_ratio: nalgebra::Matrix<f32, nalgebra::Const<1>, nalgebra::Dyn, nalgebra::VecStorage<f32, nalgebra::Const<1>, nalgebra::Dyn>> = path_cov.component_div(&path_len);

    for (i, ratio) in path_ratio.iter().enumerate() {
        gurobi_opt_var.hap_metrics[gurobi_opt_var.possible_paths_idx[i]].path_cov_ratio = Some(*ratio as f64);
    }

    if npaths > 0 {
        let mut x_binary_var_map = HashMap::new();
        let mut sum_x_binary = LinExpr::new();
        for i in 0..npaths {
            let var = add_var!(
                m,
                Binary,
                name: &format!("strain_indicator_{}", i)
            )?;
            x_binary_var_map.insert(i, var);
            sum_x_binary.add_term(1.0, var);
            m.add_constr(&format!("s_{}", i), c!(var >= (var_map[&i] - args.minimization_min_cov) * (1.0 / (2.0 * max_val))))?;
        };
        let _ = m.add_constr("total strains", c!(sum_x_binary <= npaths));
    }

    // only iter the node abundance > 0
    let valid_nodes: Vec<usize> = node_abundance_vec
        .iter()
        .enumerate()
        .filter(|&(_v, &ab)| ab > 0.0)
        .map(|(v, _)| v)
        .collect();

    let sample_valid_nodes = if args.sample_test {
        if valid_nodes.len() > 500 {
            debug!("\tFor {} species graph, subsample 500 nodes for test.", gurobi_opt_var.otu);
            sample_sorted(&valid_nodes, 500, 42)
        } else {
            valid_nodes
        }
    } else if !args.sample_test && args.sample_nodes > 0 {
        if valid_nodes.len() > args.sample_nodes {
            debug!("The {} species graph has too many nodes. Subsample {}.", gurobi_opt_var.otu, args.sample_nodes);
            sample_sorted(&valid_nodes, args.sample_nodes, 42)
        } else {
            valid_nodes
        }
    } else {
        valid_nodes
    };

    // // for debug
    // if gurobi_opt_var.otu == "Myxococcus_xanthus".to_string() {
    //     let file = File::create("valid_nodes.txt").unwrap();
    //     let mut writer = std::io::BufWriter::new(file);
    //     for (idx, &ab) in valid_nodes.iter().enumerate() {
    //         writeln!(writer, "{}\t{}", idx, node_abundance_vec[ab]).unwrap();
    //     }
    // }


    let mut y_var_vec = Vec::with_capacity(sample_valid_nodes.len());
    for &v in &sample_valid_nodes {
        let yv = add_var!(
            m,
            Continuous,
            name: &format!("y{}", v)
        )?;
        y_var_vec.push((v, yv));
    }
    m.update()?;

    // println!("{}\t{}", gurobi_opt_var.otu, y_var_vec.len());     
    let mut n_eval: usize = 0;
    for (i, (v, yv)) in y_var_vec.iter().enumerate() {
        let row = coeff_matrix.row(*v);
        let mut sumxv = LinExpr::new();
    
        for (j, val) in row.iter().enumerate() {
            if val.abs() > 1e-12 {
                if let Some(&var) = var_map.get(&j) {
                    sumxv.add_term(*val as f64, var);
                }
            }
        }
    
        let abundance = node_abundance_vec[*v];
    
        m.add_constr(&format!("y_{}_minus", v), c!( *yv >= sumxv.clone() - abundance ))?;
        m.add_constr(&format!("y_{}_plus", v), c!( *yv >= abundance - sumxv.clone() ))?;
    
        obj.add_term(1.0, *yv);
    
        n_eval += 1;
    }
    // assert_eq!(n_eval, 0);
    let obj_scaled = obj * (1.0 / n_eval as f64);
    m.set_objective(obj_scaled, grb::ModelSense::Minimize)?;
    m.update()?;
    m.set_param(param::LogToConsole, 0)?;
    m.set_param(param::Threads, args.gurobi_threads)?;
    m.set_param(param::NumericFocus, 0)?;
    m.set_param(param::PoolSearchMode, 0)?;
    m.set_param(param::PoolSolutions, 10)?;
    m.set_param(param::Method, 4)?;
    // m.set_param(param::Seed, 42)?;
    m.optimize()?;

    assert_eq!(m.status()?, Status::Optimal);
    let sols = m.get_obj_attr_batch(attr::X, x_vec.clone())?;
    let objval = m.get_attr(attr::ObjVal)?;
    debug!("\t\t{}\t{:?}\t{}", gurobi_opt_var.otu, sols, objval);

    for (i, sol) in sols.iter().enumerate() {
        gurobi_opt_var.hap_metrics[gurobi_opt_var.possible_paths_idx[i]].first_sol = Some(*sol);
    }

    let nstrains = sols.iter().filter(|&&x| x > 0.0).count();
    debug!("\t\tFirst optimization #strains / #paths = {} / {}", nstrains, npaths);

    second_filter_paths(gurobi_opt_var, &args);   
    if !gurobi_opt_var.second_opt {
        if args.debug { print!("\n"); }
        return Ok(());
    }

    debug!("\t\tSecond filter #strains / #paths = {} / {}", gurobi_opt_var.second_possible_paths_idx.len(), npaths); 

    m.reset(0)?;

    for (i, possible_path_idx) in gurobi_opt_var.possible_paths_idx.iter().enumerate() {
        if !gurobi_opt_var.second_possible_paths_idx.contains(possible_path_idx) {
            m.add_constr(&format!("x{}", i), c!(x_vec[i] == 0.0))?;
        }
    }

    m.optimize()?;

    assert_eq!(m.status()?, Status::Optimal);
    let sols2 = m.get_obj_attr_batch(attr::X, x_vec)?;
    let objval2 = m.get_attr(attr::ObjVal)?;
    debug!("\t\t{}\t{:?}\t{}", gurobi_opt_var.otu, sols2, objval2);

    let nstrains2 = sols2.iter().filter(|&&x| x > 0.0).count();
    debug!("\t\tSecond optimization #strains / #paths = {} / {}\n", nstrains2, npaths);

    for (&path_idx, sol) in gurobi_opt_var.possible_paths_idx.iter().zip(sols2) {
        if gurobi_opt_var.second_possible_paths_idx.contains(&path_idx) {
            if let Some(metric) = gurobi_opt_var.hap_metrics.get_mut(path_idx) {
                metric.second_sol = Some(sol);
            } else {
                panic!("path_idx {} out of bounds for hap_metrics (len = {})", path_idx, gurobi_opt_var.hap_metrics.len());
            }
        }
    }

    Ok(())
}


fn optimize_otu(args: &ProfileArgs, otu: &String, start: u32, end: u32, reads_cluster: &[Record]) -> Option<Vec<HapMetrics>> {
    log::debug!("Reading {} graph information, start: {}, end: {}", otu, start, end);
    let start = start - 1;
    let end = end - 1;
    let species_gfa_dir = if args.zip.is_some() {
        args.db.join("species_graph_info")
    } else {
        args.db.join("species_gfa")
    };
    
    if !species_gfa_dir.exists() {
        eprintln!("{:?} does not exist! Skipping.", species_gfa_dir);
        return None;
    }
    
    let gfa_file = match args.zip.as_deref() {
        Some("serialize") => species_gfa_dir.join(format!("{}.bin", otu)),
        Some("lz")        => species_gfa_dir.join(format!("{}.bin.lz4", otu)),
        Some("zstd")      => species_gfa_dir.join(format!("{}.bin.zst", otu)),
        Some("h5")        => species_gfa_dir.join(format!("{}.h5", otu)),
        _                 => species_gfa_dir.join(format!("{}.gfa", otu)),
    };

    // debug!("{:?}", gfa_file);
    
    let graph = match args.zip.as_deref() {
        Some("serialize") => load_from_zip_graph(&gfa_file, CompressType::Serialized)
            .map_err(|e| format!("GFA read error: {}", e)).ok()?,
        Some("lz")        => load_from_zip_graph(&gfa_file, CompressType::Lz4)
            .map_err(|e| format!("GFA read error: {}", e)).ok()?,
        Some("zstd")      => load_from_zip_graph(&gfa_file, CompressType::Zstd)
            .map_err(|e| format!("GFA read error: {}", e)).ok()?,
        #[cfg(feature = "h5")]
        Some("h5")      => load_from_zip_graph(&gfa_file, CompressType::Hdf5)
            .map_err(|e| format!("GFA read error: {}", e)).ok()?,            
        _                 => read_gfa(&gfa_file, 0)
            .map_err(|e| format!("GFA read error: {}", e)).ok()?,
    };

    // debug!("First 100 nodes_len: {:?}", &graph.nodes_len[..graph.nodes_len.len().min(100)]);

    let (unique_trio_nodes, unique_trio_nodes_len, hap2unique_trio_nodes_m) = trio_nodes_info(&graph);
    let (node_abundance_vec, trio_node_abundance_vec, node_base_cov) = get_node_abundances(otu, &graph.nodes_len, &unique_trio_nodes, &unique_trio_nodes_len, start as usize, &reads_cluster);
    let nvert = end - start + 1;
    let non_zero_count = node_abundance_vec.iter().filter(|&&x| x != 0.0).count();
    debug!("{} species node abundance > 0 number: {}", otu, non_zero_count);    
    let node_abundance_opt_vec: Vec<f64> = node_abundance_vec
        .iter()
        .map(|&x| if x > args.min_depth as f64 { x } else { 0.0 })
        .collect();

    // let hap_metrics: Vec<Vec<String>> = graph.paths.keys()
    //     .map(|hap| {
    //         let mut row = vec![hap.to_string()];
    //         row.extend(vec!["-".to_string(); 5]);
    //         row.push("0".to_string());
    //         row.push("-".to_string());
    //         row
    //     })
    //     .collect();
    // let possible_paths_idx = vec![0; graph.paths.len()];

    let mut gurobi_opt_var = GurobiOptVar {
        otu: otu.clone(),
        hap_metrics: vec![HapMetrics::default(); graph.paths.len()],
        possible_paths_idx: Vec::new(),
        second_possible_paths_idx: Vec::new(),
        orign_n_haps: 0,
        hap2trio_nodes_m_size: 0,
        same_path_flag: false,
        second_opt: false,
    };
    first_filter_paths(&mut gurobi_opt_var, &graph.paths, &hap2unique_trio_nodes_m, &trio_node_abundance_vec, &node_abundance_opt_vec, &args);
    gurobi_opt(&mut gurobi_opt_var, nvert as usize, &graph.paths, &node_abundance_vec, &node_base_cov, &graph.nodes_len, &args).unwrap();
    Some(gurobi_opt_var.hap_metrics)
}

fn abundace_constraint(species_profile_df: &DataFrame, metrics: &mut [HapMetrics]) -> Result<(), PolarsError> {
    let mut strain_absolute_abundance = Vec::new();
    for metric in metrics.iter_mut() {
        if metric.is_rescue == Some(true) && metric.first_sol.is_some() && metric.second_sol.is_some() {
            let first = metric.first_sol.unwrap();
            let second = metric.second_sol.unwrap();
            metric.second_sol = Some(first.min(second));
        }
        strain_absolute_abundance.push(metric.second_sol.unwrap_or(0.0));
    }
    let filtered_df = species_profile_df
        .clone()
        .lazy()
        .filter(col("species_taxid").eq(lit(metrics[0].otu.clone().unwrap())))
        .collect()?;

    let species_absolute_abundance = filtered_df
        .column("predicted_coverage")?
        .f64()? 
        .get(0).unwrap();

    let strain_absolute_abundance_sum: f64 = strain_absolute_abundance.iter().sum();
    let total_cov_diff = (strain_absolute_abundance_sum - species_absolute_abundance).abs() / ((strain_absolute_abundance_sum + species_absolute_abundance) / 2.0);
    for metric in metrics.iter_mut() {
        metric.total_cov_diff = Some(total_cov_diff);
    }

    if strain_absolute_abundance
        .iter()
        .cloned()
        .reduce(f64::max)
        .map_or(false, |max_val| max_val > 1.05 * species_absolute_abundance) 
    {
        let factor = species_absolute_abundance / strain_absolute_abundance_sum;
        for metric in metrics.iter_mut() {
            if !metric.is_rescue.unwrap_or(false) && metric.second_sol.is_some() {
                metric.second_sol = Some(metric.second_sol.unwrap()*factor);
            }
        }
    }

    Ok(())
}

// modify from https://stackoverflow.com/questions/73167416/creating-polars-dataframe-from-vecstruct
macro_rules! struct_to_dataframe {
    ($input:expr, [$($field:ident),+]) => {
        {
            let len = $input.len().to_owned();

            // Extract the field values into separate vectors
            $(let mut $field = Vec::with_capacity(len);)*

            for e in $input.into_iter() {
                $($field.push(e.$field.clone());)*
            }
            df! {
                $(stringify!($field) => $field,)*
            }
        }
    };
}

fn abundance_est(args: &ProfileArgs, hap_metrics_vec: &[HapMetrics]) -> Result<(), Box<dyn std::error::Error>> {
    let schema = Schema::from_iter(vec![
        Field::new("genome_ID".into(), DataType::String),
        Field::new("strain_taxid".into(), DataType::String),
        Field::new("species_taxid".into(), DataType::String),
        Field::new("organism_name".into(), DataType::String),
        Field::new("id".into(), DataType::String),
    ]);
    let genomes_metadata = LazyCsvReader::new(args.db.join("genomes_info.txt"))
        .with_has_header(true)
        .with_separator(b'\t')
        .with_schema(Some(Arc::new(schema)))
        .finish()?;
    let selected_genomes_metadata = genomes_metadata
        .select([
            col("genome_ID"),
            col("strain_taxid")
        ])
        .with_columns([
            col("genome_ID")
                .str()
                .split(lit("_"))
                .list()
                .slice(lit(0), lit(2))
                .list()
                .join(lit("_"), true)
                .alias("hap_id")
        ]);
    
    // let mut write_selected_genomes_metadata = selected_genomes_metadata.clone().collect()?;
    // save_output_to_file(&mut write_selected_genomes_metadata, "selected_genomes_metadata.txt", true)?;

    let hap_metrics_df = struct_to_dataframe!(hap_metrics_vec, [
        otu, hap_id, unique_trio_nodes_fraction, frequencies_mean, path_cov_ratio,
        first_sol, divergence, second_sol, is_rescue, total_cov_diff
    ]).unwrap();

    let hap_metrics_lazy_df = hap_metrics_df
        .lazy()
        .select([
            col("otu").alias("species_taxid"),
            col("hap_id"),
            col("unique_trio_nodes_fraction").alias("unique_trio_fraction"),
            col("frequencies_mean").alias("uniq_trio_cov_mean"),
            col("path_cov_ratio").alias("path_base_cov"),
            col("first_sol"),
            col("divergence").alias("strain_cov_diff"),
            col("second_sol").alias("predicted_coverage"),
            col("total_cov_diff")
        ]);
    
    // let mut write_hap_metrics_lazy_df = hap_metrics_lazy_df.clone().collect()?;
    // save_output_to_file(&mut write_hap_metrics_lazy_df, "write_hap_metrics_lazy_df.txt", true)?;

    let merged_df = hap_metrics_lazy_df
        .join(
            selected_genomes_metadata,
            [col("hap_id")],
            [col("hap_id")],
            JoinArgs::new(JoinType::Left),
        );

    let merged_df = merged_df.with_columns([
        (col("predicted_coverage") / col("predicted_coverage").sum())
            .alias("predicted_abundance")
    ]);

    let reordered_df = merged_df.clone().select([
        col("species_taxid"),
        col("strain_taxid"),
        col("genome_ID"),
        col("predicted_coverage"),
        col("predicted_abundance"),
        col("path_base_cov"),
        col("unique_trio_fraction"),
        col("uniq_trio_cov_mean"),
        col("first_sol"),
        col("strain_cov_diff"),
        col("total_cov_diff"),
    ]).collect()?;
    
    let mut write_reordered_df = reordered_df.clone();
    save_output_to_file(&mut write_reordered_df, "ori_strain_abundance.txt", true)?;

    let group_info = merged_df.clone()
        .group_by([col("species_taxid")])
        .agg([
            col("hap_id").count().alias("group_size"),
        ]);

    let filtered = merged_df
        .join(
            group_info,
            [col("species_taxid")],
            [col("species_taxid")],
            JoinArgs::new(JoinType::Inner),
        )
        .filter(
            // group_size > 1 || total_cov_diff <= single_cov_diff
            col("group_size").gt(lit(1))
                .or(col("total_cov_diff").lt_eq(lit(args.single_cov_diff))),
        )
        .filter(
            col("predicted_coverage")
                .gt_eq(lit(args.min_cov))
                .and(col("predicted_coverage").neq(lit(0.0))),
        )
        .with_columns([
            (col("predicted_coverage") / col("predicted_coverage").sum())
                .alias("predicted_abundance"),
        ]);    
    
    let sort_filtered_df = filtered
        .sort(["predicted_abundance"], SortMultipleOptions::new().with_order_descending(true));

    let mut final_df = if args.full {
        sort_filtered_df
            .select([
                // col("hap_id"),
                col("species_taxid"),
                col("strain_taxid"),
                col("genome_ID"),
                col("predicted_coverage"),
                col("predicted_abundance"),
                col("path_base_cov"),
                col("unique_trio_fraction"),
                col("uniq_trio_cov_mean"),
                col("first_sol"),
                col("strain_cov_diff"),
                col("total_cov_diff"),
            ])
            .collect()?
    } else {
        sort_filtered_df
            .select([
                // col("hap_id"),
                col("species_taxid"),
                col("strain_taxid"),
                col("genome_ID"),
                col("predicted_coverage").round(2),
                col("predicted_abundance"),
                col("path_base_cov").round(2),
                col("unique_trio_fraction").round(2),
                col("uniq_trio_cov_mean").round(2),
                col("first_sol").round(2),
                col("strain_cov_diff").round(2),
                col("total_cov_diff").round(2),
            ])
            .collect()?     
    };
        
    save_output_to_file(&mut final_df, &args.wd.join("strain_abundance.txt"), true)?;

    Ok(())
}

fn strain_profiling(args: &ProfileArgs, input_file: &InputFile, species_profile_df: &DataFrame, rcls_df: &DataFrame) -> Result<(), Box<dyn std::error::Error>> {
    let species_ranges = load_species_range(&args, &input_file, &species_profile_df)?;
    let read_clustered_by_species = group_reads_by_species(&rcls_df);
    let total = species_ranges.len();
    let counter = Arc::new(Mutex::new(0));
    
    let results: Vec<HapMetrics> = species_ranges
        .par_iter()
        .filter_map(|range| {
            let otu: &String = &range.species;
            let result: Option<Vec<HapMetrics>> = read_clustered_by_species
                .get(otu)
                .and_then(|reads_cluster| {
                    optimize_otu(&args, otu, range.start, range.end, reads_cluster)
                        .map(|mut metrics| {
                            let _ = abundace_constraint(&species_profile_df, &mut metrics);
                            metrics
                        } )
                });   
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
    // println!("{:?}", results);
    abundance_est(&args, &results)?;
    Ok(())
}

pub fn profile(args: ProfileArgs) -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let start_cpu = ProcessTime::now();
    
    let mut input_file = InputFile::default();
    check_args_valid(&args, &mut input_file);

    if args.species && !file_exists(&input_file.species_abund_file) {
        log::info!("- Read classification...");
        let rcls_df = rcls::rcls_profile(&input_file, args.threads)?;
        // let column_names = rcls_df.get_column_names();
        // println!("Column names: {:?}", column_names);
        if let Some(out_binning_file) = &args.out_binning_file {
            if out_binning_file != Path::new("None")  {
                let mut selected_rcls_df = rcls_df
                    .clone()
                    .lazy()
                    .select([
                        col("read_id"),
                        col("mapq"),
                        col("species"),
                        col("read_len")
                    ])
                    .collect()?;
                save_output_to_file(&mut selected_rcls_df, &out_binning_file, false)?;
            }
        }
        let filter_unmapped_rcls_df = rcls_df
            .clone()
            .lazy()
            .filter(col("species").neq(lit("U")))
            .collect()?;
    
        log::info!("- Species level profiling...");
        let species_profile_df = species_profiling(&args, &input_file, &filter_unmapped_rcls_df)?;

        if args.strain && !file_exists(&input_file.strain_abund_file) {
            log::info!("- Strain level profiling...");
            strain_profiling(&args, &input_file, &species_profile_df, &filter_unmapped_rcls_df)?;
        }
    } else if args.strain && !file_exists(&input_file.strain_abund_file) {
        log::info!("- Strain level profiling...");
        let schema = Schema::from_iter(vec![
            Field::new("read_id2".into(), DataType::String),
            Field::new("mapq2".into(), DataType::Int32),
            Field::new("species".into(), DataType::String),
            Field::new("read_len2".into(), DataType::Int32),
        ]);
        let reads_binning_df = LazyCsvReader::new(input_file.reads_binning_file.as_ref().unwrap())
            .with_has_header(false)
            .with_separator(b'\t')
            .with_schema(Some(Arc::new(schema)))
            .finish()?
            .select([col("species")])
            .collect()?;
        let gaf_df = load_gaf_file_lazy(&input_file.gaf_file.as_ref().unwrap())?;
        let rcls_df = concat_df_horizontal(
            &[gaf_df, reads_binning_df],
            false
        )?;

        // if args.debug {
        //     // Merge two GAF and read binning files, and then determine whether the read binning is out of order 
        //     // based on whether the two files have the same read_id. 
        //     // Because some fastq files have the same dual end read_id, 
        //     // using join merging may result in incorrect merging (Cartesian product).
        //     let same = rcls_df.clone().lazy().select([col("read_id").eq(col("read_id2")).all(false).alias("equal")]);
        //     let result = same.collect()?;
        //     let is_equal = result.column("equal")?.bool()?.get(0).unwrap();   
        //     let mut write_rcls_df = rcls_df.clone();
        //     save_output_to_file(&mut write_rcls_df, "reads_binning_debug.tsv", false)?;  
        //     println!("Equal: {is_equal}");       
        // }

        let filter_unmapped_rcls_df = rcls_df
            .clone()
            .lazy()
            .filter(col("species").neq(lit("U")))
            .collect()?;

        let schema = Schema::from_iter(vec![
            Field::new("species_taxid".into(), DataType::String),
            Field::new("predicted_abundance".into(), DataType::Float64),
            Field::new("predicted_coverage".into(), DataType::Float64),
        ]);
        let species_profile_df = LazyCsvReader::new(input_file.species_abund_file.as_ref().unwrap())
            .with_has_header(true)
            .with_schema(Some(Arc::new(schema)))
            .with_separator(b'\t')
            .finish()?
            .collect()?;

        strain_profiling(&args, &input_file, &species_profile_df, &filter_unmapped_rcls_df)?;

    } else {
        if args.species && !args.strain {
            info!("Species profiling abundance file: {:?} exists.", input_file.species_abund_file.as_ref().unwrap());
        } else if !args.species && args.strain {
            info!("Strain profiling abundance file: {:?} exists.", input_file.strain_abund_file.as_ref().unwrap());
        } else if args.species && args.strain {
            info!("Species and strain profiling abundance file: {:?}, {:?} both exist.", input_file.species_abund_file.as_ref().unwrap(), input_file.strain_abund_file.as_ref().unwrap());
        } 
    }

    let minutes_clock = start.elapsed().as_secs_f64() / 60.0;
    log::info!("- Profiling execution clock time: {:.2} min", minutes_clock);

    let minuted_cpu = start_cpu.elapsed().as_secs_f64() / 60.0;
    log::info!("- Profiling execution cpu time: {:.2} min", minuted_cpu);

    Ok(())
}