
use grb::expr::LinExpr;
use log::*;
use chrono::Local;
use std::fs::{File, metadata};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::collections::{BTreeMap, HashMap, HashSet};
use dashmap::DashMap;
use std::io::{BufRead, BufReader};
use regex::Regex;
use rayon::{prelude::*, range};
use polars::frame::DataFrame;
use polars::prelude::*;
use nalgebra::{DMatrix, RowDVector};
use grb::{prelude::*, Error};
use crate::{cmdline::ProfileArgs, rcls};
use rcls::save_output_to_file;
use std::io::{BufWriter, Write}; 

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

#[derive(Debug, Clone)]
pub struct Record {
    pub read_id: String,
    pub path: String,
    pub read_path_len: i64,
    pub read_start: i64,
    pub read_end: i64,
    pub species: String,
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

        records.push(Record {
            read_id: read_id.to_string(),
            path: path_series.get(i).unwrap().to_string(),
            read_path_len: read_path_len_series.get(i).unwrap(),
            read_start: read_start_series.get(i).unwrap(),
            read_end: read_end_series.get(i).unwrap(),
            species: species_series.get(i).unwrap().to_string(),
        });
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

pub struct Graph {
    pub nodes_len: Vec<i64>,
    pub paths: BTreeMap<String, Vec<usize>>,
}

pub fn read_gfa(gfa_file: &Path, previous: usize) -> Result<Graph, Box<dyn std::error::Error>> {
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

fn trio_nodes_info(graph: &Graph) -> (HashMap<(usize, usize, usize), usize>, Vec<i64>, DMatrix<u8>) {
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
    // let file = File::create("unique_trio_nodes.txt").unwrap();
    // let mut writer = BufWriter::new(file);
    // writeln!(writer, "i,j,k,value").unwrap();
    // for (&(i, j, k), &v) in unique_trio_nodes.iter() {
    //     writeln!(writer, "{},{},{},{}", i, j, k, v).unwrap();
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
        let target_len = read.read_end - read.read_start;
        let mut seen = 0;
        let mut read_nodes_len: HashMap<usize, i64> = HashMap::new();
        let mut undup_read_nodes = HashSet::new();

        // Set the initial base coverage of the node to 0    
        for &node in &read_nodes {
            read_nodes_len.insert(node, 0);
        }

        // only one node, read is shorter than the node of the graph
        if start_node == end_node && read_nodes.len() == 1 {
            *read_nodes_len.entry(start_node).or_default() += target_len;
            bases_per_node.entry(start_node).and_modify(|v| *v += target_len);

            if let Some(mut entry) = node_base_cov_info.get_mut(&start_node) {
                assert!(target_len as usize <= entry.2.len(), "read_id {} {:?} read len {} is longer than {} node len {}", read.read_id, read_nodes, target_len, start_node, entry.2.len() );
                for j in read.read_start as usize..read.read_end as usize {
                    entry.2[j] = 1;
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
                    (target_len - seen, 0)
                } else {
                    (node_len, 0)
                };

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
    // if otu.clone() == "Myxococcus_xanthus".to_string() {
    //     let file = File::create("trio_nodes_bases_count.txt").unwrap();
    //     let mut writer = BufWriter::new(file);
    //     for entry in trio_nodes_bases_count.iter() {
    //         writeln!(writer, "{},{}", entry.key(), entry.value()).unwrap();
    //     }
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

struct GurobiOptVar {
    otu: String,
    pub hap_metrics: Vec<Vec<String>>,
    pub possible_paths_idx: Vec<usize>,
}

fn filter_paths(
    // hap_metrics: &mut Vec<Vec<String>>,
    gurobi_opt_var: &mut GurobiOptVar,
    paths: &BTreeMap<String, Vec<usize>>,
    hap2trio_nodes_m: &DMatrix<u8>,
    trio_node_abundances: &Vec<f64>,
    node_abundance_vec: &Vec<f64>,
    args: &ProfileArgs,
) {
    
    // let mut filtered_haps = Vec::new();
    let orign_n_haps = paths.len();
    let hap2trio_nodes_m_size = hap2trio_nodes_m.len();
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
            
            // first metric
            gurobi_opt_var.hap_metrics[hap_idx][1] = format!("{:.2}", unique_trio_nodes_fraction);

            if args.shift {
                // filter outliers
                let filter_non_zero_frequencies = zscore_filter(&non_zero_frequencies, 3.0);
                let mean: f64 = if filter_non_zero_frequencies.is_empty() {
                    0.0
                } else {
                    filter_non_zero_frequencies.iter().sum::<f64>() / filter_non_zero_frequencies.len() as f64
                };
                let shift_unique_trio_nodes_fraction = if mean >= 1.0 {
                    let _shift_unique_trio_nodes_fraction = unique_trio_nodes_fraction + (0.8-unique_trio_nodes_fraction) * mean / 100.0;
                    if _shift_unique_trio_nodes_fraction > 0.8 {
                        0.8
                    } else {
                        _shift_unique_trio_nodes_fraction
                    }
                } else {
                    unique_trio_nodes_fraction * mean
                };
                debug!("\t\t{} unique trio node abundance > 0 ratio: {}, shift unique trio nodes fraction: {}, frequencies mean: {}", hap_id, unique_trio_nodes_fraction, shift_unique_trio_nodes_fraction, mean);
                
                if unique_trio_nodes_fraction < shift_unique_trio_nodes_fraction {
                    continue; // filter
                }

                // second metric
                gurobi_opt_var.hap_metrics[hap_idx][2] = format!("{:.2}", mean);
            } else {
                debug!("\t\t{} unique trio node abundance > 0 ratio: {}", hap_id, non_zero_frequencies.len());
                if unique_trio_nodes_fraction < args.unique_trio_nodes_fraction {
                    debug!("\t\t {}\t{}", gurobi_opt_var.otu, unique_trio_nodes_fraction);
                    continue; // filter
                }
                let filter_non_zero_frequencies = zscore_filter(&non_zero_frequencies, 3.0);
                let mean: f64 = if filter_non_zero_frequencies.is_empty() {
                    0.0
                } else {
                    filter_non_zero_frequencies.iter().sum::<f64>() / filter_non_zero_frequencies.len() as f64
                };          
                gurobi_opt_var.hap_metrics[hap_idx][2] = format!("{:.2}", mean);  
            }

            gurobi_opt_var.possible_paths_idx.push(hap_idx);
        };
        // filtered_haps.push(hap_id.clone());
    } else if orign_n_haps != 1 && hap2trio_nodes_m_size == 0 {
        let mut paths_vec = paths.values();
        let first_path = paths_vec.next().unwrap();
        let all_same = paths_vec.all(|x| x == first_path);
        if all_same {
            let non_zero_frequencies: Vec<f64> = node_abundance_vec
                .iter()
                .cloned() 
                .filter(|&x| x > 0.0)
                .collect(); 
            let mean: f64 = if non_zero_frequencies.is_empty() {
                0.0
            } else {
                non_zero_frequencies.iter().sum::<f64>() / non_zero_frequencies.len() as f64
            }; 
            gurobi_opt_var.possible_paths_idx.push(0);
        } else {
            warn!("Path more than 1, but have trio node is None.");
            gurobi_opt_var.possible_paths_idx = (0..orign_n_haps).collect();
        }

    } else if orign_n_haps == 1 {
        let non_zero_frequencies: Vec<f64> = node_abundance_vec
            .iter()
            .cloned() 
            .filter(|&x| x > 0.0)
            .collect(); 
        let mean: f64 = if non_zero_frequencies.is_empty() {
            0.0
        } else {
            non_zero_frequencies.iter().sum::<f64>() / non_zero_frequencies.len() as f64
        };
        gurobi_opt_var.possible_paths_idx.push(0);         
    } 
    // println!("end {:?}", gurobi_opt_var.possible_paths_idx);
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
        if gurobi_opt_var.possible_paths_idx.contains(&i){
            for &v in path.iter() {
                // assert!(v < nvert, "Node index {} exceeds nvert {}", v, nvert);
                // assert!(i < npaths, "Path index {} exceeds npaths {}, {:?}", i, npaths, gurobi_opt_var.possible_paths_idx);
                coeff_matrix[(v, i)] = 1.0; 
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
    let path_ratio = path_cov.component_div(&path_len);

    for (i, ratio) in path_ratio.iter().enumerate() {
        gurobi_opt_var.hap_metrics[gurobi_opt_var.possible_paths_idx[i]][3] = ratio.to_string();
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

    let valid_nodes: Vec<usize> = node_abundance_vec
        .iter()
        .enumerate()
        .filter(|&(_v, &ab)| ab > 0.0)
        .map(|(v, _)| v)
        .collect();

    // // for debug
    // if gurobi_opt_var.otu == "Myxococcus_xanthus".to_string() {
    //     let file = File::create("valid_nodes.txt").unwrap();
    //     let mut writer = std::io::BufWriter::new(file);
    //     for (idx, &ab) in valid_nodes.iter().enumerate() {
    //         writeln!(writer, "{}\t{}", idx, node_abundance_vec[ab]).unwrap();
    //     }
    // }


    let mut y_var_vec = Vec::with_capacity(valid_nodes.len());
    for &v in &valid_nodes {
        let yv = add_var!(
            m,
            Continuous,
            name: &format!("y{}", v)
        )?;
        y_var_vec.push((v, yv));
    }
    m.update()?;

    println!("{}\t{}", gurobi_opt_var.otu, y_var_vec.len());     
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
    m.set_param(param::Threads, 1)?;
    m.set_param(param::NumericFocus, 0)?;
    m.set_param(param::PoolSearchMode, 0)?;
    m.set_param(param::PoolSolutions, 10)?;
    m.set_param(param::Method, 1)?;
    m.set_param(param::Seed, 42)?;
    m.optimize()?;

    assert_eq!(m.status()?, Status::Optimal);
    let sol = m.get_obj_attr_batch(attr::X, x_vec.clone())?;
    let objval = m.get_attr(attr::ObjVal)?;
    println!("{}\t{:?}\t{}", gurobi_opt_var.otu, sol, objval);

    // m.reset(1)?;
    // m.optimize()?;

    // assert_eq!(m.status()?, Status::Optimal);
    // let sol2 = m.get_obj_attr_batch(attr::X, x_vec)?;
    // let objval2 = m.get_attr(attr::ObjVal)?;
    // println!("{}\t{:?}\t{}", gurobi_opt_var.otu, sol2, objval2);

    // println!("{:?}", sol);
    Ok(())
}

#[derive(Debug)]
pub struct OptimizeRes {
    pub otu: String,
    pub absolute_abund: f64,
}

fn optimize_otu(args: &ProfileArgs, otu: &String, start: u32, end: u32, reads_cluster: &[Record]) -> Option<Vec<OptimizeRes>> {
    log::debug!("Reading {} graph information, start: {}, end: {}", otu, start, end);
    let start = start - 1;
    let end = end - 1;
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
    let (node_abundance_vec, trio_node_abundance_vec, node_base_cov) = get_node_abundances(otu, &graph.nodes_len, &unique_trio_nodes, &unique_trio_nodes_len, start as usize, &reads_cluster);
    let nvert = end - start + 1;
    let non_zero_count = node_abundance_vec.iter().filter(|&&x| x != 0.0).count();
    debug!("{} species node abundance > 0 number: {}\n", otu, non_zero_count);    
    let node_abundance_opt_vec: Vec<f64> = node_abundance_vec
        .iter()
        .map(|&x| if x > args.min_depth as f64 { x } else { 0.0 })
        .collect();

    let hap_metrics: Vec<Vec<String>> = graph.paths.keys()
        .map(|hap| {
            let mut row = vec![hap.to_string()];
            row.extend(vec!["-".to_string(); 5]);
            row.push("0".to_string());
            row.push("-".to_string());
            row
        })
        .collect();
    // let possible_paths_idx = vec![0; graph.paths.len()];
    let possible_paths_idx = Vec::new();
    let mut gurobi_opt_var = GurobiOptVar {
        otu: otu.clone(),
        hap_metrics: hap_metrics,
        possible_paths_idx: possible_paths_idx,
    };
    filter_paths(&mut gurobi_opt_var, &graph.paths, &hap2unique_trio_nodes_m, &trio_node_abundance_vec, &node_abundance_opt_vec, &args);
    gurobi_opt(&mut gurobi_opt_var, nvert as usize, &graph.paths, &node_abundance_vec, &node_base_cov, &graph.nodes_len, &args);
    Some(vec![
        OptimizeRes { otu: otu.to_string(), absolute_abund: 10.0 }
    ])
}

pub fn profile(args: ProfileArgs) -> Result<(), Box<dyn std::error::Error>> {
    check_args_valid(&args);

    log::info!("{} - Read classification...", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let rcls_df = rcls::rcls_profile(&args)?;
    // let column_names = rcls_df.get_column_names();
    // println!("Column names: {:?}", column_names);
    let filter_unmapped_rcls_df = rcls_df
        .clone()
        .lazy()
        .filter(col("species").neq(lit("U")))
        .collect()?;

    log::info!("{} - Species level profiling...", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let species_profile_df = species_profile(&args, &filter_unmapped_rcls_df)?;

    log::info!("{} - Strain level profiling...", Local::now().format("%Y-%m-%d %H:%M:%S"));

    let species_ranges = load_species_range(&args, &species_profile_df)?;
    let read_clustered_by_species = group_reads_by_species(&rcls_df);
    let total = species_ranges.len();
    let counter = Arc::new(Mutex::new(0));
    
    let results: Vec<_> = species_ranges
        .par_iter()
        .filter_map(|range| {
            let otu: &String = &range.species;
            let result: Option<Vec<OptimizeRes>> = match read_clustered_by_species.get(otu) {
                Some(reads_cluster) => optimize_otu(&args, otu, range.start, range.end, reads_cluster),
                None => Some(vec![OptimizeRes { otu: otu.to_string(), absolute_abund: 0.0 }]),
            };    
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