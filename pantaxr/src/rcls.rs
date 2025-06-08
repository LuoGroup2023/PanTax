#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(dead_code)]
// use polars::lazy::dsl::*;
use polars::prelude::*;
// use rayon::prelude::*;
// use regex::Regex;
// use std::fs::File;
use rayon::prelude::*;
// use polars::df;
use regex::Regex;
use std::{fs::File, io::{BufRead, BufReader}, path::{self, Path, PathBuf}};
use crate::cmdline::{ProfileArgs, RclsArgs};

#[derive(Debug, Clone)]
struct SpeciesRange {
    species: String,
    start: u32,
    end: u32,
}

#[derive(Debug, Clone)]
struct GafRecord {
    read_id: String,
    length: u64,
    path: String,
    mapq: u8,
}

#[derive(Debug)]
struct ClassifiedRead {
    read_id: String,
    mapq: u8,
    species: Option<String>,
    length: u64,
}


// Step 1: Function to load species range file
fn load_species_range(file_path: &PathBuf) -> PolarsResult<Vec<(String, i64, i64)>> {
    let schema = Schema::from_iter(vec![
        Field::new("species".into(), DataType::String),
        Field::new("start".into(), DataType::Int64),
        Field::new("end".into(), DataType::Int64),
        Field::new("is_pan".into(), DataType::Int32),
    ]);
    let columns_to_select = ["species".into(), "start".into(), "end".into()];
    let species_df = CsvReadOptions::default()
        .map_parse_options(|parse_options| {
            parse_options.with_separator(b'\t')
        })
        .with_columns(Some(Arc::new(columns_to_select)))
        .with_schema(Some(Arc::new(schema)))
        .with_n_threads(Some(4))
        .with_has_header(false)
        .try_into_reader_with_file_path(Some(file_path.into()))?
        .finish()?;
    let species = species_df.column("species")?.str()?.into_no_null_iter();
    let starts = species_df.column("start")?.i64()?.into_no_null_iter();
    let ends = species_df.column("end")?.i64()?.into_no_null_iter();
    Ok(species
        .into_iter()
        .zip(starts.into_iter())
        .zip(ends.into_iter())
        .map(|((s, start), end)| (
            s.to_string(),        
            start,                 
            end,
        ))
        .collect())
}

fn read_species_range_file(path: &PathBuf) -> (Vec<SpeciesRange>, Vec<String>) {
    let file = File::open(path).expect("Cannot open species range file");
    let reader = BufReader::new(file);
    let mut ranges = Vec::new();
    let mut species_names = Vec::new();

    for line in reader.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split('\t').collect();
        let species = parts[0].to_string();
        let start = parts[1].parse().unwrap();
        let end = parts[2].parse().unwrap();
        ranges.push(SpeciesRange { species: species.clone(), start, end });
        species_names.push(species);
    }

    (ranges, species_names)
}


// Step 2: Function to load GAF file
fn load_gaf_file(file_path: &PathBuf, threads: usize) -> PolarsResult<DataFrame> {

    let columns_to_select = ["column_1".into(), "column_2".into(), "column_6".into(), "column_12".into()];
    let gaf_df = CsvReadOptions::default()
        .map_parse_options(|parse_options| {
            parse_options.with_separator(b'\t')
        })
        .with_columns(Some(Arc::new(columns_to_select)))
        .with_n_threads(Some(threads))
        .with_has_header(false)
        .try_into_reader_with_file_path(Some(file_path.into()))?
        .finish()?;
    let gaf_df_rename = gaf_df.lazy()
        .select([
            col("column_1").alias("read_id"),
            col("column_2").alias("length"),
            col("column_6").alias("path"),
            col("column_12").alias("mapq")
        ])
        .collect()?;
    // println!("{:?}", gaf_df_rename.head(Some(10)));
    Ok(gaf_df_rename)
}


pub fn load_gaf_file_lazy(file_path: &PathBuf) -> PolarsResult<DataFrame> {
    let gaf_df = LazyCsvReader::new(file_path.clone())
        .with_has_header(false)
        .with_separator(b'\t')
        .with_null_values(Some(NullValues::AllColumnsSingle("*".into())))
        .with_quote_char(None)
        .finish()?;
    let gaf_df_select = gaf_df
        .select([
            col("column_1").alias("read_id"),
            col("column_2").alias("read_len"),
            col("column_6").alias("path"),
            col("column_7").alias("read_path_len"),
            col("column_8").alias("read_start"),
            col("column_9").alias("read_end"),
            col("column_12").alias("mapq")
        ])
        .collect()?;
    // println!("{:?}", gaf_df_select.head(Some(10)));
    // println!("Number of rows: {}", gaf_df_select.height());
    // let filtered2 = gaf_df_select.clone()
    //     .lazy()
    //     .filter(col("read_id").eq(lit("S0R8221/1")))
    //     .collect()?;
    // println!("{:?}", filtered2);
    Ok(gaf_df_select)
}



fn load_gaf_file_as_structs(file_path: &PathBuf) -> std::io::Result<Vec<GafRecord>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();

    // let mut count = 0;
    for line in reader.lines() {
        // if count >= 10 {
        //     break;
        // }
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
    
        if fields.len() >= 13 {
            let record = GafRecord {
                read_id: fields[0].to_string(),
                length: fields[1].parse().unwrap_or(0),
                path: fields[5].to_string(),
                mapq: fields[11].parse().unwrap_or(0),
            };
            records.push(record);
            // count += 1;
        }
    }

    Ok(records)
}


fn classify_read(record: &GafRecord, species_ranges: &[SpeciesRange], node_re: &Regex) -> ClassifiedRead {
    
    let mut read_nodes: Vec<u32> = node_re
        .find_iter(&record.path)
        .filter_map(|m| m.as_str().parse::<u32>().ok())
        .collect();
    // println!("{:?}", read_nodes);
    if read_nodes.len() > 1 {
        let min = *read_nodes.iter().min().unwrap();
        let max = *read_nodes.iter().max().unwrap();
        read_nodes = vec![min, max];
    }
    
    let mut matched_species = None;
    for r in species_ranges {
        if read_nodes.iter().any(|&node| node >= r.start && node <= r.end) {
            matched_species = Some(r.species.clone());
            break;
        }
    }

    ClassifiedRead {
        read_id: record.read_id.clone(),
        mapq: record.mapq,
        species: matched_species,
        length: record.length,
    }
}

// Step 3: Function to process single read
fn process_single_read(
    read_name: &str,
    read_len: i64,
    node_path: &str,
    mapq: i64,
    species_info: &Vec<(String, i64, i64)>,
    regex: &Regex,
) -> (String, i64, String, i64) {
    let nodes: Vec<i64> = regex
        .find_iter(node_path)
        .filter_map(|m| m.as_str().parse::<i64>().ok())
        .collect();

    let (min, max) = match nodes.as_slice() {
        [] => (-1, -1),
        [n] => (*n, *n),
        _ => (*nodes.iter().min().unwrap(), *nodes.iter().max().unwrap()),
    };

    let matched_species = species_info
        .iter()
        .find(|(_, start, end)| min >= *start && max <= *end)
        .map(|(s, _, _)| s.clone())
        .unwrap_or_else(|| "U".to_string());
    (read_name.to_string(), mapq, matched_species, read_len)
}

fn process_single_read_simple(
    node_path: &str,
    species_info: &Vec<(String, i64, i64)>,
    regex: &Regex,
) -> String {
    let nodes: Vec<i64> = regex
        .find_iter(node_path)
        .filter_map(|m| m.as_str().parse::<i64>().ok())
        .collect();

    let (min, max) = match nodes.as_slice() {
        [] => (-1, -1),
        [n] => (*n, *n),
        _ => (*nodes.iter().min().unwrap(), *nodes.iter().max().unwrap()),
    };

    species_info
        .iter()
        .find(|(_, start, end)| min >= *start && max <= *end)
        .map(|(s, _, _)| s.clone())
        .unwrap_or_else(|| "U".to_string())
}

// Step 4: Function to process reads in parallel
fn process_reads_parallel(
    gaf_df: &DataFrame,
    species_info: &Vec<(String, i64, i64)>,
    regex: &Regex,
    threads: usize,
) -> PolarsResult<DataFrame> {
    let read_names = gaf_df.column("read_id")?.str()?.into_no_null_iter();
    let read_lens = gaf_df.column("length")?.i64()?.into_no_null_iter();
    let node_paths = gaf_df.column("path")?.str()?.into_no_null_iter();
    let mapqs = gaf_df.column("mapq")?.i64()?.into_no_null_iter();

    let results: Vec<(String, i64, String, i64)> = read_names
        .zip(read_lens)
        .zip(node_paths)
        .zip(mapqs)
        .par_bridge()
        .map(|(((name, len), path), mapq)| {
            process_single_read(name, len, path, mapq, species_info, regex)
        })
        .collect();

    let output_df = df![
        "read_id" => results.iter().map(|r| r.0.as_str()).collect::<Vec<_>>(),
        "mapq" => results.iter().map(|r| r.1).collect::<Vec<_>>(),
        "species" => results.iter().map(|r| r.2.as_str()).collect::<Vec<_>>(),
        "read_len" => results.iter().map(|r| r.3).collect::<Vec<_>>(),
    ]?;

    // let read_name_series = Series::new(
    //     "read_name".into(),
    //     results.iter().map(|r| r.0.as_str()).collect::<Vec<_>>(),
    // );

    // let mapq_series = Series::new("mapq".into(), results.iter().map(|r| r.1).collect::<Vec<_>>());
    // let species_series = Series::new(
    //     "species".into(),
    //     results.iter().map(|r| r.2.as_str()).collect::<Vec<_>>(),
    // );
    // let read_len_series = Series::new("read_len".into(), results.iter().map(|r| r.3).collect::<Vec<_>>());


    // Ok(DataFrame::new(vec![read_name_series.into(), mapq_series.into(), species_series.into(), read_len_series.into()])?)
    Ok(output_df)
}

fn process_reads_parallel_simple(
    gaf_df: &DataFrame,
    species_info: &Vec<(String, i64, i64)>,
    regex: &Regex,
    threads: usize,
) -> PolarsResult<DataFrame> {
    let node_paths = gaf_df.column("path")?.str()?.into_no_null_iter().collect::<Vec<_>>();
    assert_eq!(node_paths.len(), gaf_df.height(), "There may be null values in the path column of the gaf_df, node_paths {} should be equal to gaf_df {}", node_paths.len(), gaf_df.height());
    let matched_species: Vec<String> = node_paths
        .into_par_iter() 
        .map(|path| {
            process_single_read_simple(path, species_info, regex)
        })
        .collect();

    let output_df = gaf_df.hstack(&[Column::new("species".into(), matched_species)])?;

    Ok(output_df)
}


#[cfg(feature = "dev_test")]
fn process_single_read2(
    read_name: &str,
    read_len: i64,
    node_path: &str,
    mapq: i64,
    species_info: &Vec<(String, i64, i64)>,
    regex: &Regex,
) -> String {
    let nodes: Vec<i64> = regex
        .find_iter(node_path)
        .filter_map(|m| m.as_str().parse::<i64>().ok())
        .collect();

    let (min, max) = match nodes.as_slice() {
        [] => (-1, -1),
        [n] => (*n, *n),
        _ => (*nodes.iter().min().unwrap(), *nodes.iter().max().unwrap()),
    };

    let matched_species = species_info
        .iter()
        .find(|(_, start, end)| min >= *start && max <= *end)
        .map(|(s, _, _)| s.clone())
        .unwrap_or_else(|| "Unmapped".to_string());

        matched_species
}


#[cfg(feature = "dev_test")]
fn process_reads_parallel2(
    gaf_df: &DataFrame,
    species_info: &Vec<(String, i64, i64)>,
    regex: &Regex,
    threads: usize,
) -> PolarsResult<DataFrame> {
    let species_info_arc = Arc::new(species_info.clone());
    let regex_arc = Arc::new(regex.clone());
    gaf_df
        .clone()
        .lazy()
        .select([as_struct(vec![col("read_id"), col("length"), col("path"), col("mapq")]).map(
            move |s| {
                let gaf_struct = s.struct_()?;
                let read_id_series = gaf_struct.field_by_name("read_id")?;
                let chunk_read_id = read_id_series.str()?;
                
                let length_series = gaf_struct.field_by_name("length")?;
                let chunk_length = length_series.i64()?;
                
                let path_series = gaf_struct.field_by_name("path")?;
                let chunk_path = path_series.str()?;
                
                let mapq_series = gaf_struct.field_by_name("mapq")?;
                let chunk_mapq = mapq_series.i64()?;

                let species: StringChunked = chunk_read_id
                .into_iter()
                .zip(chunk_length.into_iter())
                .zip(chunk_path.into_iter())
                .zip(chunk_mapq.into_iter())
                .map(|(((name, len), path), mapq)| match (name, len, path, mapq) {
                    (Some(name), Some(len), Some(path), Some(mapq)) => {
                        process_single_read2(
                            name,
                            len,
                            path,
                            mapq,
                            &species_info_arc.clone(),
                            &regex_arc.clone(),
                        )
                    }
                    _ => String::from("None"),
                })
                .collect();
                Ok(Some(species.into_column()))
            }, GetOutput::from_type(DataType::String),
        )])
        .collect()
}

// Step 5: Function to save output DataFrame to file
pub fn save_output_to_file(
    df: &mut DataFrame,
    output_file: impl AsRef<Path>,
    write_with_header: bool
) -> PolarsResult<()> {
    let file_path = output_file.as_ref();
    CsvWriter::new(File::create(file_path)?)
        .with_separator(b'\t')
        .include_header(write_with_header)
        .finish(df)?;
    Ok(())
}

// use std::collections::HashMap;
// fn dataframe_equal_on_key(df1: &DataFrame, df2: &DataFrame, key: &str) -> PolarsResult<bool> {
//     let key_idx = df1
//         .get_column_names()
//         .iter()
//         .position(|name| name.as_str() == key)
//         .ok_or_else(|| PolarsError::ComputeError(format!("Key column '{}' not found", key).into()))?;

//         fn df_to_map(df: &DataFrame, key_idx: usize) -> PolarsResult<HashMap<String, Vec<AnyValue<'_>>>> {
//             let columns = df.get_columns();
//             let nrows = df.height();
//             let mut row_map = HashMap::new();
    
//             for i in 0..nrows {
//                 let row: Vec<AnyValue> = columns.iter()
//                     .map(|col| col.get(i).unwrap_or(AnyValue::Null))
//                     .collect();
//                 let key_val = row[key_idx].to_string();
//                 row_map.insert(key_val, row);
//             }
//             Ok(row_map)
//         }
    
//         let df1_rows = df_to_map(df1, key_idx)?;
//         let df2_rows = df_to_map(df2, key_idx)?;

//     Ok(df1_rows == df2_rows)
// }

use crate::profile::InputFile;
pub fn rcls_profile(input_file: &InputFile, threads: usize) -> PolarsResult<DataFrame> {
    let species_info = load_species_range(&input_file.range_file.as_ref().unwrap())?;
    let gaf_df = load_gaf_file_lazy(&input_file.gaf_file.as_ref().unwrap())?;
    let regex = Regex::new(r"\d+").unwrap();
    let output_df = process_reads_parallel_simple(&gaf_df, &species_info, &regex, threads)?;    
    Ok(output_df)
}

pub fn rcls(args: RclsArgs) -> Result<(), Box<dyn std::error::Error>> {
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    
    if args.struct_read {
        let (species_ranges, _species_names) = read_species_range_file(&args.range_file);
        let reads = load_gaf_file_as_structs(&args.input_aln_file)?;
        let node_re = Regex::new(r"\d+").unwrap();
        let results: Vec<_> = reads.par_iter()
            .map(|r| classify_read(r, &species_ranges, &node_re))
            .collect();
    
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(&args.output_name).unwrap();
    
        for r in results {
            wtr.write_record(&[
                &r.read_id,
                &r.mapq.to_string(),
                &r.species.unwrap_or_else(|| "NA".to_string()),
                &r.length.to_string(),
            ]).unwrap();
        }
        wtr.flush().unwrap();

    } else if args.lazy_read {
        let species_info = load_species_range(&args.range_file)?;
        let gaf_df = load_gaf_file_lazy(&args.input_aln_file)?;
        let regex = Regex::new(r"\d+").unwrap();
        let mut output_df = process_reads_parallel(&gaf_df, &species_info, &regex, args.threads)?;

        // let mut output_df = process_reads_parallel2(&gaf_df, &species_info, &regex, args.threads)?;

        // test rust result is equal to old python result
        // let output_df_ori = process_reads_parallel(&gaf_df, &species_info, &regex, args.threads)?;
        // let output_df_ori_clone = output_df_ori.clone();
        // let mut output_df = output_df_ori;
        // if Path::new("reads_classification.tsv").exists() {
        //     let schema = Schema::from_iter(vec![
        //         Field::new("column_1".into(), DataType::String),
        //         Field::new("column_2".into(), DataType::Int64),
        //         Field::new("column_3".into(), DataType::String),  
        //         Field::new("column_4".into(), DataType::Int64),
        //     ]);
        //     let python_reads_classification_df = LazyCsvReader::new(Path::new("reads_classification.tsv"))
        //         .with_has_header(false)
        //         .with_separator(b'\t')
        //         .with_schema(Some(Arc::new(schema)))
        //         .finish()?;
        //     let python_reads_classification_df_rename = python_reads_classification_df
        //         .select([
        //             col("column_1").alias("read_id"),
        //             col("column_2").alias("mapq"),
        //             col("column_3").alias("species"),
        //             col("column_4").alias("read_len"),
        //         ])
        //         .with_columns([
        //             col("species").fill_null(lit("U")),
        //         ])
        //         .collect()?;
        //     let are_equal = dataframe_equal_on_key(&output_df_ori_clone, &python_reads_classification_df_rename, "read_id")?;
        //     if !are_equal {
        //         println!("The dataframes are not equal!");
        //     }
        // }
        save_output_to_file(&mut output_df, &args.output_name, false)?;
    } else {
        let species_info = load_species_range(&args.range_file)?;
        let gaf_df = load_gaf_file(&args.input_aln_file, args.threads)?;
        let regex = Regex::new(r"\d+").unwrap();
        let mut output_df = process_reads_parallel(&gaf_df, &species_info, &regex, args.threads)?;
        save_output_to_file(&mut output_df, &args.output_name, false)?;
    }
    
    Ok(())
}

