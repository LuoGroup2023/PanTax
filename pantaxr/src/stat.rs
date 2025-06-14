use std::{fs::File, io::Write, path::PathBuf};
use std::fs;
use std::sync::Mutex;
use rayon::prelude::*;
use crate::cmdline::StatArgs;
use needletail::parse_fastx_file;
use polars::prelude::*;
use crate::rcls::save_output_to_file;
use log::*;
use crate::clog::init_logger;

struct Stats {
    total_number: usize,
    total_length: usize,
    gap_length: usize,
    avg_length: f64,
    n50: usize,
    n90: usize,
    max_length: usize,
    min_length: usize,
    gc_content: f64,
}

fn stat_single_fasta(filename: PathBuf) -> Stats {
    if !(filename.is_file() && filename.exists()) {
        panic!("Error: {:?} does not exists", filename);
    }

    let mut lengths = vec![];
    let mut total_gc = 0;
    let mut total_seq = 0;
    let mut total_n = 0;

    let mut reader = parse_fastx_file(filename).expect("valid path");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let seq = seqrec.seq();
        let len = seq.len();
        total_seq += len;
        total_gc += seq.iter().filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c').count();
        total_n += seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
        lengths.push(len);
    }

    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let avg = total_seq as f64 / lengths.len() as f64;
    let n50 = calc_nxx(&lengths, 0.5);
    let n90 = calc_nxx(&lengths, 0.9);
    let max = *lengths.first().unwrap_or(&0);
    let min = *lengths.last().unwrap_or(&0);
    let gc = 100.0 * total_gc as f64 / (total_seq - total_n) as f64;

    Stats {
        total_number: lengths.len(),
        total_length: total_seq,
        gap_length: total_n,
        avg_length: avg,
        n50,
        n90,
        max_length: max,
        min_length: min,
        gc_content: gc,
    }
}

fn calc_nxx(lengths: &[usize], fraction: f64) -> usize {
    let total: usize = lengths.iter().sum();
    let threshold = (total as f64 * fraction).ceil() as usize;
    let mut sum = 0;
    for &len in lengths {
        sum += len;
        if sum >= threshold {
            return len;
        }
    }
    0
}

fn print_stats_table(stats: &Stats) {
    println!("Total number (>0):         {}", stats.total_number);
    println!("Total length (bp):         {}", stats.total_length);
    println!("Gap number (bp):           {}", stats.gap_length);
    println!("Average length (bp):       {:.2}", stats.avg_length);
    println!("N50 length (bp):           {}", stats.n50);
    println!("N90 length (bp):           {}", stats.n90);
    println!("Maximum length (bp):       {}", stats.max_length);
    println!("Minimum length (bp):       {}", stats.min_length);
    println!("GC content (%):            {:.2}", stats.gc_content);
}

fn genome_metadata_stat(genome_metadata_file: PathBuf, species_len_output: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    let schema = Schema::from_iter(vec![
        Field::new("genome_ID".into(), DataType::String),
        Field::new("strain_taxid".into(), DataType::String),
        Field::new("species_taxid".into(), DataType::String),
        Field::new("organism_name".into(), DataType::String),
        Field::new("id".into(), DataType::String),
    ]);
    let genomes_metadata = LazyCsvReader::new(genome_metadata_file)
        .with_has_header(true)
        .with_separator(b'\t')
        .with_schema(Some(Arc::new(schema)))
        .with_truncate_ragged_lines(true)
        .finish()?
        .collect()?;

    let genomes_path = genomes_metadata.column("id")?.str()?.into_no_null_iter().collect::<Vec<_>>();

    let genomes_length: Vec<f64> = genomes_path
        .into_par_iter()
        .map(|path| {
            stat_single_fasta(path.into()).total_length as f64
        })
        .collect();

    let genomes_metadata = genomes_metadata.hstack(&[Column::new("genome_len".into(), genomes_length)])?;

    let mut species_avg_len_df = genomes_metadata
        .lazy()
        .group_by([col("species_taxid")])
        .agg([
            (col("genome_len").sum() / col("genome_len").len()).alias("avg_len")
        ])
        .collect()?;
    
    save_output_to_file(&mut species_avg_len_df, species_len_output, false)?;

    Ok(())
}

pub fn stat(args: StatArgs) -> Result<(), Box<dyn std::error::Error>> {

    init_logger("info");
    
    if let Some(genome) = args.input_file {
        let stats = stat_single_fasta(genome);
        print_stats_table(&stats);
    } else if let Some(genome_metadata_file) = args.genome_metadata_file {
        genome_metadata_stat(genome_metadata_file, &args.output_name)?;
        info!("- Genomes statistics results written to {}.", args.output_name.to_string_lossy());
    } else if let Some(genome_filelist) = args.input_list {
        let paths: Vec<String> = fs::read_to_string(genome_filelist).unwrap().lines().map(|s| s.to_string()).collect();
        let results = Mutex::new(vec![]);
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

        paths.par_iter().for_each(|path| {
            let stat = stat_single_fasta(path.into());
            let row = format!("{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\n",
                path,
                stat.total_number,
                stat.total_length,
                stat.gap_length,
                stat.avg_length,
                stat.n50,
                stat.n90,
                stat.max_length,
                stat.min_length,
                stat.gc_content);
            results.lock().unwrap().push(row);
        });

        let mut file = File::create(&args.output_name).unwrap();
        for line in results.into_inner().unwrap() {
            file.write_all(line.as_bytes()).unwrap();
        }
        info!("- Genomes statistics results written to {}.", args.output_name.to_string_lossy());

    }
        
    Ok(())
}