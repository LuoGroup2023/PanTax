
use crate::cmdline::*;
use rayon::prelude::*;
use dashmap::DashMap;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    sync::{Arc, Mutex},
};

#[derive(Clone)]
struct AlignmentRecord {
    line: String,
    read_id: String,
    align_10: i32,      // Number of residue matches
    align_16: f64,      // Identity
    qual_12: i32,       // Mapping quality (0-60)
    span: i32,
}

fn parse_line(line: &str) -> Option<AlignmentRecord> {
    let fields: Vec<&str> = line.trim().split('\t').collect();
    if fields.len() < 16 {
        return None;
    }

    let read_id = fields[0].to_string();
    let align_10 = fields[9].parse::<i32>().ok()?;
    let align_16 = fields[15].rsplit(':').next()?.parse::<f64>().ok()?;
    let qual_12 = fields[11].parse::<i32>().ok()?;
    let span = fields[3].parse::<i32>().ok()? - fields[2].parse::<i32>().ok()?;

    Some(AlignmentRecord {
        line: line.to_string(),
        read_id,
        align_10,
        align_16,
        qual_12,
        span,
    })
}

pub fn filter_max_alignment_mt<P: AsRef<Path>>(gaf_file: P) -> std::io::Result<()> {
    let gaf_file = gaf_file.as_ref();
    let output_path = gaf_file.with_file_name(format!(
        "{}_filtered.gaf",
        gaf_file.file_stem().unwrap().to_string_lossy()
    ));

    let reader = BufReader::new(File::open(gaf_file)?);
    let all_lines: Vec<_> = reader.lines().filter_map(Result::ok).collect();

    // Parse + Collect all valid records in parallel
    let records: Vec<_> = all_lines
        .par_iter()
        .filter_map(|line| parse_line(line))
        .collect();

    // Step 1: Find best alignment per read_id using DashMap
    let best_map: DashMap<String, (i32, f64)> = DashMap::new();

    records.par_iter().for_each(|rec| {
        best_map
            .entry(rec.read_id.clone())
            .and_modify(|e| {
                if rec.align_10 > e.0 || (rec.align_10 == e.0 && rec.align_16 > e.1) {
                    *e = (rec.align_10, rec.align_16);
                }
            })
            .or_insert((rec.align_10, rec.align_16));
    });

    // Step 2: Filter and Write
    let writer = Arc::new(Mutex::new(BufWriter::new(File::create(&output_path)?)));
    let written_ids = Arc::new(DashMap::new());

    records
        .into_par_iter()
        .filter(|rec| rec.qual_12 > 20 && rec.span > 1000)
        .for_each(|rec| {
            if let Some(&(best_10, best_16)) = best_map.get(&rec.read_id).as_deref() {
                if rec.align_10 == best_10 && rec.align_16 == best_16 {
                    // Ensure only one line per read_id
                    if written_ids.insert(rec.read_id.clone(), ()).is_none() {
                        let mut writer = writer.lock().unwrap();
                        let _ = writeln!(writer, "{}", rec.line);
                    }
                }
            }
        });

    println!("Filtered GAF file written to: {}", output_path.display());
    Ok(())
}

pub fn gaf_filter(args: GafFilterArgs) -> Result<(), Box<dyn std::error::Error>> {
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    filter_max_alignment_mt(args.input_aln_file)?;
    Ok(())
}