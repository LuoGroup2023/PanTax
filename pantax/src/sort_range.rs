
use crate::cli::Cli;
use std::{collections::HashMap, fs::File, io::{BufRead, BufReader, BufWriter, Write}, path::PathBuf, process::Command};
use std::sync::Mutex;
use rayon::prelude::*;
use anyhow::Result;

fn sort_range_only(all_multispecies: &Vec<String>, all_single_species: &Vec<String>, range_file: &PathBuf) -> Result<()> {
    let mut species_range: HashMap<String, (usize, usize, usize)> = HashMap::new();
    let reader = BufReader::new(File::open(&range_file)?);

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.trim().split('\t').collect();
        assert_eq!(fields.len(), 4);

        let species = fields[0].to_string();
        let start: usize = fields[1].parse().unwrap();
        let end: usize = fields[2].parse().unwrap();  
        let is_pan: usize = fields[3].parse().unwrap();
    
        species_range.insert(species, (start, end, is_pan));
    }

    let mut adjusted_results: Vec<(String, usize, usize, usize)> = Vec::new();
    let mut offset = 0;
    for taxid in all_multispecies.iter().chain(all_single_species) {
        let (mut start, mut end, is_pan) = species_range[taxid];
        start += offset;
        end += offset;
        adjusted_results.push((taxid.clone(), start, end, is_pan));
        offset = end;
    }

    let mut writer = BufWriter::new(File::create(range_file)?);
    for (taxid, start, end, is_pan) in adjusted_results.iter() {
        writeln!(writer, "{}\t{}\t{}\t{}", taxid, start, end, is_pan)?;
    }
    
    Ok(())
}


fn get_node_count(vg_file: &PathBuf, vg: &PathBuf) -> Option<(String, usize)> {
    if !vg_file.is_file() {
        eprintln!("Input file not found: {}", vg_file.display());
        return None;
    }

    let output = Command::new(vg)
        .args(["stats", "-N"])
        .arg(vg_file)
        .output()
        .ok()?;

    if !output.status.success() {
        eprintln!(
            "vg stats -N {} failed: {}",
            vg_file.display(),
            String::from_utf8_lossy(&output.stderr)
        );
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let node_count: usize = match stdout.trim().parse() {
        Ok(n) => n,
        Err(_) => {
            eprintln!("Failed to parse node count from output: {}", stdout.trim());
            return None;
        }
    };

    let species = vg_file
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    Some((species, node_count))
}

fn get_node_counts_from_vg_files(all_multispecies: &Vec<String>, all_single_species: &Vec<String>, gfa_build_wd: &PathBuf, vg: &PathBuf) -> Vec<(String, usize)> {
    let mut vg_files = Vec::new();
    for species in all_multispecies.iter().chain(all_single_species) {
        let vg_file = gfa_build_wd.join(format!("{}.vg", species));
        vg_files.push(vg_file);
    }

    let results = Mutex::new(vec![None; vg_files.len()]);
    vg_files
        .par_iter()
        .enumerate()
        .for_each(|(i, vg_file)| {
            if let Some((species, count)) = get_node_count(vg_file, vg) {
                results.lock().unwrap()[i] = Some((species, count));
            }
        });

    results
        .into_inner()
        .unwrap()
        .into_iter()
        .filter_map(|x| x)
        .collect()
}

fn write_species_range(
    counts: Vec<(String, usize)>,
    ge2: &Vec<String>,
    eq1: &Vec<String>,
    output_file: &PathBuf,
) {
    let mut fout = File::create(output_file).expect("Failed to create output file");
    let mut nodes_sum = 0;

    for (species, count) in counts {
        let flag = if ge2.contains(&species) {
            1
        } else if eq1.contains(&species) {
            0
        } else {
            eprintln!("Species not found in either list: {}", species);
            continue;
        };

        let start = nodes_sum + 1;
        let end = nodes_sum + count;
        writeln!(fout, "{}\t{}\t{}\t{}", species, start, end, flag).unwrap();
        nodes_sum = end;
    }
}

pub fn sort_range(all_multispecies: &Vec<String>, all_single_species: &Vec<String>, gfa_build_wd: &PathBuf, range_file: &PathBuf, args: &Cli) -> Result<()> {
    log::info!("Obtain the mapping information of graph and species...");
    if args.is_save() {    
        sort_range_only(&all_multispecies, &all_single_species, &range_file)?;
    } else {
        let counts = get_node_counts_from_vg_files(&all_multispecies, &all_single_species, &gfa_build_wd, &args.vg);
        write_species_range(counts, &all_multispecies, &all_single_species, &range_file);
    }
    log::info!("Obtain the mapping information of graph and species...done");
    Ok(())
}