
use std::{collections::{HashMap, HashSet}, fs::{read_to_string, File}, io::{BufRead, BufReader, BufWriter, Write}, path::{Path, PathBuf}, process::{exit, Command}};

use crate::cmdline::*;

use std::sync::Mutex;

use rayon::prelude::*;

fn check_file_exists(path_opt: Option<&PathBuf>, name: &str) {
    match path_opt {
        Some(path) if path.is_file() => {}
        Some(path) => {
            eprintln!("Error: {} file {:?} does not exist or is not a regular file.", name, path);
            exit(1);
        }
        None => {
            eprintln!("Error: --{} must be provided.", name);
            exit(1);
        }
    }
}

fn check_args_valid(args: &SortRangeArgs) -> std::io::Result<()> {
    if !args.pangenome_vg_files.is_file() {
        eprintln!("Error: {:?} does not exist", args.pangenome_vg_files);
        exit(1);
    }
    if args.count_only {
        check_file_exists(args.pangenome_ge2.as_ref(), "pangenome_ge2");
        check_file_exists(args.pangenome_eq1.as_ref(), "pangenome_eq1");
        check_file_exists(args.vg.as_ref(), "vg");
    } else {
        check_file_exists(args.range_file.as_ref(), "range_file");
    }

    Ok(())
}

fn sort_range_only(pangenome_vg_files: PathBuf, range_file: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
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

    let vg_list = read_to_string(pangenome_vg_files).expect("Failed to read vg list file");
    let species_vec: Vec<String> = vg_list
        .trim()
        .split_whitespace()
        .map(|s| {
            Path::new(s)
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string()
        })
        .collect();   

    assert_eq!(species_range.len(), species_vec.len());

    let mut adjusted_results: Vec<(String, usize, usize, usize)> = Vec::with_capacity(species_vec.len());
    let mut offset = 0;
    for taxid in species_vec.iter() {
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

fn read_species_set(file: &PathBuf) -> HashSet<String> {
    let f = File::open(file).expect("Cannot open species file");
    let reader = BufReader::new(f);
    reader.lines().filter_map(Result::ok).collect()
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

fn get_node_counts_from_vg_files(
    vg_list_file: &PathBuf,
    vg: &PathBuf
) -> Vec<(String, usize)> {
    let vg_list = read_to_string(vg_list_file).expect("Failed to read vg list file");
    let vg_files: Vec<PathBuf> = vg_list.trim().split_whitespace().map(|s| s.into()).collect();

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
    ge2: &HashSet<String>,
    eq1: &HashSet<String>,
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


pub fn sort_range(args: SortRangeArgs) -> Result<(), Box<dyn std::error::Error>> {
    check_args_valid(&args)?;
    if args.count_only {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
        let ge2_set = read_species_set(&args.pangenome_ge2.unwrap());
        let eq1_set = read_species_set(&args.pangenome_eq1.unwrap());
    
        let counts = get_node_counts_from_vg_files(&args.pangenome_vg_files, &args.vg.unwrap());
    
        write_species_range(counts, &ge2_set, &eq1_set, &args.range_file.unwrap());
    } else {
        sort_range_only(args.pangenome_vg_files, args.range_file.unwrap())?;
    }
    
    Ok(())
}