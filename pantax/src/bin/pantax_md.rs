
use clap::Parser;
use std::process::Command;
use anyhow::Context;
use std::path::PathBuf;
use anyhow::Result;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Parser)]
#[command(
    author = "Wenhai Zhang", 
    about = "Strain-level metagenomic profiling using pangenome graphs with PanTax",
    arg_required_else_help = true
)]

#[derive(Default, Debug)]
struct MdCli {
    // input database directories for pantax
    #[arg(short='i', long="db", num_args = 1..)]
    db_dirs: Option<Vec<PathBuf>>,    

    // database list
    #[arg(short='l', long="db-list")]
    db_list: Option<PathBuf>,

    #[arg(short='w', long="wd", default_value = "pantax_md_db")]
    wd: PathBuf,    

    /// Number of processes to run in parallel.
    #[arg(short, long, default_value_t = 32)]
    pub threads: usize,

    /// Path to vg executable. Default is "vg" assuming it's in the system PATH.
    #[arg(short, long="vg", default_value = "vg")]
    pub vg: PathBuf,
}

struct DbFiles {
    genomes_info_file: PathBuf,
    species_gfa_dir: Option<PathBuf>,
    species_graph_info_dir: Option<PathBuf>,
    gfa: PathBuf,
    species_genomes_stats_file: PathBuf,
    species_range_file: PathBuf,
}

fn check(db_dirs: Vec<PathBuf>) -> Vec<DbFiles> {
    let mut db_files_vec: Vec<DbFiles> = Vec::new();
    for db_dir in db_dirs {
        if !db_dir.exists() {
            eprintln!("Error: Database directory does not exist: {}", db_dir.display());
            std::process::exit(1);
        }
        for file in ["genomes_info.txt", "reference_pangenome.gfa", "species_genomes_stats.txt", "species_range.txt"] {
            let file_path = db_dir.join(file);
            if !file_path.exists() {
                eprintln!("Error: Required file does not exist: {}", file_path.display());
                std::process::exit(1);
            }
        }
        let species_graph_info = db_dir.join("species_graph_info");
        let species_gfa = db_dir.join("species_gfa");
        if !species_graph_info.exists() && !species_gfa.is_dir() {
            eprintln!("Error: Required directory species_graph_info or species_gfa does not exist.");
            std::process::exit(1);
        }
        let db_files = DbFiles {
            genomes_info_file: db_dir.join("genomes_info.txt"),
            species_gfa_dir: Some(db_dir.join("species_gfa")),
            species_graph_info_dir: Some(db_dir.join("species_graph_info")),
            gfa: db_dir.join("reference_pangenome.gfa"),
            species_genomes_stats_file: db_dir.join("species_genomes_stats.txt"),
            species_range_file: db_dir.join("species_range.txt"),
        };
        db_files_vec.push(db_files);
    }
    db_files_vec
}

fn final_check(wd: &PathBuf) {
    let required_files = [
        wd.join("genomes_info.txt"),
        wd.join("species_genomes_stats.txt"),
        wd.join("species_range.txt"),
        wd.join("reference_pangenome.gfa"),
    ];
    for file in required_files.iter() {
        if !file.exists() {
            eprintln!("Error: Required file does not exist: {}", file.display());
            std::process::exit(1);
        }
    }
    let species_graph_info = wd.join("species_graph_info");
    let species_gfa = wd.join("species_gfa");
    if !species_graph_info.exists() && !species_gfa.is_dir() {
        eprintln!("Error: Required directory species_graph_info or species_gfa does not exist.");
        std::process::exit(1);
    }
}

fn merge_genomes_info_files(db_files_vec: &Vec<DbFiles>, wd: &PathBuf) -> Result<()> {
    let output = File::create(wd.join("genomes_info.txt"))?;
    let mut writer = BufWriter::new(output);
    for (i, genomes_info_file) in db_files_vec.iter().map(|db_files| &db_files.genomes_info_file).enumerate()  {
        let file = File::open(genomes_info_file)?;
        let reader = BufReader::new(file);

        for (line_idx, line) in reader.lines().enumerate() {
            let line = line?;
            if i > 0 && line_idx == 0 {
                continue;
            }

            writeln!(writer, "{}", line)?;
        }
    }
    Ok(())
}

fn merge_species_genomes_stats(db_files_vec: &Vec<DbFiles>, wd: &PathBuf) -> Result<()> {
    let output = File::create(wd.join("species_genomes_stats.txt"))?;
    let mut writer = BufWriter::new(output);
    for species_genomes_stats_file in db_files_vec.iter().map(|db_files| &db_files.species_genomes_stats_file)  {
        let file = File::open(species_genomes_stats_file)?;
        let reader = BufReader::new(file);

        for (_, line) in reader.lines().enumerate() {
            let line = line?;
            writeln!(writer, "{}", line)?;
        }
    }
    Ok(())
}

fn merge_gfa_info(db_files_vec: &Vec<DbFiles>, wd: &PathBuf) -> Result<()> {
    let is_zip = db_files_vec[0].species_graph_info_dir.is_some();
    if is_zip {
        let merge_species_graph_info_dir = wd.join("species_graph_info");
        std::fs::create_dir_all(&merge_species_graph_info_dir).expect("Failed to create species_graph_info directory");
        for species_graph_info_dir in db_files_vec.iter().map(|db_files| db_files.species_graph_info_dir.as_ref().unwrap()) {
            for entry in fs::read_dir(species_graph_info_dir)? {
                let entry = entry?;
                let path = entry.path();

                if path.is_file() {
                    let file_name = path.file_name().unwrap();
                    let new_path = merge_species_graph_info_dir.join(file_name);

                    fs::copy(&path, &new_path)?;
                }
            }            
        }
    } else {
        let merge_species_gfa = wd.join("species_gfa");
        std::fs::create_dir_all(&merge_species_gfa).expect("Failed to create species_gfa directory");
        for species_gfa_dir in db_files_vec.iter().map(|db_files| db_files.species_gfa_dir.as_ref().unwrap()) {
            for entry in fs::read_dir(species_gfa_dir)? {
                let entry = entry?;
                let path = entry.path();

                if path.is_file() {
                    let file_name = path.file_name().unwrap();
                    let new_path = merge_species_gfa.join(file_name);

                    fs::copy(&path, &new_path)?;
                }
            }            
        }
    }
    Ok(())
}

fn run_shell_command(command: &str) -> Result<()> {        
    let output = Command::new("bash")
        .arg("-c")
        .arg(command)
        .output()
        .with_context(|| format!("Failed to execute command: {}", command))?;
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        log::error!("Command failed with status: {:?}", output.status);
        log::error!("Command: {}", command);
        if !stderr.trim().is_empty() {
            log::error!("Stderr: {}", stderr);
        }
        if !stdout.trim().is_empty() {
            log::error!("Stdout: {}", stdout);
        }
        
        return Err(anyhow::anyhow!(
            "Command failed with exit code: {}\nCommand: {}",
            output.status.code().unwrap_or(-1),
            command
        ));
    }

    let stderr = String::from_utf8_lossy(&output.stderr);
    // let _stdout = String::from_utf8_lossy(&output.stdout);
    
    if !stderr.trim().is_empty() {
        log::warn!("Command succeeded with stderr output:\n{}", stderr);
    }
       
    Ok(())
}

fn sort_range(db_files_vec: &Vec<DbFiles>, wd: &PathBuf) -> Result<()> {
    let mut species_range: Vec<(String, usize, usize, usize, usize)> = Vec::new();
    for (i, species_range_file) in db_files_vec.iter().map(|db_files| &db_files.species_range_file).enumerate()  {
        let reader = BufReader::new(File::open(&species_range_file)?);
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.trim().split('\t').collect();
            assert_eq!(fields.len(), 4);

            let species = fields[0].to_string();
            let start: usize = fields[1].parse().unwrap();
            let end: usize = fields[2].parse().unwrap();  
            let is_pan: usize = fields[3].parse().unwrap();
        
            species_range.push((species, start, end, is_pan, i));
        }
    }
    let mut adjusted_results: Vec<(String, usize, usize, usize)> = Vec::new();
    let mut offset = 0;
    let mut id_end = 0;
    let mut id = 0;
    for (species, mut start, mut end, is_pan, i) in species_range.iter() {
        if *i != id {
            offset = id_end;
            id = *i;
        }
        start += offset;
        end += offset;
        adjusted_results.push((species.to_string(), start, end, *is_pan));
        id_end = end;
    }
    let range_file = wd.join("species_range.txt");
    let mut writer = BufWriter::new(File::create(range_file)?);
    for (taxid, start, end, is_pan) in adjusted_results.iter() {
        writeln!(writer, "{}\t{}\t{}\t{}", taxid, start, end, is_pan)?;
    }
    
    Ok(())
}

fn merge_gfa_files(db_files_vec: &Vec<DbFiles>, wd: &PathBuf, threads: usize, vg: &PathBuf) -> Result<()> {
    let mut vg_files = Vec::new();
    for (i, gfa_file) in db_files_vec.iter().map(|db_files| &db_files.gfa).enumerate() {
        let vg_path = wd.join(format!("{}.vg", i));
        let cmd = format!(
            "{} convert -g {} -p -t {} > {}",
            vg.display(),
            gfa_file.display(),
            threads,
            vg_path.display()
        );
        run_shell_command(&cmd)?;
        vg_files.push(vg_path.to_string_lossy().to_string());
    }
    let reference_pangenome_vg = wd.join("reference_pangenome.vg");
    let reference_pangenome_gfa = wd.join("reference_pangenome.gfa");
    let vg_files_str = vg_files.join(" ");   
    let cmd = format!(
        "{} combine -p {} > {}",
        vg.display(),
        vg_files_str,
        reference_pangenome_vg.display(),
    );
    run_shell_command(&cmd)?;

    for file in vg_files.iter() {
        if let Err(e) = fs::remove_file(file) {
            eprintln!("Failed to delete {}: {}", file, e);
        }
    }

    let cmd = format!(
        "{} convert -f {} -t {} > {}",
        vg.display(),
        reference_pangenome_vg.display(),
        threads,
        reference_pangenome_gfa.display(),
    );
    
    run_shell_command(&cmd)?;

    if let Err(e) = fs::remove_file(&reference_pangenome_vg) {
        eprintln!("Failed to delete {}: {}", reference_pangenome_vg.display(), e);
    }

    Ok(())
}

fn main() {
    let args = MdCli::parse();

    let db_dirs = if let Some(db_dirs) = args.db_dirs {
        db_dirs
    } else if let Some(db_list) = args.db_list {
        // read db_list and get db_dirs
        let mut db_dirs = Vec::new();
        let file = File::open(db_list).expect("Failed to open db list file");
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.expect("Failed to read line from db list file");
            db_dirs.push(PathBuf::from(line));
        }
        db_dirs
    } else {
        eprintln!("Error: Either --db or --db-list must be provided");
        std::process::exit(1);
    };
    let db_files_vec = check(db_dirs);
    std::fs::create_dir_all(&args.wd).expect("Failed to create working directory");
    merge_genomes_info_files(&db_files_vec, &args.wd).expect("Failed to merge genomes_info.txt files");
    merge_species_genomes_stats(&db_files_vec, &args.wd).expect("Failed to merge species_genomes_stats.txt files");
    merge_gfa_info(&db_files_vec, &args.wd).expect("Failed to merge species_graph_info or species_gfa directories");
    sort_range(&db_files_vec, &args.wd).expect("Failed to sort species_range.txt file");
    merge_gfa_files(&db_files_vec, &args.wd, args.threads, &args.vg).expect("Failed to merge GFA files");
    final_check(&args.wd);
}
