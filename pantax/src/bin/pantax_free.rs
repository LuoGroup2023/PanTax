
#[cfg(not(any(
    feature = "cp",
    feature = "gb",
    feature = "glpk",
    feature = "hs",
    feature = "free"
)))]
compile_error!("You must enable exactly one solver feature");

use anyhow::Ok;
use pantax::cli::Cli;
use pantax::cli::VgLevel;
use pantax::constants::*;
use pantax::types::*;
use pantax::construct;
use pantax::index;
use pantax::alignment;
use pantax::profile;
use pantax::utils::{init_global_config, get_config, GlobalConfig, get_vg_version, compare_version};
use clap::Parser;
use anyhow::Result;
use flexi_logger::style;
use flexi_logger::{DeferredNow, Duplicate, FileSpec};
use std::fs::File;
use std::panic;
use std::path::PathBuf;
use std::fs::{create_dir_all, remove_dir_all, copy, rename};
use chrono::Local;
use std::io::{BufRead, BufReader};

fn main() -> Result<()> {
    let mut args = Cli::parse();
    let data_type = initialize_setup(&mut args);
    install_panic_cleanup(&args);
    init_global_config(&args);
    let global_config = get_config();
    let check_points = check(&args, &global_config, &data_type);
    update_default_args(&mut args, &check_points);
    let mut genomes_metadata = read_genomes_info(args.genomes_info.as_ref().unwrap())?;
    if check_points.reconstruction {
        copy(args.genomes_info.as_ref().unwrap(), args.db.join(GENOMES_INFO))?;
        construct::construct(&mut genomes_metadata, &data_type, &check_points, &args, &global_config)?;
    }
    if check_points.need_index {
        index::index(&args, &global_config)?;
    }
    if check_points.alignemnt {
        alignment::alignemnt(&args, &global_config, &data_type)?;
    }
    if check_points.profiling {
        let profiling_config = get_profiling_config(&mut args, &global_config)?;
        profile::profile(profiling_config)?;
    }
    handle_res(&args, &global_config)?;
    log::info!("PanTax Done.");
    Ok(())
}

fn handle_res(args: &Cli, global_config: &GlobalConfig) -> Result<()> {
    let pwd = std::env::current_dir()?;
    if args.species {
        if args.pantax_output.is_some() {
            let specie_abund_file = PathBuf::from(format!("{}_species_abundance.txt", args.pantax_output.as_ref().unwrap()));
            copy(&global_config.species_abund_file, &specie_abund_file)?;
        } else {
            let specie_abund_file = pwd.join(&global_config.species_abund_file.file_name().unwrap());
            copy(&global_config.species_abund_file, &specie_abund_file)?;
        }
    }
    if args.strain {
        if args.pantax_output.is_some() {
            let strain_abund_file = PathBuf::from(format!("{}_strains_abundance.txt", args.pantax_output.as_ref().unwrap()));
            copy(&global_config.strain_abund_file, &strain_abund_file)?;
        } else {
            let strain_abund_file = pwd.join(&global_config.strain_abund_file.file_name().unwrap());
            copy(&global_config.strain_abund_file, &strain_abund_file)?;
        }
        if args.test {
            let ori_strain_abund_file = pwd.join(&global_config.ori_strain_abund_file.file_name().unwrap());
            copy(&global_config.strain_abund_file, &ori_strain_abund_file)?;            
        }
    }
    if (args.species && !args.next_for_strain) || (args.strain && args.next_for_strain) {
        if let Some(read_aln) = &args.read_aln {
            rename(&global_config.gfa_mapped, PathBuf::from(format!("{}.gaf", read_aln)))?;
        }
        if !args.debug {
            if let Some(pantax_report) = &global_config.pantax_report {
                if pantax_report.exists() {
                    let dest_pantax_report = pwd.join(pantax_report.file_name().unwrap());
                    rename(pantax_report, dest_pantax_report)?;
                }
            }
            remove_dir_all(&args.tmp_dir)?;
        }
    }

    Ok(())
}

fn get_profiling_config(args: &mut Cli, global_config: &GlobalConfig) -> Result<ProfilingConfig> {
    let range_file = File::open(&global_config.range_file_db)?;
    let reader = BufReader::new(range_file);
    let mut lines = reader.lines();
    let is_single_species = lines.next().is_none();
    if args.unique_trio_nodes_fraction.is_none() {
        if args.short_read {
            args.unique_trio_nodes_fraction = Some(0.3)
        } else if args.long_read {
            args.unique_trio_nodes_fraction = Some(0.5)
        }
    }
    if args.unique_trio_nodes_count.is_none() {
        args.unique_trio_nodes_count = Some(0.46)
    }

    let mut shift= false;
    if args.shift.is_none() {
        if is_single_species { shift = true }
    } else if args.shift.is_some() {
        if args.shift.as_ref().unwrap().to_lowercase() == "true".to_string() { shift = true }
    }

    let zip = if args.save {
        Some("serialize".to_string())
    } else if args.lz {
        Some("lz".to_string())
    } else if args.zstd {
        Some("zstd".to_string())
    } else {
        None
    };

    let profiling_config = ProfilingConfig {
        db: args.db.clone(),
        wd: args.tmp_dir.clone(),
        genomes_metadata: None,
        range_file: None,
        input_aln_file: Some(global_config.gfa_mapped.clone()),
        species_len_file: None,
        reads_binning_file: None,
        species_abund_file: None,
        out_binning_file: global_config.pantax_report.clone(),
        min_species_abundance: args.min_species_abundance,
        unique_trio_nodes_fraction: args.unique_trio_nodes_fraction.unwrap(),
        unique_trio_nodes_mean_count_f: args.unique_trio_nodes_count.unwrap(),
        single_cov_ratio: args.single_cov_ratio,
        single_cov_diff: args.single_cov_diff,
        minimization_min_cov: 0.,
        min_cov: args.min_cov,
        min_depth: args.min_depth,
        species: args.species,
        strain: args.strain,
        output_dir: args.tmp_dir.clone(),
        filtered: !args.no_filter,
        sample_nodes: args.sample_nodes,
        designated_species: args.designated_species.clone(),
        mode: if args.smode.is_some() { args.smode } else { Some(2) },
        full: true,
        solver: args.solver.clone(),
        gurobi_threads: args.gurobi_threads,
        shift,
        sample_test: args.sample_test,
        zip,
        force: args.force,
        trace: false,
        debug: args.debug,
    };
    Ok(profiling_config)
}

fn read_genomes_info(genomes_info_path: &PathBuf) -> Result<Vec<GenomesInfo>> {
    let genomes_info_file = File::open(genomes_info_path)?;
    let reader = BufReader::new(genomes_info_file);
    let mut genomes_metadata = Vec::new();
    for line in reader.lines().skip(1) {
        let line = line?;
        let tokens: Vec<&str> = line.trim().split('\t').collect();
        genomes_metadata.push(
            GenomesInfo {
                genome_id: tokens[0].to_string(),
                strain_taxid: tokens[1].to_string(),
                species_taxid: tokens[2].to_string(),
                organism_name: tokens[3].to_string(),
                path_id: tokens[4].to_string()
            }
        );
    }
    Ok(genomes_metadata)
}



fn check(args: &Cli, global_config: &GlobalConfig, data_type: &DataType) -> CheckPoints {
    let mut checkpoints = CheckPoints::default();
    let reference_pangenome_gfa = &global_config.reference_pangenome_gfa_db;
    if !reference_pangenome_gfa.is_file() {
        checkpoints.reconstruction = true;
    } else if reference_pangenome_gfa.is_file() && args.create {
        log::info!("Reference_pangenome.gfa exists. Skipping.");
    }

    if checkpoints.reconstruction && args.fast_query { 
        checkpoints.fast_query = true;
    }

    if args.query_and_filter {
        checkpoints.fast_query_and_filter = true
    }

    if (args.index || data_type.needs_index(&args)) && !&global_config.index_files[2].exists() {
        checkpoints.need_index = true;
    }

    if !args.create && !args.index && !args.query_and_filter {
        if !args.species && !args.strain {
            eprintln!("No species or strain option specified. Please set --species or --strain.");
            std::process::exit(1);
        }
        if !&global_config.gfa_mapped.exists() {
            checkpoints.alignemnt = true
        }
        if args.species || args.strain {
            checkpoints.profiling = true
        }
    }

    checkpoints
}

fn update_default_args(args: &mut Cli, check_points: &CheckPoints) {
    if check_points.reconstruction {
        if let Some(long_read_type) = &args.long_read_type {
            match long_read_type.as_str().to_lowercase().as_str() {
                "ontr9" | "ontr10" => args.ani = 85.,
                // hifi, ngs, or other
                _ => {}
            }
        }
    }
}

fn my_own_format_colored(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &flexi_logger::Record,
) -> Result<(), std::io::Error> {
    let mut paintlevel = record.level();
    if paintlevel == log::Level::Info {
        paintlevel = log::Level::Debug;
    }
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        style(paintlevel).paint(record.level().to_string()),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn my_own_format(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &flexi_logger::Record,
) -> Result<(), std::io::Error> {
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        record.level(),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn initialize_setup(args: &mut Cli) -> DataType {
    let log_spec = format!("{}", args.log_level_filter().to_string());
    let filespec = FileSpec::default()
        .basename("pantax")
        .o_directory(args.log_dir.clone())
        .o_discriminant(args.log_m.clone());
    flexi_logger::Logger::try_with_str(log_spec)
        .expect("Something went wrong with logging")
        .log_to_file(filespec) // write logs to file
        .duplicate_to_stderr(Duplicate::Info) // print warnings and errors also to the console
        .format(my_own_format_colored) // use a simple colored format
        .format_for_files(my_own_format)
        .start()
        .expect("Something went wrong with creating log file");   

    let cli_args: Vec<String> = std::env::args().collect();
    log::info!("COMMAND: {}", cli_args.join(" "));
    log::info!("VERSION: {}", env!("CARGO_PKG_VERSION"));
    // log::info!("SYSTEM NAME: {}", System::name().unwrap_or(format!("Unknown")));
    // log::info!("SYSTEM HOST NAME: {}", System::host_name().unwrap_or(format!("Unknown")));

    // check tools
    // TODO: need detail error output
    args.vg = which::which(args.vg.clone()).unwrap();
    args.pangenome_building_exe = which::which(args.pangenome_building_exe.clone()).unwrap();
    args.lr_aligner = which::which(args.lr_aligner.clone()).unwrap();
    if args.lr_aligner.to_string_lossy().contains("vg") {
        args.vg = args.lr_aligner.clone();
        // println!("vg {}", args.vg.display());
        // println!("lr aligner {}", args.lr_aligner.display());
    }
    let vg_version = get_vg_version(&args).unwrap();
    args.vg_lvl = match compare_version(&vg_version, "1.71.0") {
        0 | 1 => VgLevel::V2,
        -1 => VgLevel::V1,
        _ => unreachable!("compare_version returned unexpected value"),
    };

    if !args.fast_query && args.mode == Some(1) {
        args.fast_query = true;
    }

    if args.read_files.is_none() && !args.index && !args.create {
        eprintln!("read_files need to provided except only create index");
        std::process::exit(1);
    }

    if !args.create && !args.index {
        for file in args.read_files.as_ref().unwrap() {
            if !file.exists() {
                eprintln!("Input file does not exist: {:?}", file);
                std::process::exit(1);
            }
        }
    }
    
    let file_count = if let Some(read_files) = &args.read_files {
        read_files.len()
    } else {
        0
    };

    if args.fast_query && file_count == 0 {
        eprintln!("Fast mode need reads files input. Can not find the reads files");
        std::process::exit(1);
    }

    if args.short_read && args.long_read {
        eprintln!("Cannot specify both -s and -l");
        std::process::exit(1);
    } else if !args.short_read && !args.long_read && !args.create && !args.index && !args.fast_query {
        eprintln!("Should specify either -s or -l");
        std::process::exit(1);        
    }

    let mode = if args.long_read {
        if args.paired {
            eprintln!("ERROR [pantax] Long read mode cannot be paired (-l conflicts with -p)");
            std::process::exit(1);
        }
        if file_count != 1 {
            eprintln!("Long read mode requires exactly 1 input file");
            std::process::exit(1);
        }
        DataType::LongReadSingle
    } else if args.short_read {
        if args.paired {
            if file_count == 2 {
                DataType::ShortReadPaired
            } else if file_count == 1 {
                DataType::ShortReadPairedInter
            } else {
                eprintln!("Paired-end mode requires exactly 2 input files, got {}", file_count);
                std::process::exit(1);
            }
        } else {
            if file_count != 1 {
                eprintln!("Single-end mode requires exactly 1 input file, got {}", file_count);
                std::process::exit(1);
            }
            DataType::ShortReadSingle
        }
    } else if args.index || (args.create && !args.fast_query) {
        // not important, for only create or index
        DataType::ShortReadSingle
    } else {
        eprintln!("No read file.");
        std::process::exit(1);
    };
        
    // genomes_info
    let genomes_info = match &args.genomes_info {
        Some(path) if path.exists() => path.clone(),
        _ => {
            let fallback = args.db.join(GENOMES_INFO);
            if fallback.exists() {
                fallback
            } else {
                eprintln!(
                    "genomes information file {:?} does not exist.",
                    args.genomes_info.as_ref().unwrap_or(&fallback)
                );
                std::process::exit(1);
            }
        }
    };

    args.genomes_info = Some(genomes_info);
    
    // check db
    create_dir_all(&args.db)
        .unwrap_or_else(|e| panic!("Failed to create directory {:?}: {}", &args.db, e)); 

    // check tmp
    if args.tmp_dir == PathBuf::from("pantax_db_tmp") {
        let db_basename = args.db.file_name().unwrap().to_string_lossy();
        args.tmp_dir = PathBuf::from(format!("{}_tmp", db_basename));
    } 
    create_dir_all(&args.tmp_dir)
        .unwrap_or_else(|e| panic!("Failed to create directory {:?}: {}", &args.tmp_dir, e));

    mode
}

fn install_panic_cleanup(args: &Cli) {
    if args.debug {
        return;
    }

    let tmp_dir = args.tmp_dir.clone();
    let next_for_strain = args.next_for_strain;

    panic::set_hook(Box::new(move |info| {
        let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S");

        if !next_for_strain {
            eprintln!(
                "{} - PanTax encountered an exception while running and has deleted the tmp folder.",
                timestamp
            );
            eprintln!(
                "If you need to prevent deletion of the tmp directory, please add the --debug option."
            );

            if let Err(e) = remove_dir_all(&tmp_dir) {
                eprintln!("Failed to remove tmp dir {}: {}", tmp_dir.display(), e);
            }
        } else {
            eprintln!(
                "{} - PanTax encountered an exception.",
                timestamp
            );
        }

        eprintln!("panic info: {}", info);
    }));
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn print_args() {
        let args = Cli::parse_from([
            "pantax",
            "reads.fq.gz"
        ]);
        dbg!(args);

        let default_args = Cli::default();
        dbg!(default_args);
    }
}