

use pantax::cli::Cli;
use pantax::constants::*;
use pantax::types::*;
use pantax::construct;
use clap::Parser;
use anyhow::Result;
use flexi_logger::style;
use flexi_logger::{DeferredNow, Duplicate, FileSpec};
use std::fs::File;
use std::panic;
use std::path::PathBuf;
use std::fs::{create_dir_all, remove_dir_all};
use chrono::Local;
use std::io::{BufRead, BufReader};


fn main() -> Result<()> {
    let mut args = Cli::parse();
    let data_type = initialize_setup(&mut args);
    install_panic_cleanup(&args);
    let check_points = check(&args);
    update_default_args(&mut args, &check_points);
    let mut genomes_metadata = read_genomes_info(args.genomes_info.as_ref().unwrap())?;
    construct::construct(&mut genomes_metadata, &data_type, &check_points, &args)?;
    Ok(())
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



fn check(args: &Cli) -> CheckPoints {
    let mut checkpoints = CheckPoints::default();
    let reference_pangenome_gfa = args.db.join(REFERENCE_PANGENOME_GFA);
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

    if args.read_files.is_none() && !args.index {
        eprintln!("read_files need to provided except only create index");
        std::process::exit(1);
    }

    for file in args.read_files.as_ref().unwrap() {
        if !file.exists() {
            eprintln!("Input file does not exist: {:?}", file);
            std::process::exit(1);
        }
    }

    let file_count = args.read_files.as_ref().unwrap().len();

    if args.short_read && args.long_read {
        eprintln!("Cannot specify both -s and -l");
        std::process::exit(1);
    } else if !args.short_read && !args.long_read {
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
    } else {
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
    };
        


    // genomes_info
    if let Some(genomes_info) = &args.genomes_info {
        if !genomes_info.exists() {
            if !args.db.join(GENOMES_INFO).exists() {
                eprintln!("genomes information file {:?} does not exist.", genomes_info);
                std::process::exit(1);                
            } else {
                args.genomes_info = Some(args.db.join(GENOMES_INFO));
            }
        }
    }
    

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