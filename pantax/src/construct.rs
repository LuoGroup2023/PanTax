
use crate::cli::Cli;
use crate::types::*;
use crate::constants::*;
use crate::utils::{output_genomes_info, run_shell_command, GlobalConfig};
use crate::task_scheduling::task_schedule_multi;
use crate::build_eq1::{parallel_build_gfa, parallel_convert_and_zip};
use crate::sort_range::sort_range;
use crate::stat::stat;
use sylph::api::contain_new;
use std::{collections::HashSet, path::Path};
use std::fs::{self, create_dir_all, File, remove_dir_all, remove_file};
use anyhow::{Result, anyhow};
// use rayon::prelude::*;
use needletail::parse_fastx_file;
use fxhash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use rayon::prelude::*;

pub fn construct(genomes_metadata: &mut Vec<GenomesInfo>, data_type: &DataType, check_points: &CheckPoints, args: &Cli, global_config: &GlobalConfig) -> Result<()> {
    if check_points.fast_query {
        log::info!("Sylph query and filter with ANI threshold {}...", args.ani);
        create_dir_all(args.tmp_dir.join("sylph_db"))?;
        if !&global_config.query_done.exists() {
            let reference_genomes_paths: Vec<String> = genomes_metadata.iter().map(|x| x.path_id.clone()).collect();
            let read_files: Vec<String> = args.read_files.as_ref().unwrap().iter().map(|x| x.to_string_lossy().into_owned()).collect();
            let is_interleaved = match data_type {
                DataType::ShortReadPairedInter => true,
                _ => false
            };
            // println!("read_files: {:?}", read_files);
            // println!("is_interleaved: {}", is_interleaved);
            let out_prefix = Some(global_config.query_res.to_string_lossy().into_owned());
            let mut query_res = contain_new(&reference_genomes_paths, &read_files, &out_prefix, is_interleaved);
            query_res.retain(|x| x.2 >= args.ani);
            // for (_, s, _) in query_res.iter_mut() {
            //     if let Some(first) = s.split_whitespace().next() {
            //         *s = first.to_string();
            //     }
            // }
            // let seq2genome: FxHashMap<String, String> = genomes_metadata
            //     .par_iter()
            //     .flat_map(|genomes_info| process_genome(&genomes_info.genome_id, &genomes_info.path_id))
            //     .collect();
            let retained_genomes_path: Vec<String> = query_res.iter().map(|x| x.0.clone()).collect::<HashSet<String>>().into_iter().collect();
            genomes_metadata.retain(|x| retained_genomes_path.contains(&x.path_id));
            output_genomes_info(&genomes_metadata, &args.tmp_dir.join("sylph_db/filter_genomes_info.txt"))?;
            output_genomes_info(&genomes_metadata, &args.db.join("genomes_info.txt"))?;
            File::create(&global_config.query_done)?;
            log::info!("Sylph query and filter with {}...done", args.ani);
        } else {
            log::info!("Sylph query and filter with {} has been done. Skipping.", args.ani);
        }
        if check_points.fast_query_and_filter && !args.debug {
            remove_dir_all(&args.tmp_dir)?;
        }
    }
    let pangeno_dir = &global_config.pangeno_dir;
    std::fs::create_dir_all(&pangeno_dir)?;
    let gfa_build_wd = &global_config.gfa_build_wd;
    std::fs::create_dir_all(&gfa_build_wd)?;
    let single_species_file = &global_config.single_species_file;
    let (single_species2genome, multi_species2genomes) = prepare_species_genomes(genomes_metadata, pangeno_dir, single_species_file)?;
    // if args.debug {
    //     let _ = write_multi_species2genomes(&multi_species2genomes, &pangeno_dir);
    // }
    log::info!("Species pangenome need to build: {}", multi_species2genomes.len());

    // multispecies
    let all_multispecies = if !multi_species2genomes.is_empty() {
        task_schedule_multi(&args, &multi_species2genomes)?
    } else {
        log::warn!("No species has more than 2 genomes!");
        Vec::new()
    };
    
    // single species
    let mut all_single_species = Vec::new();
    if !single_species2genome.is_empty() {
        let gfa_paths = parallel_build_gfa(&single_species2genome, gfa_build_wd, CHUNK_SIZE);
        let done_path = &global_config.gfa_done_path;
        let mut done_writer = BufWriter::new(File::create(done_path).expect("Cannot write done file"));
        for path in &gfa_paths {
            writeln!(done_writer, "{}", path).unwrap();
            let species = PathBuf::from(path).file_stem().unwrap().to_string_lossy().to_string();
            all_single_species.push(species);
        }
        if !args.debug {
            remove_dir_all(pangeno_dir)?;
            remove_file(single_species_file)?;
        }
        parallel_convert_and_zip(&gfa_paths, args)?;    
    } else {
        log::info!("No species has 1 genome.");
    }
    let all_species: Vec<&String> = all_multispecies
        .iter()
        .chain(all_single_species.iter()).collect();
    vg_combine(&all_species, gfa_build_wd, &global_config.reference_pangenome_vg, &args.vg)?;
    vg_convert_pan(&args.vg, args.threads, &global_config.reference_pangenome_vg, &global_config.reference_pangenome_gfa)?;
    sort_range(&all_multispecies, &all_single_species, gfa_build_wd, &global_config.range_file, &args)?;
    if !args.debug {
        remove_file(&global_config.reference_pangenome_vg)?;
    }
    handle_build_res(&args, &global_config, all_species)?;
    stat(&global_config.genomes_metadata, &global_config.species_genomes_stats)?;
    log::info!("Building reference pangenome completely.");
    if args.create && !args.index {
        if !args.debug {
            remove_dir_all(&args.tmp_dir)?;
        } 
        std::process::exit(0);
    }
    Ok(())
}

fn handle_build_res(args: &Cli, global_config: &GlobalConfig, all_species: Vec<&String>) -> Result<()> {
    if args.is_save() {
        fs::create_dir_all(&global_config.species_graph_info)?;
        move_dir_all(&global_config.zip_wd, &global_config.species_graph_info)?;        
    } else {
        fs::create_dir_all(&global_config.species_gfa)?;
        for species in all_species {
            let target = &global_config.gfa_build_wd.join(format!("{}.gfa", species));
            let dest = &global_config.species_gfa.join(target.file_name().unwrap());
            fs::rename(target, dest)?;
        }
    }
    fs::rename(&global_config.reference_pangenome_gfa, &global_config.reference_pangenome_gfa_db)?;
    fs::rename(&global_config.range_file, &global_config.range_file_db)?;
    Ok(())
}

fn move_dir_all(source: &Path, target: &Path) -> Result<()> {
    for entry in fs::read_dir(source)? {
        let entry = entry?;
        let dest = target.join(entry.file_name());
        
        if entry.path().is_dir() {
            fs::create_dir_all(&dest)?;
            move_dir_all(&entry.path(), &dest)?;
            fs::remove_dir(&entry.path())?;
        } else {
            fs::rename(entry.path(), dest)?;
        }
    }
    Ok(())
}

fn vg_combine(all_species: &Vec<&String>, gfa_build_wd: &PathBuf, reference_pangenome_vg: &PathBuf, vg: &PathBuf) -> Result<()> {
    // let mut command = Command::new(vg);
    // command.arg("combine").arg("-p");
    // for species in all_species {
    //     let vg_file = gfa_build_wd.join(format!("{}.vg", species));
    //     command.arg(&vg_file);
    // }
    // let output_file = File::create(reference_pangenome_vg)
    //     .with_context(|| format!("Failed to create output file: {:?}", reference_pangenome_vg))?;
    // command.stdout(output_file);
    // command.arg(">").arg(reference_pangenome_vg);
    // dbg!(&command);
    // let status = command
    //     .status()
    //     .with_context(|| "Failed to execute vg combine command")?;
    
    // if !status.success() {
    //     return Err(anyhow!("vg combine failed with exit code: {}", status));
    // }
    let mut vg_files: Vec<String> = Vec::new();
    for species in all_species {
        let vg_file = gfa_build_wd.join(format!("{}.vg", species));
        vg_files.push(vg_file.to_string_lossy().to_string());
    }
    let vg_files_str = vg_files.join(" ");   
    let cmd = format!(
        "{} combine -p {} > {}",
        vg.display(),
        vg_files_str,
        reference_pangenome_vg.display(),
    );
    run_shell_command(&cmd)?;
    if !reference_pangenome_vg.exists() {
        return Err(anyhow!("Output file was not created: {:?}", reference_pangenome_vg));
    }
    Ok(())
}

fn vg_convert_pan(vg: &PathBuf, threads: usize, reference_pangenome_vg: &PathBuf, reference_pangenome_gfa: &PathBuf) -> Result<()> {
    let cmd = format!(
        "{} convert -f {} -t {} > {}",
        vg.display(),
        reference_pangenome_vg.display(),
        threads,
        reference_pangenome_gfa.display(),
    );
    
    run_shell_command(&cmd)?;
    if !reference_pangenome_gfa.exists() {
        return Err(anyhow!("GFA file not created: {:?}", reference_pangenome_gfa));
    }
    Ok(())
}

fn _process_genome(genome: &String, file_path: &String) -> FxHashMap<String, String> {
    let mut seq_ids = FxHashMap::default();

    let reader = parse_fastx_file(Path::new(file_path));
    if let Err(e) = reader {
        eprintln!("Failed to open file {}: {}", file_path, e);
        return seq_ids;
    }

    let mut reader = reader.unwrap();
    while let Some(record) = reader.next() {
        match record {
            Ok(rec) => {
                let seq_id = String::from_utf8_lossy(rec.id()).split_whitespace().next().unwrap().to_string();
                seq_ids.insert(seq_id, genome.to_string());
            }
            Err(e) => eprintln!("Error reading record in {}: {}", file_path, e),
        }
    }

    seq_ids
}

fn prepare_species_genomes(genomes_metadata: &Vec<GenomesInfo>, pangeno_dir: &PathBuf, single_species_file: &PathBuf) -> Result<(FxHashMap<String, PathBuf>, FxHashMap<String, Vec<PathBuf>>)> {
    let mut species_groups: FxHashMap<String, Vec<&GenomesInfo>> = FxHashMap::default();
    for genome in genomes_metadata {
        species_groups
            .entry(genome.species_taxid.clone())
            .or_insert_with(Vec::new)
            .push(genome);
    }

    let (multi_species, single_species): (FxHashMap<_, _>, FxHashMap<_, _>) = species_groups
        .into_iter()
        .partition(|(_, genomes)| genomes.len() > 1);

    let mut multi_species2genomes: FxHashMap<String, Vec<PathBuf>> = FxHashMap::default();
    if !multi_species.is_empty() {
        for (species_taxid, genomes) in multi_species.iter() {
            let file_path = pangeno_dir.join(format!("{}.txt", species_taxid));
            let mut file = File::create(&file_path)?;
            
            for genome in genomes {
                let genome_path = PathBuf::from(&genome.path_id);
                if !genome_path.exists() {
                    eprintln!("{} does not exist, please check the path.", genome.path_id);
                    continue;
                }
                if !is_fasta(&genome_path) {
                    eprintln!("{} is not fasta format, skip it, please check the path.", genome.path_id);
                    continue;                    
                }
                multi_species2genomes.entry(species_taxid.clone()).or_insert_with(Vec::new).push(genome_path);
                writeln!(file, "{}", genome.path_id)?;
            }
        }
    }

    let mut single_species2genome: FxHashMap<String, PathBuf> = FxHashMap::default();
    if !single_species.is_empty() {
        let mut file = File::create(&single_species_file)?;
        
        for (species_taxid, genomes) in single_species.iter() {
            assert_eq!(genomes.len(), 1);
            if let Some(genome) = genomes.first() {
                let genome_path = PathBuf::from(&genome.path_id);
                if !genome_path.exists() {
                    eprintln!("{} does not exist, please check the path.", genome.path_id);
                    continue;
                }
                if !is_fasta(&genome_path) {
                    eprintln!("{} is not fasta format, skip it, please check the path.", genome.path_id);
                    continue;                    
                }
                single_species2genome.insert(species_taxid.clone(), genome_path);
                writeln!(file, "{}\t{}", species_taxid, genome.path_id)?;
            }
        }
    }

    Ok((single_species2genome, multi_species2genomes))
}

fn is_fasta<P: AsRef<Path>>(file: P) -> bool {
    let path = file.as_ref();

    if let Some(file_name) = path.to_str() {
        file_name.ends_with(".fa") ||
            file_name.ends_with(".fna") ||
            file_name.ends_with(".fasta") ||
            file_name.ends_with(".fa.gz") ||
            file_name.ends_with(".fna.gz") ||
            file_name.ends_with(".fasta.gz")
    } else {
        eprintln!("{:?} can not convert to utf-8 string", path);
        false
    }

}

fn _write_multi_species2genomes(multi_species2genomes: &FxHashMap<String, Vec<PathBuf>>, output_file_dir: &PathBuf) -> Result<()> {
    let _results: Vec<Result<()>> = multi_species2genomes
        .par_iter()
        .map(|(species_taxid, genome_paths)| {
            let output_file_path = output_file_dir.join(format!("{}.txt", species_taxid));
            if output_file_path.exists() {
                return Ok(());
            }
            let file = File::create(&output_file_path)?;
            let mut writer = BufWriter::new(file);
            
            for path in genome_paths {
                writeln!(writer, "{}", path.display())?;
            }
            
            writer.flush()?;
            Ok(())
        })
        .collect();
    Ok(())
}