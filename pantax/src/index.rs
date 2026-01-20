
use std::fs;
use crate::cli::{Cli, VgLevel};
use crate::utils::{run_shell_command, GlobalConfig};
use anyhow::{Result, anyhow};


pub fn index(args: &Cli, global_config: &GlobalConfig) -> Result<()> {
    log::info!("Create index...");
    let index_files = &global_config.index_files;
    match args.vg_lvl {
        VgLevel::V1 => if args.auto {
            let cmd = format!(
                "{} autoindex -p {} -w giraffe -g {} -t {}",
                &args.vg.display(),
                &global_config.autoindex_prefix.display(),
                &global_config.reference_pangenome_gfa_db.display(),
                args.threads,
            );
            run_shell_command(&cmd)?;
            if !index_files[2].exists() {
                return Err(anyhow!(".min file not created: {:?}", index_files[2]));
            }
        } else {
            let cmd1 = format!(
                "{} gbwt -g {} --gbz-format -G {} --num-jobs {}",
                &args.vg.display(),
                index_files[0].display(),
                &global_config.reference_pangenome_gfa_db.display(),
                args.threads,
            );
            
            run_shell_command(&cmd1)?;
            if !index_files[0].exists() {
                return Err(anyhow!(".gbz index file not created: {:?}", index_files[0]));
            }

            let cmd2 = format!(
                "{} index -j {} {} -t {}",
                &args.vg.display(),
                index_files[1].display(),
                index_files[0].display(),
                args.threads,
            );
            
            run_shell_command(&cmd2)?;
            if !index_files[1].exists() {
                return Err(anyhow!(".dist file not created: {:?}", index_files[1]));
            }

            let cmd3 = format!(
                "{} minimizer -d {} -o {} {}",
                &args.vg.display(),
                index_files[1].display(),
                index_files[2].display(),
                index_files[0].display(),
            );
            
            run_shell_command(&cmd3)?;
            if !index_files[2].exists() {
                return Err(anyhow!(".min file not created: {:?}", index_files[2]));
            }
        }
        VgLevel::V2 => {
            enum Flag {
                Short,
                Long
            }
            let flag = if args.long_read {
                Flag::Long
            } else {
                log::warn!("Vg index does not specify reads type (short or long), default is for short.");
                Flag::Short
            };

            match flag {
                Flag::Short => {
                    if args.auto {
                        let cmd = format!(
                            "{} autoindex -p {} -w sr-giraffe -g {} -t {}",
                            &args.vg.display(),
                            &global_config.autoindex_prefix.display(),
                            &global_config.reference_pangenome_gfa_db.display(),
                            args.threads,
                        );
                        run_shell_command(&cmd)?;
                        if !index_files[2].exists() {
                            return Err(anyhow!(".min file not created: {:?}", index_files[2]));
                        }   
                    } else {
                        let cmd1 = format!(
                            "{} gbwt -g {} --gbz-format -G {} --num-jobs {}",
                            &args.vg.display(),
                            index_files[0].display(),
                            &global_config.reference_pangenome_gfa_db.display(),
                            args.threads,
                        );
                        
                        run_shell_command(&cmd1)?;
                        if !index_files[0].exists() {
                            return Err(anyhow!(".gbz index file not created: {:?}", index_files[0]));
                        }

                        let cmd2 = format!(
                            "{} index -j {} {} -t {}",
                            &args.vg.display(),
                            index_files[1].display(),
                            index_files[0].display(),
                            args.threads,
                        );
                        
                        run_shell_command(&cmd2)?;
                        if !index_files[1].exists() {
                            return Err(anyhow!(".dist file not created: {:?}", index_files[1]));
                        }

                        let cmd3 = format!(
                            "{} minimizer -d {} -o {} -Z {} {}",
                            &args.vg.display(),
                            index_files[1].display(),
                            index_files[2].display(),
                            index_files[3].display(),
                            index_files[0].display(),
                        );
                        
                        run_shell_command(&cmd3)?;
                        if !index_files[2].exists() {
                            return Err(anyhow!(".min file not created: {:?}", index_files[2]));
                        }
                                        
                    }
                }
                Flag::Long => {
                    log::info!("Long read index construction default use vg autoindex.");
                    let cmd = format!(
                        "{} autoindex -p {} -w lr-giraffe -g {} -t {}",
                        &args.vg.display(),
                        &global_config.autoindex_prefix.display(),
                        &global_config.reference_pangenome_gfa_db.display(),
                        args.threads,
                    );
                    run_shell_command(&cmd)?;
                    if !index_files[2].exists() {
                        return Err(anyhow!(".min file not created: {:?}", index_files[2]));
                    }                  
                }
            }
        
        } 
    } 
    log::info!("Create index...done");
    if args.index {
        if !args.debug {
            fs::remove_dir_all(&args.tmp_dir)?;
        }
        std::process::exit(0);
    }
    Ok(())
}