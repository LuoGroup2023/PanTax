
use std::fs;

use crate::cli::Cli;
use crate::utils::{run_shell_command, GlobalConfig};
use crate::types::DataType;
use crate::gaf_filter::filter_max_alignment_mt;
use anyhow::{Result, anyhow};

pub fn alignemnt(args: &Cli, global_config: &GlobalConfig, mode: &DataType) -> Result<()> {
    let index_files = &global_config.index_files;
    match mode {
        DataType::ShortReadSingle => {
            let cmd = format!(
                "{} giraffe -Z {} -m {} -d {} -f {} -t {} --named-coordinates -o gaf > {}",
                &args.vg.display(),
                index_files[0].display(),
                index_files[2].display(),
                index_files[1].display(),
                &args.read_files.as_ref().unwrap()[0].display(),
                args.threads,
                &global_config.gfa_mapped.display()
            );
            
            run_shell_command(&cmd)?;
            if !&global_config.gfa_mapped.exists() {
                return Err(anyhow!("GFA file not created: {:?}", &global_config.gfa_mapped));
            }
            Ok(())
        }
        DataType::ShortReadPaired => {
            let cmd = format!(
                "{} giraffe -Z {} -m {} -d {} -f {} -f {} -t {} --named-coordinates -o gaf > {}",
                &args.vg.display(),
                index_files[0].display(),
                index_files[2].display(),
                index_files[1].display(),
                &args.read_files.as_ref().unwrap()[0].display(),
                &args.read_files.as_ref().unwrap()[1].display(),
                args.threads,
                &global_config.gfa_mapped.display()
            );
            
            run_shell_command(&cmd)?;
            if !&global_config.gfa_mapped.exists() {
                return Err(anyhow!("GFA file not created: {:?}", &global_config.gfa_mapped));
            }
            Ok(())            
        }
        DataType::ShortReadPairedInter => {
            let cmd = format!(
                "{} giraffe -Z {} -m {} -d {} -i -f {} -t {} --named-coordinates -o gaf > {}",
                &args.vg.display(),
                index_files[0].display(),
                index_files[2].display(),
                index_files[1].display(),
                &args.read_files.as_ref().unwrap()[0].display(),
                args.threads,
                &global_config.gfa_mapped.display()
            );
            
            run_shell_command(&cmd)?;
            if !&global_config.gfa_mapped.exists() {
                return Err(anyhow!("GFA file not created: {:?}", &global_config.gfa_mapped));
            }
            Ok(())            
        }
        DataType::LongReadSingle => {
            
            let precise = match &args.long_read_type {
                Some(long_read_type) => {
                    match long_read_type.as_str().to_lowercase().as_str() {
                        "ontr9" | "clr" => 0.75,
                        "hifi" => 0.9,
                        "ontr10" => 0.8,
                        _ => args.precise_clipping
                    }
                }
                None => args.precise_clipping,
            };
            let cmd = format!(
                "{} -g {} -f {} -a {} -x vg -t {} --precise-clipping {}",
                args.lr_aligner.display(),
                &global_config.reference_pangenome_gfa_db.display(),
                &args.read_files.as_ref().unwrap()[0].display(),
                &global_config.gfa_mapped.display(),
                args.threads,
                precise
            );
            
            run_shell_command(&cmd)?;
            if !&global_config.gfa_mapped.exists() {
                return Err(anyhow!("GFA file not created: {:?}", &global_config.gfa_mapped));
            }

            let filter_gaf = filter_max_alignment_mt(&global_config.gfa_mapped)?;
            fs::rename(filter_gaf, &global_config.gfa_mapped)?;
            Ok(())            
        }
    }
}

