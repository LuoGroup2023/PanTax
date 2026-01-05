
use std::fs;
use crate::cli::Cli;
use crate::utils::{run_shell_command, GlobalConfig};
use anyhow::{Result, anyhow};


pub fn index(args: &Cli, global_config: &GlobalConfig) -> Result<()> {
    let index_files = &global_config.index_files;
    let cmd1 = format!(
        "{} gbwt -g {} --gbz-format -G {} --num-jobs {}",
        &args.vg.display(),
        index_files[0].display(),
        &global_config.reference_pangenome_gfa_db.display(),
        args.threads,
    );
    
    run_shell_command(&cmd1)?;
    if !index_files[0].exists() {
        return Err(anyhow!("GFA file not created: {:?}", index_files[0]));
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
        return Err(anyhow!("GFA file not created: {:?}", index_files[1]));
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
        return Err(anyhow!("GFA file not created: {:?}", index_files[2]));
    }

    if args.index {
        if !args.debug {
            fs::remove_dir_all(&args.tmp_dir)?;
        }
        std::process::exit(0);
    }
    Ok(())
}