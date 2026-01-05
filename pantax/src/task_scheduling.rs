
use crate::cli::Cli;
use crate::fastixe::process_all_fasta_and_merge;
use crate::zip::{CompressType, zip};
use crate::utils::get_config;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use glob::glob;
use anyhow::{Context, Result, anyhow};
// use async_trait::async_trait;

use fxhash::{FxHashMap, FxHashSet};
use tokio::process::Command;
use tokio::sync::{mpsc, Semaphore, Mutex};
use tokio::io::AsyncWriteExt;
// use tokio::time::sleep;
use rust_htslib::faidx::build;

use futures::future::select_all;
// use std::time::Duration;

#[derive(Debug, Clone)]
pub struct TaskConfig {
    wd: PathBuf,
    zip_wd: PathBuf,
    multi_species2genomes: FxHashMap<String, Vec<PathBuf>>,
    vg: PathBuf,
    pangenome_building_exe: PathBuf,
    reference: Option<PathBuf>,
    parallel: bool,
    pub save: bool,
    range_file: PathBuf,
    pub lz4: bool,
    pub zstd: bool,
    threads: usize,
    pub debug: bool,
    verbose: bool,
}

fn initialize_and_get_task_config(args: &Cli, multi_species2genomes: &FxHashMap<String, Vec<PathBuf>>) -> Result<TaskConfig> {
    log::info!("Pangenome building first step: building pangenome for species with more than two genomes using {}", args.pangenome_building_exe.file_name().unwrap().to_string_lossy());
    let global_config = get_config();
    let exe_name = args.pangenome_building_exe
        .to_string_lossy()
        .to_lowercase();
    
    if !exe_name.contains("pggb") && !exe_name.contains("cactus-pangenome") {
        return Err(anyhow!("PanTax does not support {}", exe_name));
    }
    let wd = global_config.gfa_build_wd.clone();

    if args.force {
        if wd.exists() {
            log::info!("Choose to forcibly rebuild all pangenomes, deleting directory {:?}", wd);
            fs::remove_dir_all(&wd)?;
        }
        let zip_wd = PathBuf::from(format!("{}2", wd.display()));
        if zip_wd.exists() {
            fs::remove_dir_all(&zip_wd)?;
        }
        fs::create_dir_all(&wd)?;
        log::info!("Old {:?} has been deleted.", wd);
    } else {
        if !wd.exists() {
            log::debug!("Create {:?} for pangenome construction.", wd);
            fs::create_dir_all(&wd)?;
        }
    }

    // Prepare save directory if needed
    let zip_wd = global_config.zip_wd.clone();
    if args.is_save() {
        fs::create_dir_all(&zip_wd)?;
    }

    let parallel = if args.no_parallel { false } else { true };
    let save = args.save;
    // let force = args.force;
    let debug = args.debug;
    let verbose = args.verbose;
    
    // let sleep = Duration::from_secs(2);
    let range_file = global_config.range_file.clone();
    let tc = TaskConfig {
        wd,
        zip_wd,
        // pan_species: pangeno_dir,
        multi_species2genomes: multi_species2genomes.clone(),
        vg: args.vg.clone(),
        pangenome_building_exe: args.pangenome_building_exe.clone(),
        reference: args.reference.clone(),
        parallel,
        save,
        range_file,
        lz4: args.lz,
        zstd: args.zstd,
        threads: args.threads,
        debug,
        verbose,
    };
    Ok(tc)
}




struct TaskScheduler {
    config: TaskConfig,
    pan_species_files: Vec<PathBuf>,
    finished_species: Arc<Mutex<Vec<String>>>,
    species2reference: FxHashMap<String, PathBuf>,
    // tasks: Vec<Task>,
    task_groups: Vec<SpeciesTaskGroup>,
    core_allocator: Arc<Semaphore>,
    progress_tracker: Arc<Mutex<ProgressTracker>>,
}

#[derive(Debug, Clone)]
struct SpeciesResult {
    species: String,
    success: bool,
}

struct ProgressTracker {
    total_species: usize,
    completed_species: FxHashSet<String>,
    failed_species: FxHashSet<String>,
    running_species: FxHashSet<String>,
}

impl ProgressTracker {
    fn new() -> Self {
        Self {
            total_species: 0,
            completed_species: FxHashSet::default(),
            failed_species: FxHashSet::default(),
            running_species: FxHashSet::default(),
        }
    }

    fn set_total_species(&mut self, total_species: usize) {
        self.total_species = total_species;
    }
    
    fn add_running_species(&mut self, species: &str) {
        self.running_species.insert(species.to_string());
    }
    
    fn complete_task(&mut self, species: &str, success: bool) {
        self.running_species.remove(species);
        if success {
            self.completed_species.insert(species.to_string());
        } else {
            self.failed_species.insert(species.to_string());
        }
    }
    
    fn get_progress(&self) -> (usize, usize, usize, f64) {
        let completed = self.completed_species.len();
        let failed = self.failed_species.len();
        let running = self.running_species.len();
        let percentage = if self.total_species > 0 {
            (completed as f64 / self.total_species as f64) * 100.0
        } else {
            0.0
        };
        (completed, failed, running, percentage)
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct Task {
    task_type: TaskType,
    species: String,
    id: usize,
    cores: usize,
}

#[derive(Debug, Clone)]
struct SpeciesTaskGroup {
    species: String,
    tasks: Vec<TaskType>,
    total_cores: usize,
}

#[derive(Debug, Clone)]
enum TaskType {
    PangenomeBuild {
        tool: PangenomeTool,
        species_taxid: String,
        threads: usize,
        genomes_num: usize,
    },
    ConvertToVG {
        species_taxid: String,
        threads: usize,
    },
    ZipGFA {
        gfa_path: PathBuf,
        output_dir: PathBuf,
        range_file: PathBuf,
        threads: usize,
        debug: bool,
    },
    Cleanup {
        dir_to_clean: PathBuf,
    },
}

#[derive(Debug, Clone)]
enum PangenomeTool {
    Pggb,
    CactusPangenome { reference: Option<PathBuf> },
}

impl TaskScheduler {
    async fn new(config: TaskConfig) -> Result<Self> {
        let mut scheduler = Self {
            config: config.clone(),
            pan_species_files: Vec::new(),
            // finished_species: Vec::new(),
            finished_species: Arc::new(Mutex::new(Vec::new())),
            species2reference: FxHashMap::default(),
            // tasks: Vec::new(),
            task_groups: Vec::new(),
            core_allocator: Arc::new(Semaphore::new(config.threads)),
            // executor: Arc::new(ShellCommandExecutor),
            progress_tracker: Arc::new(Mutex::new(ProgressTracker::new())),
        };
        
        scheduler.initialize().await?;
        Ok(scheduler)
    }

    async fn initialize(&mut self) -> Result<()> {        
        // Load finished pangenomes
        let finished_file = self.config.wd.join("finished_pangenome.txt");
        if finished_file.exists() {
            let file = File::open(&finished_file)?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(species) = line {
                    // self.finished_species.push(species);
                    self.finished_species.lock().await.push(species);
                }
            }
        } else {
            File::create(&finished_file)?;
        }
        
        log::info!("Already finished species: {}", self.finished_species.lock().await.len());
        
        // Load reference genomes if provided
        if let Some(ref reference_path) = self.config.reference {
            let file = File::open(reference_path)?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line?;
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 2 {
                    self.species2reference.insert(
                        parts[0].to_string(),
                        PathBuf::from(parts[1]),
                    );
                }
            }
        }
        
        // Calculate fold factor for thread allocation
        let mut used_threads = 0;
        let mut fold = 1.0;
        
        for (_species_taxid, genome_paths) in &self.config.multi_species2genomes {
            let threads = self.get_max_cores_based_on_genomes(genome_paths.len(), 8, fold);
            used_threads += threads;
            if used_threads > self.config.threads {
                break;
            }
        }
        
        if used_threads < self.config.threads {
            fold = self.config.threads as f64 / used_threads as f64;
        }
        
        self.generate_task_groups(fold).await?;
        
        let mut tracker = self.progress_tracker.lock().await;
        *tracker = ProgressTracker::new();
        tracker.set_total_species(self.task_groups.len());        
        Ok(())
    }

    fn get_max_cores_based_on_genomes(&self, num_genomes: usize, base_threads: usize, fold: f64) -> usize {
        let threads = match num_genomes {
            n if n < 5 => base_threads,
            n if (5..=10).contains(&n) => base_threads * 2,
            n if (11..=25).contains(&n) => base_threads * 4,
            n if (26..=50).contains(&n) => base_threads * 6,
            n if (51..=100).contains(&n) => base_threads * 8,
            _ => 64,
        };
        
        let adjusted = (threads as f64 * fold).ceil() as usize;
        adjusted.min(self.config.threads)
    }

    async fn generate_task_groups(&mut self, fold: f64) -> Result<()> {
        let exe_name = self.config.pangenome_building_exe
            .to_string_lossy()
            .to_lowercase();
        let zip_wd = &self.config.zip_wd;
        for (species_taxid, genome_paths) in &self.config.multi_species2genomes {
            if self.finished_species.lock().await.contains(species_taxid) {
                continue;
            }
            let genomes_num = genome_paths.len();
            if genomes_num < 2 {
                log::debug!("Skipping species {} with only {} genome(s)", species_taxid, genomes_num);
                continue;
            }
            let threads = self.get_max_cores_based_on_genomes(genomes_num, 8, fold);
            let pangenome_tool = if exe_name.contains("pggb") {
                PangenomeTool::Pggb
            } else if exe_name.contains("cactus-pangenome") {
                let reference = self.species2reference.get(species_taxid.as_str()).cloned();
                PangenomeTool::CactusPangenome { reference }
            } else {
                return Err(anyhow!("Unsupported pangenome tool: {}", exe_name));
            };  
            
            let mut tasks = Vec::new();
            tasks.push(TaskType::PangenomeBuild {
                tool: pangenome_tool,
                species_taxid: species_taxid.clone(),
                threads,
                genomes_num,
            });
            tasks.push(TaskType::ConvertToVG {
                species_taxid: species_taxid.clone(),
                threads,
            });
            if self.config.save || self.config.lz4 || self.config.zstd {
                tasks.push(TaskType::ZipGFA {
                    gfa_path: self.config.wd.join(format!("{}.gfa", species_taxid)),
                    output_dir: zip_wd.clone(),
                    range_file: self.config.range_file.clone(),
                    threads,
                    debug: self.config.debug,
                });
            }
            if !self.config.debug {
                tasks.push(TaskType::Cleanup {
                    dir_to_clean: self.config.wd.join(species_taxid),
                });
            }
            let task_group = SpeciesTaskGroup {
                species: species_taxid.clone(),
                tasks,
                total_cores: threads.max(1), 
            };
            self.task_groups.push(task_group);
        } 
        log::info!("Generated {} task groups to execute", self.task_groups.len());
        Ok(())
    }
    
    async fn execute_task_direct(&self, task: Task) -> Result<String> {
        match task.task_type {
            TaskType::PangenomeBuild { tool, species_taxid, threads, genomes_num } => {
                self.execute_pangenome_build(tool, species_taxid, threads, genomes_num).await
            }
            TaskType::ConvertToVG { species_taxid, threads } => {
                self.execute_convert_to_vg(species_taxid, threads).await
            }
            TaskType::ZipGFA { gfa_path, output_dir, range_file, threads, debug } => {
                self.execute_zip_gfa(gfa_path, output_dir, range_file, threads, debug).await
            }
            TaskType::Cleanup { dir_to_clean } => {
                self.execute_cleanup(dir_to_clean).await
            }
        }
    }

    async fn execute_pangenome_build(
        &self,
        tool: PangenomeTool,
        species_taxid: String,
        threads: usize,
        genomes_num: usize,
    ) -> Result<String> {
        match tool {
            PangenomeTool::Pggb => {
                self.execute_pggb(species_taxid, threads, genomes_num).await
            }
            PangenomeTool::CactusPangenome { reference } => {
                self.execute_cactus_pangenome(species_taxid, threads).await
            }
        }
    }

    async fn execute_pggb(&self, species_taxid: String, threads: usize, genomes_num: usize) -> Result<String> {
        let species_wd = self.config.wd.join(&species_taxid);
        let time_log = if self.config.verbose {
            format!("/usr/bin/time -v -o {}/{}_pangenome_building_time.log ", 
                    species_wd.display(), species_taxid)
        } else {
            String::new()
        };
        
        let pggb_cmd = format!(
            "{}{} -i {}/{}_merged.fa.gz -o {}/species{}_pangenome_building -t {} -p 90 -n {} -v > /dev/null 2>&1",
            time_log,
            self.config.pangenome_building_exe.display(),
            species_wd.display(),
            species_taxid,
            species_wd.display(),
            species_taxid,
            threads,
            genomes_num
        );
        
        self.run_shell_command(&pggb_cmd).await?;
        Ok("pggb".to_string())
    }

    async fn execute_cactus_pangenome(&self, species_taxid: String, threads: usize) -> Result<String> {
        let species_wd = self.config.wd.join(&species_taxid);
        let genome2id = species_wd.join("genome2id.tsv");
        let pan_species_file = self.config.wd.join(format!("species_pangenome/{}.txt", species_taxid));
        if !pan_species_file.exists() && !self.config.debug {
            eprintln!("This function is experimental, please run it with --debug.");
            std::process::exit(1);
        }
        let mut commands = Vec::new();
        // Create genome2id.tsv
        let genome2id_cmd = format!(
            "awk -F'/' '{{file=$NF; sub(/\\..*/, \"\", file); print file \"\\t\" $0}}' {} > {}",
            pan_species_file.display(),
            genome2id.display()
        );
        commands.push(genome2id_cmd);
        
        let time_log = if self.config.verbose {
            format!("/usr/bin/time -v -o {}/{}_pangenome_building_time.log ", 
                    species_wd.display(), species_taxid)
        } else {
            String::new()
        };  
        if let Some(reference) = self.species2reference.get(species_taxid.as_str()) {
            let cactus_cmd = format!(
                "{} {} {}/js {} --outDir {}/species{}_pangenome_building --outName {} --reference {} --maxCores {} --mapCores {} > {}/{}_mc.log 2>&1",
                time_log,
                self.config.pangenome_building_exe.display(),
                species_wd.display(),
                genome2id.display(),
                species_wd.display(),
                species_taxid,
                species_taxid,
                reference.display(),
                threads,
                std::cmp::max(threads / 4, 6),
                species_wd.display(),
                species_taxid
            );
            commands.push(cactus_cmd);
        } else {
            commands.push(format!("reference=$(awk 'NR==1 {{print $1}}' {})", genome2id.display()));
            
            let cactus_cmd = format!(
                "{} {} {}/js {} --outDir {}/species{}_pangenome_building --outName {} --reference $reference --maxCores {} --mapCores {} > {}/{}_mc.log 2>&1",
                time_log,
                self.config.pangenome_building_exe.display(),
                species_wd.display(),
                genome2id.display(),
                species_wd.display(),
                species_taxid,
                species_taxid,
                threads,
                std::cmp::max(threads / 4, 6),
                species_wd.display(),
                species_taxid
            );
            commands.push(cactus_cmd);
            
            let gunzip_cmd = format!(
                "gunzip {}/species{}_pangenome_building/{}.gfa.gz",
                species_wd.display(),
                species_taxid,
                species_taxid
            );
            commands.push(gunzip_cmd);
        }   
        let full_command = commands.join("; ");
        self.run_shell_command(&full_command).await?;
        Ok("mc".to_string())
        
    }

    async fn execute_convert_to_vg(&self, species_taxid: String, threads: usize) -> Result<String> {
        let gfa_pattern = self.config.wd.join(format!("{}/species{}_pangenome_building/*.gfa", species_taxid, species_taxid));
        let vg_path = self.config.wd.join(format!("{}.vg", species_taxid));
        let cmd = format!(
            "{} convert -g {} -p -t {} > {}",
            self.config.vg.display(),
            gfa_pattern.display(),
            threads,
            vg_path.display()
        );
        
        self.run_shell_command(&cmd).await?;
        if !vg_path.exists() {
            return Err(anyhow!("VG file not created: {:?}", vg_path));
        }
        self.move_gfa_files(&self.config.wd, &species_taxid)?;
        Ok("convert_to_vg".to_string())
    }

    async fn execute_zip_gfa(&self, gfa_path: PathBuf, output_dir: PathBuf, range_file: PathBuf, threads: usize, _debug: bool) -> Result<String> {
        // dbg!(&self.config);
        let compress_type = CompressType::from_config(&self.config)?; 
        zip(&gfa_path, &output_dir, &range_file, threads, &compress_type)?;
        fs::remove_file(gfa_path)?;
        Ok("zip".to_string())
    }

    async fn run_shell_command(&self, command: &str) -> Result<String> {        
        let output = Command::new("bash")
            .arg("-c")
            .arg(command)
            .output()
            .await
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
        let stdout = String::from_utf8_lossy(&output.stdout);
        
        if !stderr.trim().is_empty() {
            log::warn!("Command succeeded with stderr output:\n{}", stderr);
        }
        
        if self.config.verbose && !stdout.trim().is_empty() {
            log::debug!("Command stdout: {}", stdout);
        }
        
        Ok(stdout.trim().to_string())
    }

    async fn execute_cleanup(&self, dir_to_clean: PathBuf) -> Result<String> {        
        if dir_to_clean.exists() {
            tokio::fs::remove_dir_all(&dir_to_clean).await?;
        }
        
        Ok("cleanup".to_string())
    }


    fn move_gfa_files(&self, wd: &Path, species_taxid: &String) -> Result<()> {
        let source_pattern = format!("{}/{}/species{}_pangenome_building/*.gfa", 
            wd.display(), species_taxid, species_taxid);
    
        let dest_path = wd.join(format!("{}.gfa", species_taxid));
        
        let files: Vec<PathBuf> = glob(&source_pattern)?
            .filter_map(Result::ok)
            .collect();
        
        if files.is_empty() {
            anyhow::bail!("No .gfa files found matching pattern: {}", source_pattern);
        }
        
        if files.len() == 1 {
            let src = &files[0];
            fs::rename(src, &dest_path)?;
            // println!("Moved {} to {}", src.display(), dest_path.display());
        } else {
            eprintln!("It is not possible to get more than one gfa file.");
            std::process::exit(1);
        }
        
        Ok(())
    }

    async fn execute_task_group(&mut self, group: SpeciesTaskGroup) -> Result<bool> {
        log::debug!("Starting task group for species: {}", group.species);
        
        let permit = self.core_allocator.acquire_many(group.total_cores as u32).await?;
        
        {
            let mut tracker = self.progress_tracker.lock().await;
            tracker.add_running_species(&group.species);
        }
        
        let mut all_success = true;
        
        for task_type in group.tasks {
            let task = Task {
                task_type: task_type.clone(),
                species: group.species.clone(),
                id: 0, 
                cores: match &task_type {
                    TaskType::PangenomeBuild { threads, .. } => *threads,
                    _ => 1,
                },
            };
            
            // log::info!("Executing task for species {}: {:?}", group.species, task_type);
            
            match self.execute_task_direct(task).await {
                Ok(_output) => {
                    // log::debug!("Task completed successfully: {}", output);
                    
                    // If the critical task fails, exit early.
                    if matches!(task_type, TaskType::PangenomeBuild { .. }) {
                        // check gfa
                        if let Err(e) = self.check_gfa_generation(&group.species).await {
                            log::error!("GFA generation check failed for {}: {}", group.species, e);
                            all_success = false;
                            break;
                        }
                    }
                }
                Err(e) => {
                    log::error!("Task failed for species {}: {}", group.species, e);
                    all_success = false;
                    
                    // If the critical task fails, exit early.
                    if matches!(task_type, TaskType::PangenomeBuild { .. }) {
                        break;
                    }
                    // Non-critical task failures allow continuation, but are marked as failed.
                }
            }
        }
        
        drop(permit);
        
        // last check
        if all_success {
            all_success = self.check_species_completion(&group.species).await;
        }
        
        log::debug!("Task group for species {} completed: {}", 
                group.species, if all_success { "success" } else { "failed" });
        
        Ok(all_success)
    }

    async fn check_gfa_generation(&self, species: &str) -> Result<()> {
        let species_wd = self.config.wd.join(species);
        let gfa_pattern = format!(
            "{}/species{}_pangenome_building/*.gfa",
            species_wd.display(), species
        );
        
        let gfa_files: Vec<PathBuf> = glob(&gfa_pattern)?
            .filter_map(Result::ok)
            .collect();
        
        if gfa_files.is_empty() {
            anyhow::bail!("No GFA files generated by pggb for species {}", species);
        }
        
        // log::debug!("Found {} GFA files for species {}", gfa_files.len(), species);
        Ok(())
    }

    async fn check_species_completion(&mut self, species: &str) -> bool {
        let gfa_path = self.config.wd.join(format!("{}.gfa", species));
        let vg_path = self.config.wd.join(format!("{}.vg", species));
        
        let compressed_exists = if self.config.save {
            let ext = if self.config.lz4 {
                "bin.lz4"
            } else if self.config.zstd {
                "bin.zst"
            } else {
                "bin"
            };
            let compressed_path = PathBuf::from(format!("{}2/{}.{}", self.config.wd.display(), species, ext));
            compressed_path.exists()
        } else {
            true
        };
        
        let gfa_exists = gfa_path.exists() || compressed_exists;
        let vg_exists = vg_path.exists();
        
        log::debug!("Checking completion for species {}: GFA exists={}, VG exists={}", 
                species, gfa_exists, vg_exists);
        
        if gfa_exists && vg_exists {
            self.finished_species.lock().await.push(species.to_string());
            log::info!("Species {} completed successfully!", species);
            true
        } else {
            // log::warn!("Species {} incomplete: GFA={}, VG={}", 
            //         species, gfa_exists, vg_exists);
            false
        }
    }

    async fn _run_parallel(&self) -> Result<()> {
        let total_species = self.task_groups.len();
        if total_species == 0 {
            log::info!("No species to process");
            return Ok(());
        }
        
        log::info!("Starting parallel execution with {} species", total_species);
        
        let (result_tx, mut result_rx) = mpsc::channel::<SpeciesResult>(100);
        let progress_tracker = Arc::clone(&self.progress_tracker);
        
        let result_processor_task: tokio::task::JoinHandle<anyhow::Result<(usize, usize)>> = tokio::spawn({
            let progress_tracker = progress_tracker.clone();
            let finished_file = self.config.wd.join("finished_pangenome.txt");
            
            async move {
                let mut successes = 0;
                let mut failures = 0;
                
                let mut file = tokio::fs::OpenOptions::new()
                    .append(true)
                    .open(&finished_file)
                    .await?;
                
                while let Some(result) = result_rx.recv().await {
                    {
                        let mut tracker = progress_tracker.lock().await;
                        tracker.complete_task(&result.species, result.success);
                    }
                    
                    if result.success {
                        successes += 1;
                        file.write_all(format!("{}\n", result.species).as_bytes()).await?;
                        file.flush().await?; 
                        // log::debug!("Wrote species {} to finished file", result.species);
                    } else {
                        failures += 1;
                        log::warn!("Species {} failed", result.species);
                    }
                }
                
                Ok((successes, failures))
            }
        });
        
        let progress_monitor_task = {
            let progress_tracker = Arc::clone(&self.progress_tracker);
            tokio::spawn(async move {
                monitor_progress_smart(progress_tracker, total_species).await;
            })
        };
        
        let task_group_futures: Vec<_> = self.task_groups
            .iter()
            .cloned()
            .map(|task_group| {
                let result_tx = result_tx.clone();
                let mut scheduler = self.clone_for_task();
                let progress_tracker = progress_tracker.clone();
                
                tokio::spawn(async move {
                    {
                        let mut tracker = progress_tracker.lock().await;
                        tracker.add_running_species(&task_group.species);
                    }
                    
                    let success = scheduler.execute_task_group(task_group.clone()).await?;
                    
                    let result = SpeciesResult {
                        species: task_group.species,
                        success,
                    };
                    
                    result_tx.send(result).await?;
                    Ok::<(), anyhow::Error>(())
                })
            })
            .collect();
        
        let task_group_results = futures::future::join_all(task_group_futures).await;
        
        for (i, result) in task_group_results.iter().enumerate() {
            if let Err(e) = result {
                log::error!("Task group {} failed with error: {}", i, e);
            }
        }
        
        drop(result_tx);
        
        progress_monitor_task.await?;
        
        let (successes, failures) = result_processor_task.await??;
        
        log::info!("{}", "=".repeat(60));
        log::info!("PANGENOME CONSTRUCTION SUMMARY");
        log::info!("{}", "=".repeat(60));
        log::info!("Total species: {}", total_species);
        log::info!("Successfully completed: {} ({:.1}%)", 
                successes, (successes as f64 / total_species as f64) * 100.0);
        log::info!("Failed: {} ({:.1}%)", 
                failures, (failures as f64 / total_species as f64) * 100.0);
        log::info!("{}", "=".repeat(60));
        
        if failures > 0 {
            log::warn!("{} species failed during execution", failures);
        }
        
        Ok(())
    }

    async fn run_parallel_smart(&self) -> Result<()> {
        let total_species = self.task_groups.len();
        if total_species == 0 {
            log::info!("No species to process");
            return Ok(());
        }
        
        log::info!("Starting smart parallel execution with {} cores, {} species", 
                self.config.threads, total_species);
        
        let (result_tx, mut result_rx) = mpsc::channel::<SpeciesResult>(100);
        let progress_tracker = Arc::clone(&self.progress_tracker);

        let result_processor_task: tokio::task::JoinHandle<anyhow::Result<(usize, usize)>> = tokio::spawn({
            let progress_tracker = progress_tracker.clone();
            let finished_file = self.config.wd.join("finished_pangenome.txt");
            
            async move {
                let mut successes = 0;
                let mut failures = 0;
                
                let mut file = tokio::fs::OpenOptions::new()
                    .append(true)
                    .open(&finished_file)
                    .await?;
                
                while let Some(result) = result_rx.recv().await {
                    {
                        let mut tracker = progress_tracker.lock().await;
                        tracker.complete_task(&result.species, result.success);
                    }
                    
                    if result.success {
                        successes += 1;
                        file.write_all(format!("{}\n", result.species).as_bytes()).await?;
                        file.flush().await?; 
                        // self.finished_species.push(result.species.clone());
                        log::debug!("Wrote species {} to finished file", result.species);
                    } else {
                        failures += 1;
                        log::warn!("Species {} failed", result.species);
                    }
                }
                
                Ok((successes, failures))
            }
        });

        let progress_monitor_task = {
            let progress_tracker = Arc::clone(&self.progress_tracker);
            
            tokio::spawn(async move {
                monitor_progress_smart(progress_tracker, total_species).await;
            })
        };
        
        let mut task_queue: Vec<SpeciesTaskGroup> = self.task_groups.clone();
        task_queue.sort_by(|a, b| b.total_cores.cmp(&a.total_cores));
        
        let mut running_tasks = Vec::new();
        let mut running_species = Vec::new();
        let mut running_used_cores = Vec::new();
        let mut available_cores = self.core_allocator.available_permits();
        while !task_queue.is_empty() || !running_tasks.is_empty() {
            while !task_queue.is_empty() {
                let task_group = task_queue.last().unwrap();
                let num_cores = task_group.total_cores;
                if num_cores <= available_cores {
                    let task_group = task_queue.pop().unwrap();
                    log::debug!("Submitted task {} with {} cores", task_group.species, task_group.total_cores);
                    let species = task_group.species.clone();
                    let task = self.submit_single_task(
                        task_group, 
                        result_tx.clone(), 
                        progress_tracker.clone()
                    ).await?;
                    available_cores -= num_cores;
                    log::debug!("available_cores: {}", available_cores);
                    running_tasks.push(task);
                    running_species.push(species);
                    running_used_cores.push(num_cores);
                } else {
                    break;
                }
            }
            
            if !running_tasks.is_empty() && self.config.debug {
                log::debug!("{:?} species running...", running_species);
            }
            
            if !running_tasks.is_empty() {
                let (_, remove_idx, remaining_futures) = select_all(running_tasks).await;
                running_tasks = remaining_futures;
                running_species.remove(remove_idx);
                available_cores += running_used_cores.remove(remove_idx);
            }            
        }
        
        drop(result_tx);
        
        progress_monitor_task.await?;
        
        let (successes, failures) = result_processor_task.await??;
        
        log::info!("{}", "=".repeat(60));
        log::info!("EXECUTION COMPLETE - FINAL SUMMARY");
        log::info!("{}", "=".repeat(60));
        log::info!("Total species: {}", total_species);
        log::info!("Successfully completed: {} ({:.1}%)", 
                successes, (successes as f64 / total_species as f64) * 100.0);
        log::info!("Failed: {} ({:.1}%)", 
                failures, (failures as f64 / total_species as f64) * 100.0);
        log::info!("{}", "=".repeat(60));
        
        Ok(())
    }

    async fn submit_single_task(
        &self,
        task_group: SpeciesTaskGroup,
        result_tx: mpsc::Sender<SpeciesResult>,
        progress_tracker: Arc<Mutex<ProgressTracker>>,
    ) -> Result<tokio::task::JoinHandle<anyhow::Result<()>>, anyhow::Error> {
        let mut scheduler = self.clone_for_task();
        let species = task_group.species.clone();
        
        let handle = tokio::spawn(async move {
            // Mark species as running
            {
                let mut tracker = progress_tracker.lock().await;
                tracker.add_running_species(&species);
            }
            
            // Execute the task group
            let success = scheduler.execute_task_group(task_group).await?;
            
            // Send result
            result_tx.send(SpeciesResult {
                species,
                success,
            }).await?;
            
            Ok(())
        });
        
        Ok(handle)
    }

    async fn run_serial(&self) -> Result<()> {
        Ok(())
    }

    fn clone_for_task(&self) -> Self {
        Self {
            config: self.config.clone(),
            pan_species_files: self.pan_species_files.clone(),
            finished_species: self.finished_species.clone(),
            species2reference: self.species2reference.clone(),
            task_groups: Vec::new(),
            core_allocator: Arc::clone(&self.core_allocator),
            progress_tracker: Arc::clone(&self.progress_tracker),
        }
    }
    
    async fn run(&self) -> Result<()> {
        if self.task_groups.is_empty() {
            log::warn!("No task groups to execute");
            return Ok(());
        }
        
        log::info!("Starting execution with {} species", self.task_groups.len());
        
        if self.config.parallel {
            // self.run_parallel().await
            self.run_parallel_smart().await
        } else {
            self.run_serial().await
        }
    }
}

impl Clone for TaskScheduler {
    fn clone(&self) -> Self {
        Self {
            config: self.config.clone(),
            pan_species_files: self.pan_species_files.clone(),
            finished_species: self.finished_species.clone(),
            species2reference: self.species2reference.clone(),
            task_groups: self.task_groups.clone(),
            core_allocator: Arc::clone(&self.core_allocator),
            progress_tracker: Arc::clone(&self.progress_tracker),
        }
    }
}

fn prepare_genomes_before_pggb(config: &TaskConfig) -> Result<()> {
    for (species_taxid, genome_paths) in config.multi_species2genomes.iter() {
        // bgzip output. file path: wd/species/species_merged.fa.gz
        let output_file_dir = config.wd.join(species_taxid);
        if !output_file_dir.exists() {
            fs::create_dir_all(&output_file_dir)?;
        }
        let output_file_path = config.wd.join(format!("{}/{}_merged.fa.gz", species_taxid, species_taxid));
        process_all_fasta_and_merge(genome_paths, &output_file_path)?;
        build(&output_file_path).expect(&format!("Failed to build {} FASTA index", species_taxid.as_str()));
    }
    Ok(())
} 

async fn monitor_progress_smart(
    progress_tracker: Arc<Mutex<ProgressTracker>>,
    total_species: usize,
) {
    let update_threshold = 5.0; // update 5%
    let mut last_reported_percentage = 0.;
    log::debug!("Progress monitoring started for {} species", total_species);
    loop {
        let (completed, failed, running, current_percentage) = {
            let tracker = progress_tracker.lock().await;
            tracker.get_progress()
        };
        
        let should_report = 
            // // first report
            // (last_reported_percentage < 0.0 && completed > 0) ||
            // up to threshold
            (current_percentage - last_reported_percentage >= update_threshold) ||
            // all task finished
            (completed + failed >= total_species);

            // (total_species > 20 && completed - last_completed_count >= 10);
        // log::info!("completed {}, should_report {}", completed, should_report);
        if should_report {
            let success_rate = if completed + failed > 0 {
                (completed as f64 / (completed + failed) as f64) * 100.0
            } else {
                0.0
            };
            
            log::info!("Progress: {:.1}% | Success: {}/{} | Running: {} | Success rate: {:.1}%", 
                current_percentage, completed, completed + failed, running, success_rate);
            
            last_reported_percentage = current_percentage;
            
            if completed + failed >= total_species {
                if current_percentage < 100.0 {
                    log::info!("Progress: 100.0% | All {} species processed", total_species);
                }
                break;
            }
        }
        tokio::time::sleep(std::time::Duration::from_secs(2)).await;
    }
}


#[tokio::main]
pub async fn task_schedule_multi(args: &Cli, multi_species2genomes: &FxHashMap<String, Vec<PathBuf>>) -> Result<Vec<String>> {
    let config = initialize_and_get_task_config(args, multi_species2genomes)?;
    prepare_genomes_before_pggb(&config)?;
    let scheduler = TaskScheduler::new(config).await?;
    scheduler.run().await?;
    let finished_species = scheduler.finished_species.lock().await.clone();
    Ok(finished_species)
}