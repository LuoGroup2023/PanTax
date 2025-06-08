#!/usr/bin/env python3
import concurrent.futures
import subprocess, argparse, sys
import time, asyncio
from toolkits import Logger, delete_directory_contents, is_file_non_empty
from pathlib import Path
import pandas as pd

script_path = Path(__file__).resolve()
script_dir = script_path.parent 
usage = "Task scheduling strategy for pangenome construction."
"""
PGGB building genome time test: 8 threads-10 min, 16 threads-5 min, 32 threads-5 min, 64 threads-5 min for 10 genomes
Plan: < 5 genomes-4 threads, 5-10 genomes-8 threads, 10-25 genomes-16 threads, 25-50 genomes-32 threads, 50-100 genomes-64 threads
"""


class TaskScheduling:

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.support_pangenome_building_exe = ["pggb", "cactus-pangenome"]
        self.set_para()
        if self.reference:
            self.get_reference_genome()
        self.get_all_cmd()

    def set_para(self):   
        if self.debug.lower() == "true" or self.debug.lower() == "on":
            self.debug = True
            self.is_show_detailed_log = "true"
        else:
            self.debug = False

        is_show_detailed_log = self.is_show_detailed_log.lower()
        if is_show_detailed_log == "true" or is_show_detailed_log == "on":
            self.log = Logger(level="debug").logger
            self.is_show_detailed_log = True
        elif is_show_detailed_log == "false" or is_show_detailed_log == "off":
            self.log = Logger(level="info").logger
            self.is_show_detailed_log = False
        else:
            raise ValueError(f"The option {self.is_show_detailed_log} is not valid.")        

        exe = Path(self.pangenome_building_exe).stem.lower()
        if not any(exe in _exe for _exe in self.support_pangenome_building_exe):
                raise ValueError(f"PanTax does not support {Path(self.pangenome_building_exe).stem}")

        self.log.info(f"Pangenome building first step: building pangenome for species with more than two genomes using {Path(self.pangenome_building_exe).name}.")
        wd = Path(self.wd)

        force = self.force.lower()
        if force == "true" or force == "on":
            if wd.is_dir():
                self.log.info(f"Choose to forcibly rebuild all pangenomes, deleting directory {self.wd}.")        
                delete_directory_contents(self.wd)
                self.log.info(f"old {self.wd} has been deleted.")
            wd.mkdir(exist_ok=True)
        elif force == "false" or force == "off":
            if wd.exists():
                self.log.info("Resume to build pangenome.")
            else:
                self.log.debug(f"Create {self.wd} for pangenome construction.")
                wd.mkdir()
        else:
            raise ValueError(f"The option {self.force} is not valid.")

        # check pangenome need to be built
        pan_species = Path(self.pan_species)
        self.pan_species_files = list(pan_species.glob("*.txt"))
        self.pangenomes_need_to_build = len(self.pan_species_files)
        self.log.info(f"Species pangenome need to build: {self.pangenomes_need_to_build}.")

        # check finished pangenomes
        self.finished_pangenome = wd / "finished_pangenome.txt"
        if self.finished_pangenome.exists():
            finished_pangenome_df = pd.read_csv(self.finished_pangenome, sep="\t", header=None, dtype=object)
            self.finished_species = finished_pangenome_df.iloc[:,0].tolist()
            self.finished_pangenome_num = len(self.finished_species)
        else:
            self.finished_pangenome.touch()
            self.finished_species = []
            self.finished_pangenome_num = 0

        used_threads = 0
        self.fold = 1
        for pan_species_file in self.pan_species_files:    
            with open(pan_species_file, "r") as f:
                genomes = [line.strip() for line in f]
            genomes_num = len(genomes)
            threads = self.get_max_cores_based_on_genomes(genomes_num)
            used_threads += threads
            if used_threads > self.threads:
                break
        if used_threads < self.threads:
            self.fold = self.threads / used_threads

        if self.save.lower() == "true":
            self.save = True
        else:
            self.save = False

        if not self.sleep:
            self.sleep = 30

    def get_max_cores_based_on_genomes(self, num_genomes, base_threads=8):
        """Return max cores based on the number of genomes."""
        if num_genomes < 5:
            return min(self.threads, base_threads*self.fold)
        elif 5 <= num_genomes <= 10:
            return min(self.threads, base_threads*2*self.fold)
        elif 10 < num_genomes <= 25:
            return min(self.threads, base_threads*4*self.fold)
        elif 25 < num_genomes <= 50:
            return min(self.threads, base_threads*6*self.fold)
        elif 50 < num_genomes <= 100:
            return min(self.threads, base_threads*8*self.fold)
        else:
            # Default to 64 cores if over 100 genomes
            return self.threads

    def get_reference_genome(self):
        reference_genomes_df = pd.read_csv(self.reference, sep="\t", header=None, dtype=object)
        self.species2reference = dict(zip(reference_genomes_df.iloc[:, 0], reference_genomes_df.iloc[:, 1]))

    def get_all_cmd(self):
        # build pangenome
        self.all_cmd = []
        for pan_species_file in self.pan_species_files:
            cmd = []
            species = Path(pan_species_file).stem
            if species in self.finished_species: continue
            species_wd_dir = f"{self.wd}/{species}"
            if Path(species_wd_dir).exists():
                cmd1 = f"rm -rf {self.wd}/{species}; mkdir {self.wd}/{species}"
            else:
                cmd1 = f"mkdir {self.wd}/{species}"
            cmd.append(cmd1)
            with open(pan_species_file, "r") as f:
                genomes = [line.strip() for line in f]
            # print(pan_species_file)
            # print(genomes)
            genomes_num = len(genomes)
            threads = int(self.get_max_cores_based_on_genomes(genomes_num))
            if "pggb" in self.pangenome_building_exe.lower():
                cmd_tmp = f"{self.pantaxr} fastixe -l {pan_species_file} -b -o {self.wd}/{species} -e {species}_merged.fa --up"
                # print(cmd_tmp)
                cmd.append(cmd_tmp)
                # cmd_tmp_list = []
                # new_genomes_path = []
                # for genome in genomes:
                #     genome_name = Path(genome).name
                #     genome_name = "_".join(genome_name.split("_")[:2])
                #     if genome.endswith("gz"):
                #         gunzip_genome_name = Path(genome).name.replace(".gz", "")
                #         cmd_tmp_list.append(f"gunzip -c {genome} > {self.wd}/{species}/{gunzip_genome_name}; {self.fastix} {self.wd}/{species}/{gunzip_genome_name} -p '{genome_name}#1#' > {self.wd}/{species}/{genome_name}.fa")
                #     else:
                #         cmd_tmp_list.append(f"{self.fastix} {genome} -p '{genome_name}#1#' > {self.wd}/{species}/{genome_name}.fa") 
                #     new_genomes_path.append(f"{self.wd}/{species}/{genome_name}.fa")
                # cmd_tmp = "\n".join(cmd_tmp_list)
                # cmd2 = f"echo '{cmd_tmp}' | xargs -I{{}} -P {threads} bash -c '{{}}'"
                # cmd.append(cmd2)
                # all_new_genomes_path = " ".join(new_genomes_path)
                # cmd3 = f"cat {all_new_genomes_path} | bgzip -c -@ {threads} > {self.wd}/{species}/{species}_merged.fa.gz"
                # cmd.append(cmd3)
                # cmd4 = f"samtools faidx {self.wd}/{species}/{species}_merged.fa.gz"
                # cmd.append(cmd4)

            if self.is_show_detailed_log:
                time_log = f"/usr/bin/time -v -o {self.wd}/{species}/{species}_pangenome_building_time.log "
            else:
                time_log = ""
            if "pggb" in self.pangenome_building_exe.lower():
                cmd5 = f"{time_log}{self.pangenome_building_exe} -i {self.wd}/{species}/{species}_merged.fa.gz -o {self.wd}/{species}/species{species}_pangenome_building -t {threads} -p 90 -n {genomes_num} -v > /dev/null 2>&1"
                cmd.append(cmd5)
            elif "cactus-pangenome" in self.pangenome_building_exe.lower():
                # cmd5_0 = f"find {self.wd}/{species} -type f -name '*.fa.gz' | awk -F'/' '{{file=$NF; sub(/\\.fa\\.gz$/, \"\", file); print file \"\\t\" $0}}' > {self.wd}/{species}/genome2id.tsv"
                cmd5_0 = f"awk -F'/' '{{file=$NF; sub(/\\..*/, \"\", file); print file \"\\t\" $0}}' {pan_species_file} > {self.wd}/{species}/genome2id.tsv"
                cmd.append(cmd5_0)
                if self.reference:
                    reference = self.species2reference[species]
                    cmd5_1 = f"{time_log}{self.pangenome_building_exe} {self.wd}/{species}/js {self.wd}/{species}/genome2id.tsv --outDir {self.wd}/{species}/species{species}_pangenome_building --outName {species} --reference {reference} --maxCores {threads} --mapCores {max(int(threads/4), 6)} > {self.wd}/{species}/{species}_mc.log 2>&1"
                    cmd.append(cmd5_1)
                else:    
                    cmd5_1 = f"reference=$(awk 'NR==1 {{print $1}}' {self.wd}/{species}/genome2id.tsv)"
                    cmd.append(cmd5_1)
                    cmd5_2 = f"{time_log}{self.pangenome_building_exe} {self.wd}/{species}/js {self.wd}/{species}/genome2id.tsv --outDir {self.wd}/{species}/species{species}_pangenome_building --outName {species} --reference $reference --maxCores {threads} --mapCores {max(int(threads/4), 6)} > {self.wd}/{species}/{species}_mc.log 2>&1"
                    cmd.append(cmd5_2)
                    cmd5_3 = f"gunzip {self.wd}/{species}/species{species}_pangenome_building/{species}.gfa.gz"
                    cmd.append(cmd5_3)
            cmd6 = f"{self.vg} convert -g {self.wd}/{species}/species{species}_pangenome_building/*gfa -p -t {threads} > {self.wd}/{species}.vg; mv {self.wd}/{species}/species{species}_pangenome_building/*gfa {self.wd}/{species}.gfa"
            cmd.append(cmd6)
            if self.save:
                cmd_tmp = f"python {script_dir}/get_h5_from_gfa.py {self.wd}/{species}.gfa {self.wd}/{species}.h5; rm {self.wd}/{species}.gfa"
                cmd.append(cmd_tmp)
            if not self.debug:
                cmd7 = f"rm -rf {self.wd}/{species}"
                cmd.append(cmd7)
            cmd8 = f"echo {species}"
            cmd.append(cmd8)
            self.all_cmd.append(("; ".join(cmd), threads, species))
            # print("; ".join(cmd))
        # Sort the commands based on core requirements, for optimal allocation
        self.all_cmd.sort(key=lambda x: x[1], reverse=True)

    def run_command(self, cmd, num_cores):
        """Run a shell command using the specified number of cores."""
        try:
            self.log.debug(f"Running: {cmd} with {num_cores} cores")
            result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return result.stdout.decode().strip()
        except subprocess.CalledProcessError as e:
            self.log.error(f"Error running command: {cmd}\n{e.stderr.decode()}")
            raise RuntimeError(f"Command failed: {cmd}")

    # async def submit_task(self, cmd, num_cores):
    #     """Asynchronously call run_command in the event loop."""
    #     result = await asyncio.to_thread(self.run_command, cmd, num_cores)
    #     return result        

    async def async_run_command(self, cmd, num_cores):
        """Asynchronous wrapper for running shell commands."""
        # self.log.debug(f"Running: {cmd} with {num_cores} cores")
        # Run the command asynchronously
        process = await asyncio.create_subprocess_shell(
            cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()

        if process.returncode != 0:
            raise Exception(f"Error running command: {cmd}\n{stderr.decode()}")
        
        return stdout.decode().strip()
    
    async def schedule_tasks(self):
        """Schedule and execute all commands."""
        available_cores = self.threads
        last_percentage = 0
        pending_tasks = set()  # Set of tasks currently in progress
        pending_species = set()
        command_cores = {}
        while self.all_cmd or pending_tasks:
             # Submit new tasks until all cores are used
            while self.all_cmd and available_cores > 0:
                cmd, num_cores, species = self.all_cmd[0]
                if available_cores >= num_cores:
                    task = asyncio.create_task(self.async_run_command(cmd, num_cores))
                    pending_tasks.add(task)
                    pending_species.add(species)
                    command_cores[species] = num_cores
                    self.all_cmd.pop(0)
                    available_cores -= num_cores
                    # self.log.debug(f"Submitted task {cmd} with {num_cores} cores")
                    self.log.debug(f"Submitted task {species} with {num_cores} cores")
                else:
                    break  # Not enough cores available
            self.log.debug(f"available_cores: {available_cores}")
            self.log.debug(f"{list(command_cores.keys())} species running")
            # Wait for at least one task to complete
            done, pending_tasks = await asyncio.wait(pending_tasks, return_when=asyncio.FIRST_COMPLETED)
            
            # Process completed tasks
            for task in done:
                try:
                    result = task.result()  # Get the result of the task
                    result_gfa = Path(self.wd) / f"{result}.gfa"
                    result_vg = Path(self.wd) / f"{result}.vg"
                    result_h5 = Path(self.wd) / f"{result}.h5"
                    if not ((is_file_non_empty(result_gfa) or is_file_non_empty(result_h5)) and is_file_non_empty(result_vg)):
                        raise IOError(f"{result} gfa or vg files are either missing or empty.")
                    self.finished_pangenome_num += 1
                    current_percentage = self.finished_pangenome_num * 100 / self.pangenomes_need_to_build
                    if current_percentage > last_percentage + 5:
                        self.log.info(f"Species pangenome building percentage: {round(current_percentage, 1)}%")
                        last_percentage = current_percentage

                    with open(self.finished_pangenome, "a") as f:
                        f.write(result + "\n")
                    pending_species.remove(result)
                    num_cores = command_cores.pop(result)
                    # Release cores
                    available_cores += num_cores
                    self.log.debug(f"Completed {result}, freed {num_cores} cores")
                except Exception as e:
                    self.log.error(f"Error in task execution: {e}")

            # Check if there are still commands left and available cores
            if not pending_tasks and not self.all_cmd:
                await asyncio.sleep(1)

    def parallel_run(self):
        current_percentage = self.finished_pangenome_num * 100 / self.pangenomes_need_to_build
        self.log.info(f"Species pangenome building percentage: {round(current_percentage, 1)}%")
        if current_percentage == 100:
            self.log.info("All species pangenome has been built. Skipping this step.")
        else:
            start_time = time.time()
            self.log.info(f"Started executing pangenome building commands.")
            """Start the asynchronous task scheduler."""
            asyncio.run(self.schedule_tasks())
            end_time = time.time()
            total_time = end_time - start_time
            hours, remainder = divmod(total_time, 3600)
            minutes, seconds = divmod(remainder, 60)
            self.log.info(f"Completed all pangenome building commands for the species more than two genomes. Total execution wall clock time: {int(hours)}h {int(minutes)}min {int(seconds)}s")

    def run(self):
        current_percentage = self.finished_pangenome_num * 100 / self.pangenomes_need_to_build
        self.log.info(f"Species pangenome building percentage: {round(current_percentage, 1)}%")
        if current_percentage == 100:
            self.log.info("All species pangenome has been built. Skipping this step.")
        else:
            start_time = time.time()
            last_percentage = 0
            self.log.info(f"Started executing pangenome building commands.")
            for i, para in enumerate(self.all_cmd):
                cmd, num_cores, species = para
                result = self.run_command(cmd, num_cores)
                self.finished_pangenome_num += 1
                current_percentage = self.finished_pangenome_num * 100 / self.pangenomes_need_to_build
                if current_percentage > last_percentage + 5:
                    self.log.info(f"Species pangenome building percentage: {round(current_percentage, 1)}%")
                    last_percentage = current_percentage
                assert result == species
                with open(self.finished_pangenome, "a") as f:
                    f.write(result + "\n")                
            end_time = time.time()
            total_time = end_time - start_time
            hours, remainder = divmod(total_time, 3600)
            minutes, seconds = divmod(remainder, 60)
            self.log.info(f"Completed all pangenome building commands for the species more than two genomes. Total execution wall clock time: {int(hours)}h {int(minutes)}min {int(seconds)}s")


    def limited_core_pool(self):
        """Run commands with limited total cores using a ProcessPoolExecutor."""

        current_percentage = self.finished_pangenome_num * 100 / self.pangenomes_need_to_build
        self.log.info(f"Species pangenome building percentage: {round(current_percentage, 1)}%")
        if current_percentage == 100:
            self.log.info("All species pangenome has been built. Skipping this step.")
        else:
            start_time = time.time()
            self.log.info(f"Started executing pangenome building commands.")
            # Track core allocation and futures
            available_cores = self.threads
            futures = []
            command_cores = {}
            last_percentage = 0
            flag = 0
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.threads) as executor:
                # first submit
                allocated_tasks = []
                for i, (cmd, num_cores, species) in enumerate(self.all_cmd):
                    if available_cores >= num_cores:
                        future = executor.submit(self.run_command, cmd, num_cores)
                        command_cores[future] = num_cores
                        futures.append(future)
                        available_cores -= num_cores  
                        allocated_tasks.append(i)
                    else:
                        break
                for i in reversed(allocated_tasks):
                    self.all_cmd.pop(i)

                while self.all_cmd or futures:
                    flag += 10
                    for future in concurrent.futures.as_completed(futures[:]):
                        flag += 1
                        self.log.debug(f"futures: {len(futures)} {futures}")
                        self.log.debug(f"flag: {flag}")
                        result = future.result()  
                        with open(self.finished_pangenome, "a") as f:
                            f.write(result + "\n")
                        self.finished_pangenome_num += 1
                        current_percentage = self.finished_pangenome_num * 100 / self.pangenomes_need_to_build
                        if current_percentage > last_percentage + 5:
                            self.log.info(f"Species pangenome building percentage: {round(current_percentage, 1)}%")
                            last_percentage = current_percentage

                        used_cores = command_cores.pop(future)
                        available_cores += used_cores
                        futures.remove(future)
                        self.log.debug(f"Completed {result}, freed {used_cores} cores")
                        self.log.debug(f"futures: {len(futures)} {futures}")
                        self.log.debug(f"flag: {flag}")
                        # Each time the memory is released, it is determined whether a new command can be submitted
                        while self.all_cmd and available_cores > 0:
                            cmd, num_cores, species = self.all_cmd[0]
                            if available_cores >= num_cores:
                                future = executor.submit(self.run_command, cmd, num_cores)
                                command_cores[future] = num_cores
                                available_cores -= num_cores 
                                futures.append(future)
                                self.all_cmd.pop(0)  
                        self.log.debug(f"futures: {len(futures)} {futures}")
                    if not futures and not self.all_cmd:
                        time.sleep(1)
            end_time = time.time()
            total_time = end_time - start_time
            hours, remainder = divmod(total_time, 3600)
            minutes, seconds = divmod(remainder, 60)
            self.log.info(f"Completed all pangenome building commands for the species more than two genomes. Total execution wall clock time: {int(hours)}h {int(minutes)}min {int(seconds)}s")

def main():
    parser = argparse.ArgumentParser(prog="multi_tasks_parallel.py", description=usage)
    parser.add_argument("wd", type=str, help="Pangenome building directory.")
    parser.add_argument("pan_species", type=str, help="Pangenome species information directory.")
    # parser.add_argument("fastix", type=str, help="Fastix executable file.")
    parser.add_argument("pantaxr", type=str, help="pantaxr executable file.")
    parser.add_argument("vg", type=str, help="Vg executable file.")
    parser.add_argument("-e", "--pangenome_building_exe", type=str, default="pggb", help="Pangenome building executable file(PGGB, Minigraph-Cactus).")
    parser.add_argument("-r", "--reference", type=str, help="Reference genomes for each species(Minigraph-Cactus need).")
    parser.add_argument("-p", "--parallel", dest="parallel", type=str, default="True", help="Parallel task.")
    parser.add_argument("-g", "--save", dest="save", type=str, default="False", help="Save GFA file information to h5 file.")
    parser.add_argument("-s", "--sleep", type=int, help="Parallel loop sleep time.")
    parser.add_argument("-t", "--threads", dest="threads", type=int, required=True, help="Max threads.")
    parser.add_argument("-f", "--force", dest="force", type=str, default="False", help="Force all to rebuild.")
    parser.add_argument("-d", "--debug", dest="debug", type=str, default="False", help="Debug mode(not delete work directory).")
    parser.add_argument("-v", "--verbose", dest="is_show_detailed_log", type=str, default="False", help="Show detailed log.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    # print("\nProgram settings:\n")
    # for arg in vars(args):
    #     print(arg, "=", getattr(args, arg))
    # print()
    task_scheduling = TaskScheduling(**vars(args))
    if args.parallel.lower() == "true" or args.parallel.lower() == "on":
        # parallel_scheduling.limited_core_pool()
        task_scheduling.parallel_run()
    else:
        task_scheduling.run()


if __name__ == "__main__":
    sys.exit(main())

