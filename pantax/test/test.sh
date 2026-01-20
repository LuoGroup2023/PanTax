
pantax=/path/to/pantax

cd example/ngs
#------------------------------------------
# the command to directly perform profiling
#------------------------------------------
$pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --species --strain -t 32

#------------------------------------------
# the commands to perform profiling step by step
#------------------------------------------
# create database only
$pantax -f ../example_genomes_info.txt --create -t 32

# create index only
$pantax -f ../example_genomes_info.txt --index -t 32

# perform profiling
$pantax -s -p -r short_reads.fq.gz --species --strain -t 32

#------------------------------------------
# the command to directly perform profiling with fast mode
#------------------------------------------
$pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --species --strain -t 32 --fast

#------------------------------------------
# other commands
#------------------------------------------
# create and index
$pantax -f ../example_genomes_info.txt --create --index

# fast mode construction and index
$pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --create --index --fast 

# If you need to process multiple samples and run in fast mode, 
# it is recommended to build the genome index files in advance and pass them via parameters.
awk -F'\t' 'NR>1 {print $5}' ../example_genomes_info.txt > genomes.txt
sylph sketch -l genomes.txt -o genomes
$pantax -f ../example_genomes_info.txt -s -p -r short_reads.fq.gz --create --index --fast --syldb genomes.syldb

# long read
# same as short read but use -l instead of -s
cd example/hifi
#------------------------------------------
# the command to directly perform profiling
#------------------------------------------
$pantax -f ../example_genomes_info.txt -l -r long_reads.fq.gz --species --strain -t 32

# use vg as long read aligner, need specify long read type with --lt (ontr10 or hifi)
# please use vg >= 1.71, not test lower version
$pantax -f ../example_genomes_info.txt -l -r long_reads.fq.gz -t 32 --species --strain --lr-aligner /path/to/vg --lt hifi

# If you encounter difficulties installing or running Gurobi, you can use the open-source solver version of PanTax (pantax_free).
# pantax      ------> support gurobi (default), highs, cbc, glpk
# pantax_free ------> support highs (default), cbc, glpk 
# pantax_gb   ------> support gurobi 
# pantax_cp   ------> support cplex  (Note that it need to build from source)
pantax_free -f ../example_genomes_info.txt -l -r long_reads.fq.gz --species --strain -t 32

