
# GTDB genome ANI estimation and remove
pantax-rg -r /home/work/wenhai/metaprofiling/bacteria_GTDB/data/gtdb_genomes_all/GTDB_complete/2024-12-23_10-47-46/files --gm /home/work/wenhai/test/genome_process_test/test_gtdb_metadata.tsv --compute --cluster

# specify species
pantax-rg -r /home/work/wenhai/metaprofiling/bacteria_GTDB/data/gtdb_genomes_all/GTDB_complete/2024-12-23_10-47-46/files --gm /home/work/wenhai/test/genome_process_test/test_gtdb_metadata.tsv --compute --cluster -sc "Enterococcus_F_willemsii"

# specify example_genomes_info
pantax-rg --custom /home/work/wenhai/test/genome_process_test/example_genomes_info.txt --compute  --cluster

# refseq
pantax-rg -r /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/genomeData-20230803 -s /home/work/wenhai/test/genome_process_test/test_rs_complete_assembly_summary.txt --remove --db rs --compute --cluster

# refseq specify species
pantax-rg -r /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/genomeData-20230803 -s /home/work/wenhai/test/genome_process_test/test_rs_complete_assembly_summary.txt --remove --db rs --compute --cluster -sc 632,1358

# refseq specify species, complete genome
pantax-rg -r /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/genomeData-20230803 --remove --db rs -sc 1282 -l complete

# specify genome with list, hcls cluster method
pantax-rg -c genomes_info_1282_test.txt --cluster-method hcls

# specify species, hcls cluster method
pantax-rg -r /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/genomeData-20230803 -s /home/work/wenhai/test/genome_process_test/test_rs_complete_assembly_summary.txt --remove --db rs -sc 562 --cluster-method hcls

# specify ani matrix 
pantax-rg -c /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/single_species_strain_level_1282_all/database_build/pantax4/genomes_info_1282_species.txt --matrix /home/work/wenhai/metaprofiling/bacteria_refgenome_NCBIdata/single_species_strain_level_1282_all/database_build/pantax4/ani_res.matrix  --cluster-method hcls