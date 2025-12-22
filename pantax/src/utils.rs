
use crate::types::*;
use std::io::{Write, BufWriter};
use std::fs::File;
use std::path::PathBuf;
use anyhow::Result;

pub fn output_genomes_info(genomes_info: &Vec<GenomesInfo>, genomes_info_path: &PathBuf) -> Result<()> {
    let file = File::create(genomes_info_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "genome_ID\tstrain_taxid\tspecies_taxid\torganism_name\tid")?;
    for genome_info in genomes_info.iter() {
        writeln!(writer, "{}\t{}\t{}\t{}\t{}", genome_info.genome_id, genome_info.strain_taxid, genome_info.species_taxid, genome_info.organism_name, genome_info.path_id)?;
    }
    Ok(())
}