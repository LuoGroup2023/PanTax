

#[derive(Default, Debug)]
pub struct CheckPoints {
    pub reconstruction: bool,
    pub fast_query: bool,
    pub fast_query_and_filter: bool,
}

#[derive(Clone)]
pub struct GenomesInfo {
    pub genome_id: String,
    pub strain_taxid: String,
    pub species_taxid: String,
    pub organism_name: String,
    pub path_id: String,
}

#[derive(Debug, PartialEq)]
pub enum DataType {
    ShortReadSingle,
    ShortReadPaired,
    ShortReadPairedInter,
    LongReadSingle,
}
