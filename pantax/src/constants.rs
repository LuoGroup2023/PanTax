pub const TS_DASHES_BLANK_COLONS_DOT_BLANK: &str = "%Y-%m-%d %H:%M:%S%.3f";
pub const GENOMES_INFO: &str = "genomes_info.txt";
pub const CHUNK_SIZE: usize = 1024;

#[cfg(feature = "gb")]
pub const DEFAULT_SOLVER: &str = "gurobi";

#[cfg(feature = "cp")]
pub const DEFAULT_SOLVER: &str = "cplex";

#[cfg(feature = "glpk")]
pub const DEFAULT_SOLVER: &str = "glpk";

#[cfg(feature = "hs")]
pub const DEFAULT_SOLVER: &str = "highs";

#[cfg(feature = "free")]
pub const DEFAULT_SOLVER: &str = "highs";

// glpk
// https://github.com/firedrakeproject/glpk/blob/master/src/glpk.h
pub const GLP_MIN: i32 = 1;
pub const GLP_MAX: i32 = 2;

pub const GLP_FR: i32 = 1;  // Free variable
pub const GLP_LO: i32 = 2;  // Variable with lower bound
pub const GLP_UP: i32 = 3;  // Variable with upper bound
pub const GLP_DB: i32 = 4;  // Double bounded variable
pub const GLP_FX: i32 = 5;  // Fixed variable

pub const GLP_CV: i32 = 1;  // Continuous variable
pub const GLP_IV: i32 = 2;  // Integer variable
pub const GLP_BV: i32 = 3;  // Binary variable

pub const GLP_ON: i32 = 1;
pub const GLP_OFF: i32 = 0;

pub const CLI_HEADINGS: [&str; 7] = [
    "Basic options",
    "Database creation",
    "Index construction(for vg giraffe)",
    "Read alignment",
    "Profiling",
    "General options",
    "Log options"
];