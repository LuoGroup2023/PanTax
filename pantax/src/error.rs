use grb::Error as GurobiError;
use highs::HighsStatus;
use std::fmt::Debug;
#[derive(Debug)]
pub enum CbcError {
    Status(i32),
}

pub enum GlpkError {
    Status(i32),
}

impl std::fmt::Debug for GlpkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GlpkError::Status(code) => write!(f, "GLPK returned status code {}", code),
        }
    }
}

impl std::fmt::Display for GlpkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl std::fmt::Display for CbcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CbcError::Status(code) => write!(f, "CBC returned status code {}", code),
        }
    }
}

impl std::error::Error for GlpkError {}

#[derive(Debug)]
pub enum SolverError {
    Gurobi(GurobiError),
    Cbc(CbcError),
    Glpk(GlpkError),
    Highs(HighsStatus),
    Other(String),
}

impl std::error::Error for SolverError {}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolverError::Gurobi(e) => write!(f, "Gurobi error: {}", e),
            SolverError::Cbc(e) => write!(f, "CBC error: {:?}", e),
            SolverError::Glpk(e) => write!(f, "Glpk error: {:?}", e),
            SolverError::Highs(e) => write!(f, "Highs error: {:?}", e),
            SolverError::Other(msg) => write!(f, "Other solver error: {}", msg),
        }
    }
}

impl From<GurobiError> for SolverError {
    fn from(e: GurobiError) -> Self {
        SolverError::Gurobi(e)
    }
}

impl From<CbcError> for SolverError {
    fn from(e: CbcError) -> Self {
        SolverError::Cbc(e)
    }
}

impl From<GlpkError> for SolverError {
    fn from(e: GlpkError) -> Self {
        SolverError::Glpk(e)
    }
}

impl From<HighsStatus> for SolverError {
    fn from(e: HighsStatus) -> Self {
        SolverError::Highs(e)
    }
}

