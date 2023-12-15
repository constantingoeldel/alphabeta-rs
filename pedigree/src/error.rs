use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("No nodes could be parsed from the nodelist")]
    NodeParsing,
    #[error("Could not open node file {0}")]
    NodeFile(#[from] std::io::Error),
    #[error("The edge file contains the node {0} that does not exist in the node file")]
    EdgelistContainingNonExistentNode(String),
    #[error("The time difference {0} could not be parsed")]
    TimeDifferenceUnparseable(#[from] std::num::ParseIntError),
    #[error("The edge file does not contain a 'from' column and values")]
    NoFrom,
    #[error("The edge file does not contain a 'to' column and values")]
    NoTo,
    #[error("The nodelist line does not contain a filename\n{0}")]
    NoFile(String),
    #[error("The nodelist line does not contain a node name\n{0}")]
    NoName(String),
    #[error("The nodelist line does not contain a time generation\n{0}")]
    NoGeneration(String),
    #[error("The generation {0} could not be parsed")]
    GenerationUnparseable(String),
    #[error("The nodelist line does not contain a methylation indicator like 'Y' or 'N'\n{0}")]
    NoMethylation(String),
    #[error("The methylation data file {0} could not be opened")]
    MethylationFile(PathBuf),
    #[error("The methylation site {0} could not be parsed")]
    MethylationParsingError(String),
}
