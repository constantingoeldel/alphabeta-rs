use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("No nodes could be parsed from the nodelist")]
    NodeParsing,
    #[error("Could not open node file {0}")]
    NodeFile(#[from] std::io::Error),
}
