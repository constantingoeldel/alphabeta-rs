use std::{io, path::PathBuf};
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Argument error")]
    Argument(#[from] clap::Error),

    #[error("Could not find the specified file or directory! Does it exist? \nPath: {0}")]
    File(PathBuf),

    #[error("File error {0}")]
    FileSystem(#[from] io::Error),
    #[error("Unable to extract CG site from line")]
    CGSite,

    #[error("Unable to convert: Are you passing a valid number? {0}")]
    NumberConversion(#[from] std::num::ParseIntError),

    #[error("Unable to convert: Are you passing a valid number? {0}")]
    FloatConversion(#[from] std::num::ParseFloatError),

    #[error("{0}")]
    Simple(&'static str),
}
