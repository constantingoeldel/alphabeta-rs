mod config;
mod extract;
mod files;
mod genes;
mod plot;
mod setup;
mod windows;

pub use config::Metaprofile;
pub use config::Subcommands;
pub use extract::extract;
pub type Return<T> = Result<T, Error>;

use std::{io, path::PathBuf};
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Argument error")]
    Argument(#[from] clap::Error),

    #[error("Could not find any files in the methylome directory. Please check your input. Files with .tsv or .fn extensions are ignored.")]
    NoMethylomeFiles,

    #[error("Could not parse a single annotation from the annotation file. Please check your input or add a parser implemenation for your data format.")]
    NoGenesFound,

    #[error("Could not find the specified file or directory! Does it exist? \nPath: {0}")]
    File(PathBuf),

    #[error("Methylome error: {0}")]
    Methylome(#[from] methylome::Error),

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
    #[error("Methylation site could not be parsed: Wrong format")]
    MethlyationSiteFormat,

    #[error("Chromosome could not be parsed from this string: {0}")]
    Chromosome(String),
}

#[macro_export]
macro_rules! print_dev {
    ($($rest:tt)*) => {
        #[cfg(debug_assertions)]
        std::println!($($rest)*)
    }
}
