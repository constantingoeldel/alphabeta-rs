use clap::Parser;
use pedigree::Pedigree;
use std::path::PathBuf;

fn main() {
    let args = PedigreeArgs::parse();
    let pedigree =
        Pedigree::build(args.nodes, args.edges, args.posterior_max_filter, None).unwrap();
}

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
#[non_exhaustive]
pub struct PedigreeArgs {
    /// Relative or absolute path to an edgelist, see /data for an example
    #[arg(long, short, default_value_os_t = PathBuf::from("./edgelist.txt"), value_parser = validate_default_edgelist_existence)]
    pub edges: std::path::PathBuf,

    /// Relative or absolute path to a nodelist, see /data for an example
    #[arg(long, short, default_value_os_t = PathBuf::from("./nodelist.txt"), value_parser = validate_default_nodelist_existence)]
    pub nodes: std::path::PathBuf,
    /// Minimum posterior probability for a singe basepair read to be included in the estimation
    #[arg(long, short, default_value_t = 0.99)]
    pub posterior_max_filter: f64,
    /// Relative or absolute path to an output directory, must exist, EXISTING FILES WILL BE OVERWRITTEN
    #[arg(long, short, default_value_os_t = PathBuf::from("."), value_parser = validate_default_output_dir)]
    pub output: std::path::PathBuf,
}

#[derive(Debug, thiserror::Error)]
enum Error {
    #[error("Please provide a valid output directory. By default, we will try {0}, which does not exist.")]
    OutputDir(String),

    #[error("Please provide a valid path to an edgelist file. By default, we will try {0}, which does not exist.")]
    Edgelist(String),

    #[error("Please provide a valid path to a nodelist file. By default, we will try {0}, which does not exist.")]
    Nodelist(String),
    // TODO handle drawing error
    // #[error("Plotting error {0}")]
    // Plot(#[from] plotters::drawing::DrawingAreaErrorKind<_>),
}

fn validate_default_output_dir(s: &str) -> Result<PathBuf, Error> {
    if PathBuf::from(s).exists() {
        println!(
            "Using output directory: {} ✅",
            // Display full path
            PathBuf::from(s).canonicalize().unwrap().display()
        );
        Ok(PathBuf::from(s))
    } else {
        Err(Error::OutputDir(s.to_string()))
    }
}

fn validate_default_edgelist_existence(s: &str) -> Result<PathBuf, Error> {
    if PathBuf::from(s).exists() {
        println!("Using edgelist: {} ✅", PathBuf::from(s).display());
        Ok(PathBuf::from(s))
    } else {
        Err(Error::Edgelist(s.to_string()))
    }
}

fn validate_default_nodelist_existence(s: &str) -> Result<PathBuf, Error> {
    if PathBuf::from(s).exists() {
        println!("Using nodelist: {} ✅", PathBuf::from(s).display());
        Ok(PathBuf::from(s))
    } else {
        Err(Error::Nodelist(s.to_string()))
    }
}
