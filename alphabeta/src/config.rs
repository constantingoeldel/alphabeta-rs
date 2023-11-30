use clap::{command, Parser};
use std::path::PathBuf;
use std::sync::OnceLock;

#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
#[non_exhaustive]
pub struct AlphaBeta {
    /// Number of iterations to run for Nelder-Mead optimization, even 100 is enough
    #[arg(short, long, default_value_t = 1000)]
    pub iterations: usize,

    /// Relative or absolute path to an edgelist, see /data for an example
    #[arg(long, short, default_value_os_t = PathBuf::from("./edgelist.txt"), value_parser = validate_default_file_existence)]
    pub edges: std::path::PathBuf,

    /// Relative or absolute path to a nodelist, see /data for an example
    #[arg(long, short, default_value_os_t = PathBuf::from("./nodelist.txt"), value_parser = validate_default_file_existence)]
    pub nodes: std::path::PathBuf,
    /// Minimum posterior probability for a singe basepair read to be included in the estimation
    #[arg(long, short, default_value_t = 0.99)]
    pub posterior_max_filter: f64,
    /// Relative or absolute path to an output directory, must exist, EXISTING FILES WILL BE OVERWRITTEN
    #[arg(long, short, default_value_os_t = PathBuf::from("."), value_parser = validate_default_output_dir)]
    pub output: std::path::PathBuf,
}

pub static ALPHABETA_ARGS: OnceLock<AlphaBeta> = OnceLock::new();

pub fn get() -> &'static AlphaBeta {
    ALPHABETA_ARGS
        .get()
        .expect("Config must be initialized at this point")
}

pub fn set(args: AlphaBeta) {
    ALPHABETA_ARGS.get_or_init(|| args);
}

fn validate_default_output_dir(s: &str) -> Result<PathBuf, String> {
    if PathBuf::from(s).exists() {
        println!(
            "Using default output directory: {}",
            // Display full path
            PathBuf::from(s).canonicalize().unwrap().display()
        );
        Ok(PathBuf::from(s))
    } else {
        Err(format!(
            "Please provide a valid output directory. By default, we fill try {s}, which does not exist."
        ))
    }
}

fn validate_default_file_existence(s: &str) -> Result<PathBuf, String> {
    if PathBuf::from(s).exists() {
        println!("Using default file: {}", PathBuf::from(s).display());
        Ok(PathBuf::from(s))
    } else {
        Err(format!(
            "Please provide a valid file path. By default, we fill try {s}, which does not exist."
        ))
    }
}

impl AlphaBeta {
    pub fn default(output_dir: PathBuf, iterations: usize) -> Self {
        Self {
            edges: output_dir.join("edgelist.txt"),
            nodes: output_dir.join("nodelist.txt"),
            output: output_dir,
            posterior_max_filter: 0.99,
            iterations,
        }
    }
}
