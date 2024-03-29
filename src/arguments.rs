use clap::{command, Parser, Subcommand};
use std::path::PathBuf;
use std::time::SystemTime;

/// simple tool to separate a methylome by position within a gene
#[derive(Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
#[non_exhaustive]
pub struct Windows {
    /// Path of directory containing the methlyome files from which to extract the CG-sites
    #[arg(short, long)]
    pub methylome: PathBuf,

    /// Path of the annotation file containing information about beginning and end of gbM-genes
    #[arg(short, long)]
    pub genome: PathBuf,

    /// Size of the window in percent of the gbM-gene length or in basepair number if --absolute is supplied
    #[arg(short, long, default_value_t = 5)]
    pub window_size: u32,

    /// Size of the step between the start of each window. Default value is window-size, so no overlapp happens
    #[arg(long, short('s'), default_value_t = 0)]
    pub window_step: u32,

    /// Path of the directory where extracted segments shall be stored
    #[arg(short, long, default_value_os_t = PathBuf::from("."), value_parser = validate_default_output_dir)]
    pub output_dir: PathBuf,

    /// Use absolute length in base-pairs for window size instead of percentage of gene length
    #[arg(short, long, default_value_t = false)]
    pub absolute: bool,

    /// Number of basepairs to include upstream and downstream of gene
    #[arg(short, long, default_value_t = 2048)]
    pub cutoff: u32,

    /// Invert strands, to switch from 5' to 3' and vice versa
    #[arg(short, long, default_value_t = false)]
    pub invert: bool,

    /// Name of the run to be used when storing the result in Postgres
    #[arg(long, default_value_t = format!("Anonymous Run {}", SystemTime::now().duration_since(SystemTime::UNIX_EPOCH).unwrap().as_secs()))]
    pub name: String,

    /// Overwrite existing content in output directory? If false (default) it will reuse existing windows
    #[arg(long, short, default_value_t = true)]
    pub force: bool,
    /// Let the cutoff be the gene length instead of a fixed number.
    /// So if the gene is 1000 bp long, the cutoff will be 1000 bp instead of 2048 bp (the default).
    /// This option takes preference over the cutoff option.    
    #[arg(long, default_value_t = false)]
    pub cutoff_gene_length: bool,

    /// Number of iterations to run for Nelder-Mead optimization, both in the model estimation and bootstrap phase
    /// Default is 100, you can probably get away with 10
    #[arg(long, default_value_t = 100)]
    pub iterations: usize,

    #[command(subcommand)]
    pub command: Option<Subcommands>,
}

#[derive(Subcommand, Debug, Clone)]
#[clap(rename_all = "lower")]
pub enum Subcommands {
    /// Enable AlphaBeta estimation on the extracted windows
    AlphaBeta(AlphaBeta),
}
impl Default for Windows {
    fn default() -> Self {
        Windows {
            invert: false,
            absolute: false,
            cutoff: 2048,
            genome: PathBuf::from("./genome"),
            methylome: PathBuf::from("./methylome"),
            output_dir: PathBuf::from("./"),
            window_size: 5,
            window_step: 1,
            command: Some(Subcommands::AlphaBeta(AlphaBeta::default(
                PathBuf::from("."),
                100,
            ))),
            name: String::new(),
            force: false,
            cutoff_gene_length: false,
            iterations: 100,
        }
    }
}

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
