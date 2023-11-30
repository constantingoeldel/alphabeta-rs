mod chromosome;
mod methylation;
mod strand;

pub use chromosome::Chromosome;
pub use methylation::MethylationSite;
pub use methylation::MethylationStatus;
pub use strand::Strand;

use std::{io, path::PathBuf};
#[derive(Debug, thiserror::Error)]
pub enum Error {
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
    #[error("Methylation site could not be parsed: Wrong format")]
    MethlyationSiteFormat,

    #[error("Chromosome could not be parsed from this string: {0}")]
    Chromosome(String),
}

pub type Return<T> = std::result::Result<T, Error>;

// TODO: Move to a more appropriate place

/// Calculate the steady state UU level
pub fn p_uu_est(alpha: f64, beta: f64) -> f64 {
    (beta * ((1.0 - beta).powi(2) - (1.0 - alpha).powi(2) - 1.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0))
}

/// Calculate the steady state MM level
pub fn p_mm_est(alpha: f64, beta: f64) -> f64 {
    (alpha * ((1.0 - alpha).powi(2) - (1.0 - beta).powi(2) - 1.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0))
}

/// Calculate the steady state methylation level
pub fn steady_state(alpha: f64, beta: f64) -> f64 {
    let pi_2 = (4.0 * alpha * beta * (alpha + beta - 2.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));

    p_mm_est(alpha, beta) + 0.5 * pi_2
}
