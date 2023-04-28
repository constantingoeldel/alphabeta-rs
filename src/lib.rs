pub mod ABneutral;
pub mod BootModel;
pub mod alphabeta;
pub mod arguments;
pub mod divergence;
pub mod error;
pub mod extract;
pub mod files;
pub mod genes;
pub mod macros;
pub mod methylation_site;
pub mod pedigree;
pub mod progress;
pub mod setup;
pub mod structs;
pub mod windows;

extern crate blas_src;

use anyhow::Result;
use error::Error;
