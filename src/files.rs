use std::{
    fs::{self, File},
    io::{self, BufRead},
    path::{Path, PathBuf},
};

use anyhow::bail;

use crate::error::Error;

use crate::*;

pub fn open_file(path: &PathBuf) -> Result<File> {
    let file = File::open(path).or(Err(Error::File(path.to_owned())))?;
    Ok(file)
}

pub fn lines_from_file(path: &PathBuf) -> Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(path).or(Err(Error::File(path.to_owned())))?;
    Ok(io::BufReader::new(file).lines())
}

pub fn load_methylome(methylome: &PathBuf) -> Result<Vec<PathBuf>> {
    let methylome_dir = fs::read_dir(methylome).or(Err(Error::File(methylome.to_owned())))?;
    let methylome_files: Vec<PathBuf> = methylome_dir
        .map(|f| f.as_ref().unwrap().path())
        .filter(|path| {
            path.extension().is_some()
                && !path.extension().unwrap().to_str().unwrap().contains("tsv")
                && !path.extension().unwrap().to_str().unwrap().contains("fn")
        }) // Filter out tsv and fn files, which are often nodelist/edgelist files.
        .collect();
    if methylome_files.is_empty() {
        bail!("Could not find any files in the methylome directory. Please check your input. Files with .tsv or .fn extensions are ignored.")
    } else {
        Ok(methylome_files)
    }
}
/// Highly sacriligous
pub fn file_name(path: &Path) -> String {
    let name = path.components().nth_back(0).unwrap();
    name.as_os_str().to_string_lossy().into()
}
