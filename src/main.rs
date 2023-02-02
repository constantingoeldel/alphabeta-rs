mod ABneutral;
mod BootModel;
mod divergence;
mod pedigree;
mod structs;
extern crate blas_src;

pub use divergence::*;
pub use structs::*;

fn main() {
    let pedigree = Pedigree::from_file("/home/cgoeldel/epigenomics/alphabeta/pedigree.txt");

    let model = ABneutral::run(pedigree, 0.5, 0.7, 0.75, 2).expect("Model failed");
    let result = BootModel::run(pedigree, 0.5, 0.7, 0.75, model, 1000).expect("Bootstrap failed");
    println!("{result:?}");
}
