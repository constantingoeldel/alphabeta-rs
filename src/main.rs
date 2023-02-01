mod ABneutral;
mod divergence;
mod pedigree;
extern crate blas_src;
fn main() {
    let pedigree =
        pedigree::read_pedigree_from_file("/home/cgoeldel/epigenomics/alphabeta/pedigree.txt");

    let result = ABneutral::run(pedigree, 0.75, 0.5, 0.7, 2);
    println!("{result:?}");
}
