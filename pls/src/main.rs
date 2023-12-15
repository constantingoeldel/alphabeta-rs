use polars::prelude::*;

fn main() {
    println!("Hello, world!");
    let tic = std::time::Instant::now();
    let df = CsvReader::from_path("../data/methylome/G0.txt")
        .unwrap()
        .finish()
        .unwrap();
    println!("Elapsed time: {:?}", tic.elapsed());
    println!("{}", df);
}
