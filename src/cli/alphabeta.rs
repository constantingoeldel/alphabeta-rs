use std::path::Path;

use alphabeta::*;
use clap::Parser;

fn main() {
    let mut args = Args::parse();

    let (pedigree, p0uu) =
        Pedigree::build(&args.nodelist, &args.edgelist, args.posterior_max_filter)
            .expect("Error while building pedigree: ");
    pedigree
        .to_file(Path::new("./data/pedigree_wt.txt"))
        .unwrap();

    let model =
        ABneutral::run(&pedigree, p0uu, p0uu, 1.0, args.iterations, None).expect("Model failed");
    let result = BootModel::run(&pedigree, &model, p0uu, p0uu, 1.0, args.iterations, None)
        .expect("Bootstrap failed");

    println!("##########");
    println!("Results:\n");
    println!("{model}");
    println!("{result}");
    println!("##########");
    args.output.push("results.txt");
    model.to_file(&args.output, &result).unwrap();
}
