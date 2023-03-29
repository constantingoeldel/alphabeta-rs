use alphabeta::alphabeta::run;
use alphabeta::alphabeta::steady_state;
use alphabeta::{arguments::AlphaBeta as Args, progress};

use clap::Parser;

fn main() {
    let args = Args::parse();

    let (multi, _) = progress::multi(1);

    let result = run(args.clone(), &multi);

    match result {
        Err(e) => println!("Error: {e}"),
        Ok((model, deviations, pedigree)) => {
            println!("##########");
            println!("Results:\n");
            println!("{model}");
            println!("{deviations}");
            println!(
                "Estimated steady state {}",
                steady_state(model.alpha, model.beta)
            );
            println!("##########");
            pedigree
                .to_file(&args.output.join("pedigree.txt"))
                .expect("Failed to write pedigree");
            model
                .to_file(&args.output.join("pedigree.txt"), &deviations)
                .expect("Failed to write results");
        }
    }
}
