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
        Ok((model, analysis, pedigree, obs_steady_state)) => {
            println!("##########");
            println!("Results:\n");
            println!("{model}");
            println!("{analysis}");
            println!(
                "Estimated steady state {}",
                steady_state(model.alpha, model.beta)
            );
            println!("Observed steady state methylation {obs_steady_state}");
            println!("##########");
            pedigree
                .to_file(&args.output.join("pedigree.txt"))
                .expect("Failed to write pedigree");
            analysis
                .to_file(&args.output.join("analysis.txt"))
                .expect("Failed to write results");
        }
    }
}
