use alphabeta::AlphaBeta as Args;
use clap::Parser;
use methylome::steady_state;

fn main() {
    let args = Args::parse();

    let result = alphabeta::run(args.clone());

    match result {
        Err(e) => println!("{e}"),
        Ok((model, analysis, pedigree, obs_steady_state)) => {
            println!();
            println!();
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
                .to_file(args.output.join("pedigree.txt"))
                .expect("Failed to write pedigree");
            analysis
                .to_file(args.output.join("analysis.txt"))
                .expect("Failed to write results");
        }
    }
}
