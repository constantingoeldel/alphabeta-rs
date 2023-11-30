use alphabeta::AlphaBeta as Args;
use clap::Parser;
use methylome::steady_state;

fn main() {
    let args = Args::parse();

    let (multi, _) = progress_bars::multi(1);
    // PROGRESS_BARS.get_or_init(|| multi);

    let result = alphabeta::run(args.clone(), &multi);

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
                .to_file(args.output.join("pedigree.txt"))
                .expect("Failed to write pedigree");
            analysis
                .to_file(args.output.join("analysis.txt"))
                .expect("Failed to write results");
        }
    }
}
