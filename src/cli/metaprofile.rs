use std::fs;

use alphabeta::{
    alphabeta::steady_state, arguments::Windows as Args, extract::extract, genes::Region, plot,
    progress, structs::Analysis,
};
use clap::Parser;

fn main() {
    let args = Args::parse();
    println!("Starting run {}", args.name);

    let result = extract(args.clone());

    match (result, args.alphabeta) {
        (Err(e), _) => println!("Error: {e}"),
        (Ok(_), false) => println!("Done"),
        (Ok((max_gene_length, distribution)), true) => {
            alphabeta_multiple(args, max_gene_length, distribution)
        }
    }
}

fn alphabeta_multiple(args: Args, max_gene_length: u32, distribution: Vec<i32>) {
    let regions = vec![
        (Region::Upstream, args.cutoff),
        (Region::Gene, max_gene_length),
        (Region::Downstream, args.cutoff),
    ];
    let mut results = Vec::new();
    let total_steps = if args.absolute {
        (max_gene_length + 2 * args.cutoff) / args.window_step
    } else {
        (3 * 100) / args.window_step
    };

    let (multi, pb) = progress::multi(total_steps as u64);
    for region in regions {
        let max = if args.absolute { region.1 } else { 100 };

        for window in (0..max).step_by(args.window_step as usize) {
            pb.inc(1);

            let args = alphabeta::arguments::AlphaBeta::default(
                args.output_dir
                    .join(region.0.to_string())
                    .join(window.to_string()),
            );

            let alphabeta_result = alphabeta::alphabeta::run(args, &multi);
            match alphabeta_result {
                Err(e) => println!("Error: {e}"),
                Ok((model, deviations, _, obs_meth_lvl)) => {
                    results.push((model, deviations, region.0.clone(), obs_meth_lvl))
                }
            }
        }
    }
    pb.finish();
    let mut print = String::from("run;window;cg_count;region;alpha;beta;alpha_error;beta_error;1/2*(alpha+beta);pred_steady_state;obs_steady_state\n");

    for (i, ((model, analysis, region, obs_meth_lvl), d)) in
        results.iter().zip(distribution.iter()).enumerate()
    {
        print += &format!(
            "{};{};{};{};{};{};{};{};{};{};{}\n",
            args.name,
            i,
            d,
            region,
            model.alpha,
            model.beta,
            analysis.alpha,
            analysis.beta,
            0.5 * (model.alpha + model.beta),
            steady_state(model.alpha, model.beta),
            obs_meth_lvl
        )
    }

    println!("{print}");

    fs::write(args.output_dir.join("results.txt"), print).expect("Could not save results to file.");
    // let db = db::connect()
    //     .await
    //     .expect("Could not connect to database: Did you provide a connection string?");
    // import_results(&db, args.name, results).await.expect("Could not save results to a database. Your data is stored in files in each directory");
    let analyses = results
        .iter()
        .map(|r| r.1.clone())
        .collect::<Vec<Analysis>>();
    plot::metaplot(&analyses, &args).expect("Could not plot results");
}

// #[cfg(test)]
// mod test {
//     use std::path::PathBuf;

//     use alphabeta::arguments::Windows;

//     #[test]
//     fn end_to_end() {
//         let args = Windows {
//             name: String::from("test"),
//             genome: PathBuf::from("data/annotation.bed"),
//             methylome: PathBuf::from("data/methylome"),
//             output_dir: String::from("data/output_metaplot"),
//             alphabeta: false,
//         };
//         }
//     }
// }
