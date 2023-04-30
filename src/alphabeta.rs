use anyhow::anyhow;
use indicatif::MultiProgress;

use crate::{
    arguments::AlphaBeta as Args,
    pedigree::Pedigree,
    progress::specific,
    structs::{Model, StandardDeviations},
    *,
};

pub fn run(args: Args, bars: &MultiProgress) -> Result<(Model, StandardDeviations, Pedigree, f64)> {
    println!("Building pedigree...");
    let (pedigree, p0uu) = Pedigree::build(&args.nodes, &args.edges, args.posterior_max_filter)
        .map_err(|e| anyhow!("Error while building pedigree: {}", e))?;

    let (pb_neutral, pb_boot) = specific(bars, args.iterations);

    let model = ABneutral::run(
        &pedigree,
        p0uu,
        p0uu,
        1.0,
        args.iterations,
        Some(pb_neutral.clone()),
    )
    .map_err(|e| anyhow!("Model failed: {}", e))?;
    let result = BootModel::run(
        &pedigree,
        &model,
        p0uu,
        p0uu,
        1.0,
        args.iterations,
        Some(pb_boot.clone()),
    )
    .map_err(|e| anyhow!("Bootstrap failed: {}", e))?;
    bars.remove(&pb_neutral);
    bars.remove(&pb_boot);

    Ok((model, result, pedigree, 1.0 - p0uu))
}

pub fn steady_state(alpha: f64, beta: f64) -> f64 {
    let pi_1 = (alpha * ((1.0 - alpha).powi(2) - (1.0 - beta).powi(2) - 1.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));

    let pi_2 = (4.0 * alpha * beta * (alpha + beta - 2.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));

    pi_1 + 0.5 * pi_2
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{assert_close, assert_within_10_percent, pedigree::Pedigree, ABneutral, BootModel};

    #[test] // Recommended to run with --release
    fn end_to_end() {
        let (pedigree, p0uu) = Pedigree::build(
            Path::new("./data/nodelist.txt"),
            Path::new("./data/edgelist.txt"),
            0.99,
        )
        .expect("Error while building pedigree: ");

        let model = ABneutral::run(&pedigree, p0uu, p0uu, 1.0, 1000, None).expect("Model failed");
        let result = BootModel::run(&pedigree, &model, p0uu, p0uu, 1.0, 200, None)
            .expect("Bootstrap failed");
        println!("{result}");
        assert_within_10_percent!(model.alpha, 5.7985750419976e-05);
        assert_within_10_percent!(model.beta, 0.00655710970515347);
        assert_within_10_percent!(p0uu, 0.991008120326199);
        assert_close!(model.alpha, 5.7985750419976e-05);
        assert_close!(model.beta, 0.00655710970515347);
        assert_close!(p0uu, 0.991008120326199);

        // assert_close_10_percent!(result.alpha, 1.19025049000535e-05);
        // assert_close_10_percent!(result.beta, 0.00138056339729951);
        // assert_close_10_percent!(result.alpha_beta, 0.594234575518036);
    }
}
