use anyhow::anyhow;
use indicatif::MultiProgress;

use crate::{
    analysis::{Analysis, RawAnalysis},
    arguments::AlphaBeta as Args,
    pedigree::Pedigree,
    progress::specific,
    structs::Model,
    *,
};

type ObsSteadyState = f64;

/// Run AlphaBeta
///
/// Returns:
///
/// * `Model` - The best model found by the algorithm
/// * `Analysis` - The analysis of the model, done by bootstrapping
/// * `Pedigree` - The pedigree used for the analysis
/// * The observed steady state methylation level
pub fn run(
    args: Args,
    bars: &MultiProgress,
) -> Result<(Model, Analysis, RawAnalysis, Pedigree, ObsSteadyState)> {
    println!("Building pedigree...");
    let (pedigree, p0uu) = Pedigree::build(&args.nodes, &args.edges, args.posterior_max_filter)
        .map_err(|e| anyhow!("Error while building pedigree: {}", e))?;

    let (pb_neutral, pb_boot) = specific(bars, args.iterations);

    let (model, pred_div, residuals) = ab_neutral::run(
        &pedigree,
        p0uu,
        p0uu,
        1.0,
        args.iterations,
        Some(&pb_neutral),
    )
    .map_err(|e| anyhow!("Model failed: {}", e))?;
    let (analysis, raw_analysis) = boot_model::run(
        &pedigree,
        &model,
        pred_div,
        residuals,
        p0uu,
        p0uu,
        1.0,
        args.iterations,
        Some(&pb_boot),
        &args.output,
    )
    .map_err(|e| anyhow!("Bootstrap failed: {}", e))?;
    bars.remove(&pb_neutral);
    bars.remove(&pb_boot);

    Ok((model, analysis, raw_analysis, pedigree, 1.0 - p0uu))
}

/// Calculate the steady state UU level
pub fn p_uu_est(alpha: f64, beta: f64) -> f64 {
    (beta * ((1.0 - beta).powi(2) - (1.0 - alpha).powi(2) - 1.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0))
}

/// Calculate the steady state MM level
pub fn p_mm_est(alpha: f64, beta: f64) -> f64 {
    (alpha * ((1.0 - alpha).powi(2) - (1.0 - beta).powi(2) - 1.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0))
}

/// Calculate the steady state methylation level
pub fn steady_state(alpha: f64, beta: f64) -> f64 {
    let pi_2 = (4.0 * alpha * beta * (alpha + beta - 2.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));

    p_mm_est(alpha, beta) + 0.5 * pi_2
}

#[cfg(test)]
mod tests {
    // use super::run;
    // use crate::{
    //     ab_neutral, arguments::AlphaBeta, assert_close, assert_within_10_percent, boot_model,
    //     pedigree::Pedigree,
    // };
    // use std::path::{Path, PathBuf};

    // #[test] // Recommended to run with --release
    // fn end_to_end() {
    //     let (pedigree, p0uu) = Pedigree::build(
    //         Path::new("./data/nodelist.txt"),
    //         Path::new("./data/edgelist.txt"),
    //         0.99,
    //     )
    //     .expect("Error while building pedigree: ");
    //     dbg!(&pedigree);
    //     let output_dir = PathBuf::from("./");
    //     let (model, pred_div, residuals) =
    //         ab_neutral::run(&pedigree, p0uu, p0uu, 1.0, 1000, None).expect("Model failed");
    //     let _res = boot_model::run(
    //         &pedigree,
    //         &model,
    //         pred_div,
    //         residuals,
    //         p0uu,
    //         p0uu,
    //         1.0,
    //         200,
    //         None,
    //         &output_dir,
    //     )
    //     .expect("Bootstrap failed");
    //     dbg!(&model);
    //     // Reference data from https://github.com/jlab-code/AlphaBeta/blob/master/inst/extdata/dm/Col_CG_global_estimates_ABneutral.Rdata
    //     assert_within_10_percent!(model.beta, 0.0017179248);
    //     assert_within_10_percent!(model.alpha, 2.298873e-04);
    //     assert_within_10_percent!(p0uu, 0.8811696);
    //     assert_close!(model.alpha, 2.298873e-04);
    //     assert_close!(model.beta, 0.0017179248);
    //     assert_close!(p0uu, 0.8811696);
    // }

    // #[test]
    // fn random_window() {
    //     let pb = indicatif::MultiProgress::new();
    //     let res = run(AlphaBeta {
    //         edges: PathBuf::from("/mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/wt/gene/50/edgelist.txt"),
    //         nodes: PathBuf::from("/mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/wt/gene/50/nodelist.txt"),
    //         iterations: 100,
    //         posterior_max_filter: 0.99,
    //         output: PathBuf::from("/mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/wt/gene/50/")
    //     }, &pb);

    //     assert!(res.is_ok());

    //     println!("{:?}", res.unwrap().1);
    // }
}
