extern crate openblas_src;

mod ab_neutral;
mod analysis;
mod boot_model;
mod config;
mod divergence;
pub(crate) mod optimizer;
mod plot;

pub use analysis::Analysis;
pub use config::AlphaBeta;
use optimizer::Model;
use pedigree::Pedigree;
use progress_bars::MultiProgress;

type ObsSteadyState = f64;
pub type Return<T> = Result<T, Error>;

/// Run AlphaBeta
///
/// Returns:
///
/// * `Model` - The best model found by the algorithm
/// * `Analysis` - The analysis of the model, done by bootstrapping
/// * `Pedigree` - The pedigree used for the analysis
/// * The observed steady state methylation level
pub fn run(
    args: AlphaBeta,
    bars: &MultiProgress,
) -> Return<(Model, Analysis, Pedigree, ObsSteadyState)> {
    config::set(args);

    println!("Building pedigree...");

    let (pedigree, p0uu) = Pedigree::build(
        &config::get().nodes,
        &config::get().edges,
        config::get().posterior_max_filter,
    )?;

    let (pb_neutral, pb_boot) = progress_bars::specific(bars, config::get().iterations);

    let (model, pred_div, residuals) =
        ab_neutral::run(&pedigree, p0uu, p0uu, 1.0, Some(&pb_neutral))?;
    let analysis = boot_model::run(
        &pedigree,
        &model,
        pred_div,
        residuals,
        p0uu,
        p0uu,
        1.0,
        Some(&pb_boot),
    )?;
    bars.remove(&pb_neutral);
    bars.remove(&pb_boot);

    Ok((model, analysis, pedigree, 1.0 - p0uu))
}

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Argument error")]
    Argument(#[from] clap::Error),

    #[error("Error while building pedigree: {0}")]
    Pedigree(#[from] pedigree::Error),

    #[error("Faild to run optimization algorithm: {0}")]
    Optimization(#[from] argmin::core::Error),
    // TODO handle drawing error
    // #[error("Plotting error {0}")]
    // Plot(#[from] plotters::drawing::DrawingAreaErrorKind<_>),
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
