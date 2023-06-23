use indicatif::ProgressBar;
use rand::{distributions::Slice, thread_rng, Rng};
use rayon::prelude::*;
use std::{ops::Div, path::Path, sync::Mutex};

use argmin::{core::Executor, solver::neldermead::NelderMead};
use ndarray::{array, s, Array1, Array2, ArrayView1, Axis};
use ndarray_stats::interpolate::Linear;
use ndarray_stats::Quantile1dExt;
use noisy_float::types::n64;

use crate::{
    analysis::{Analysis, RawAnalysis},
    pedigree::Pedigree,
    structs::{Model, PredictedDivergence, Problem, Progress, Residuals},
    *,
};

pub fn run(
    pedigree: &Pedigree,
    params: &Model,
    pred_div: PredictedDivergence,
    residuals: Residuals,
    p0uu: f64,
    eqp: f64,
    eqp_weight: f64,
    n_boot: u64,
    pb: Option<&ProgressBar>,
    output_dir: &Path,
) -> Result<Analysis, Box<dyn std::error::Error>> {
    let alternative_pb = Progress::new("BootModel", n_boot).0;
    let pb = pb.unwrap_or(&alternative_pb);

    let p0mm = 1.0 - p0uu;
    let p0um = 0.0;

    assert_eq!(p0mm + p0uu + p0um, 1.0);

    // Alpha, Beta, Weight, Intercept, pr_mm, pr_um, pr_uu
    let results = Mutex::new(Array2::<f64>::zeros((0, 7)));

    // Optimization loop
    (0..n_boot).into_par_iter().for_each(|_i| {
        // pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)
        let rng = thread_rng();
        let residual_dist = Slice::new(&residuals).unwrap();
        let residual_sample: Vec<&f64> = rng
            .sample_iter(&residual_dist)
            .take(pedigree.len_of(Axis(0)))
            .collect();
        assert!(pred_div.len() == residual_sample.len());
        let div_ops: Vec<f64> = pred_div
            .iter()
            .zip(residual_sample.iter())
            .map(|(a, b)| a + *b)
            .collect();

        let mut pedigree = pedigree.clone();
        pedigree.slice_mut(s![.., 3]).assign(&Array1::from(div_ops));

        let problem = Problem {
            pedigree,
            eqp_weight,
            eqp,
            p_mm: p0mm,
            p_um: p0um,
            p_uu: p0uu,
        };
        // Run Nelder-Mead optimization
        // Use the previous result as the initial guess, supplement with random values close-by
        let nm = NelderMead::new(vec![
            params.to_vec(),
            params.vary().to_vec(),
            params.vary().to_vec(),
            params.vary().to_vec(),
            params.vary().to_vec(),
        ]);

        let res = Executor::new(problem, nm)
            .configure(|state| {
                state
                    // .param(vec![alpha, beta, weight, intercept])
                    .max_iters(1000)
            })
            .run()
            .expect("Failed to run Nelder-Mead optimization");

        let m = Model::from_vec(&res.state.best_param.unwrap());

        let pr_mm = m.est_mm();
        let pr_um = m.est_um();
        let pr_uu = m.est_uu();
        let r = array![m.alpha, m.beta, m.weight, m.intercept, pr_mm, pr_um, pr_uu];
        results
            .lock()
            .unwrap()
            .push(Axis(0), r.view())
            .expect("Insertion failed");
        // let c = counter.fetch_add(1, Ordering::Relaxed);
        pb.inc(1);
        //  println!("Progress: {}%", ((c * 100) as f32 / (n_boot) as f32));
    });
    pb.finish();

    let results = results.into_inner().unwrap();

    plot::bootstrap(
        results.column(0).to_vec(),
        results.column(1).to_vec(),
        output_dir,
    )?;

    let raw_analysis = RawAnalysis(results);
    let analysis = raw_analysis.into();

    Ok(analysis)
}
