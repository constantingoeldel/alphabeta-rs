use indicatif::ProgressBar;
use rayon::prelude::*;
use std::sync::Mutex;

use argmin::{core::Executor, solver::neldermead::NelderMead};
use ndarray::{array, Array1, Array2, Axis};

use crate::{
    pedigree::Pedigree,
    structs::{Model, Problem, Progress, StandardDeviations},
    *,
};

pub fn run(
    pedigree: &Pedigree,
    params: &Model,
    p0uu: f64,
    eqp: f64,
    eqp_weight: f64,
    n_boot: u64,
    pb: Option<ProgressBar>,
) -> Result<StandardDeviations, Box<dyn std::error::Error>> {
    let pb = pb.unwrap_or_else(|| Progress::new("BootModel", n_boot).0);

    let p0mm = 1.0 - p0uu;
    let p0um = 0.0;

    assert_eq!(p0mm + p0uu + p0um, 1.0);

    // Alpha, Beta, Weight, Intercept, pr_mm, pr_um, pr_uu
    let results = Mutex::new(Array2::<f64>::zeros((0, 7)));
    // let counter = AtomicU32::new(0);

    // Optimization loop
    (0..n_boot).into_par_iter().for_each(|_i| {
        let problem = Problem {
            pedigree: pedigree.clone(),
            eqp_weight,
            eqp,
            p_mm: p0mm,
            p_um: p0um,
            p_uu: p0uu,
        };
        // Run Nelder-Mead optimization
        let nm = NelderMead::new(vec![
            params.vary().to_vec(),
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

    dbg!(results.shape());

    plot::bootstrap(
        results.column(0).to_slice().unwrap(),
        results.column(1).to_slice().unwrap(),
    )?;
    // Standard deviations
    let sd_alpha = results.column(0).std(1.0);
    let sd_beta = results.column(1).std(1.0);
    // Alpha - Beta
    let sd_alpha_beta = results
        .column(1)
        .iter()
        .zip(results.column(0).iter())
        .map(|(b, a)| b / a)
        .collect::<Array1<f64>>();

    let sd_weight = results.column(2).std(1.0);
    let sd_intercept = results.column(3).std(1.0);

    let sd_pr_mm_inf = results.column(4).std(1.0);
    let sd_pr_um_inf = results.column(5).std(1.0);
    let sd_pr_uu_inf = results.column(6).std(1.0);

    // Quantiles not implemented yet

    Ok(StandardDeviations {
        alpha: sd_alpha,
        beta: sd_beta,
        alpha_beta: sd_alpha_beta.std(1.0),
        weight: sd_weight,
        intercept: sd_intercept,
        p_mm: sd_pr_mm_inf,
        p_um: sd_pr_um_inf,
        p_uu: sd_pr_uu_inf,
    })
}
