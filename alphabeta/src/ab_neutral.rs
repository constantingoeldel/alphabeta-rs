use std::sync::Mutex;

use crate::divergence::Divergence;
use crate::optimizer::{Model, PredictedDivergence, Problem, Residuals};
use crate::{config, Error, Return};
use argmin::{core::Executor, solver::neldermead::NelderMead};
use pedigree::{DivergenceBetweenSamples, Pedigree};

use rayon::prelude::*;

pub fn run(
    divergence: &DivergenceBetweenSamples,
    p0uu: f64,
    eqp: f64,
    eqp_weight: f64,
) -> Return<(Model, PredictedDivergence, Residuals)> {
    let pb = progress_bars::Progress::new("ABneutral", config::get().iterations);
    let n_starts = config::get().iterations;
    let p0mm = 1.0 - p0uu;
    let p0um = 0.0;
    let max_divergence = *divergence
        .column(3)
        .iter()
        .max_by(|a, b| a.partial_cmp(b).expect("There is a NaN value in the divergence column"))
        .expect("Maximum divergence could not be calculated, most likely because there are no valid rows in the pedigree");

    assert_eq!(p0mm + p0uu + p0um, 1.0);

    let results: Mutex<Vec<Model>> = Mutex::new(Vec::new());
    // let counter = AtomicU32::new(0);

    // Optimization loop
    let res: Result<Vec<_>, Error> = (0..n_starts)
        .into_par_iter()
        .map(|_| -> Return<()> {
            // Draw random starting values

            let problem = Problem {
                divergence: divergence.clone(),
                eqp_weight,
                eqp,
                p_mm: p0mm,
                p_um: p0um,
                p_uu: p0uu,
            };
            // Run Nelder-Mead optimization
            let nm = NelderMead::new(vec![
                Model::new(max_divergence).to_vec(),
                Model::new(max_divergence).to_vec(),
                Model::new(max_divergence).to_vec(),
                Model::new(max_divergence).to_vec(),
                Model::new(max_divergence).to_vec(),
            ]);

            let res = Executor::new(problem, nm)
                .configure(|state| {
                    state
                        // .param(vec![alpha, beta, weight, intercept])
                        .max_iters(10000)
                })
                .run()?;

            let m = Model::from_vec(&res.state.best_param.expect(
                "Best parameters must have been found if the optimization ran successfully",
            ));

            // let predicted_mm = (m.alpha * ((1.0 - m.alpha).powi(2) - (1.0 - m.beta).powi(2) - 1.0))
            //     / ((m.alpha + m.beta) * ((m.alpha + m.beta - 1.0).powi(2) - 2.0));
            // let predicted_um = (4.0 * m.alpha * m.beta * (m.alpha + m.beta - 2.0))
            //     / ((m.alpha + m.beta) * ((m.alpha + m.beta - 1.0).powi(2) - 2.0));
            // let predicted_uu = (m.beta * ((1.0 - m.beta).powi(2) - (1.0 - m.alpha).powi(2) - 1.0))
            //     / ((m.alpha + m.beta) * ((m.alpha + m.beta - 1.0).powi(2) - 2.0));
            results.lock().unwrap().push(m);
            // let c = counter.fetch_add(1, Ordering::SeqCst);
            pb.inc(1);
            //  println!("Progress: {}%", ((c * 100) as f32 / (n_starts) as f32));
            Ok(())
        })
        .collect();

    res?;
    pb.finish();

    let mut results = results.into_inner().unwrap();
    // Calculating the least squares error for all results and selecting the best one
    results.sort_by(|a, b| {
        let divergence = divergence.clone();
        let divergence_a =
            Divergence::calc(&divergence, p0mm, p0um, p0uu, a.alpha, a.beta, a.weight);
        let divergence_b =
            Divergence::calc(&divergence, p0mm, p0um, p0uu, b.alpha, b.beta, b.weight);

        let lse_a = divergence_a
            .dt1t2
            .iter()
            .zip(divergence.column(3))
            .map(|(div, ped)| (ped - a.intercept - div).powi(2))
            .sum::<f64>();
        let lse_b = divergence_b
            .dt1t2
            .iter()
            .zip(divergence.column(3))
            .map(|(div, ped)| (ped - b.intercept - div).powi(2))
            .sum::<f64>();
        lse_a
            .partial_cmp(&lse_b)
            .expect("There is a NaN value in the LSE comparison") // Sort ascending
    });

    // Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
    // Caution: Calculating predicted divergence based on lowest LSQ model: check the biology!", "\n")

    let best: &Model = &results[0];

    let div = Divergence::calc(
        divergence,
        p0mm,
        p0um,
        p0uu,
        best.alpha,
        best.beta,
        best.weight,
    );

    let mut delta_t = Vec::new();
    for row in divergence.rows() {
        delta_t.push(row[2] - row[3] - 2.0 * row[1])
    }

    let mut predicted_divergence = Vec::new();

    for (i, _row) in divergence.rows().into_iter().enumerate() {
        predicted_divergence.push(
            best.intercept + div.dt1t2[i], /*+ delta_t[i] * best.alpha*/
        );
    }

    let mut residuals = Vec::new();

    for (i, row) in divergence.rows().into_iter().enumerate() {
        residuals.push(row[3] - predicted_divergence[i]);
    }

    // Generating theoretical fit

    // Not needed for now

    Ok((best.to_owned(), predicted_divergence, residuals))
}

#[cfg(test)]
mod tests {
    // use super::*;
    // Takes quite a while to run
    //     #[test]
    //     fn it_calculates_model() {
    //         let pedigree = Pedigree::from_file("./data/pedigree.txt");

    //         let result = ABneutral::run(&pedigree, 0.75, 0.5, 0.7, 1).expect("Model failed");

    //         let r = Model::default();

    //         assert_close!(result.alpha, r.alpha);
    //         assert_close!(result.beta, r.beta);
    //         assert_close!(result.weight, r.weight);
    //         assert_close!(result.intercept, r.intercept);
    //     }
}
