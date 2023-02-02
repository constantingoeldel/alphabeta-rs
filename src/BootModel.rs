use argmin::{core::Executor, solver::neldermead::NelderMead};

use crate::*;

pub fn run(
    pedigree: Pedigree,
    eqp: f64,
    eqp_weight: f64,
    p0uu: f64,
    params: Model,
    n_boot: i32,
) -> Result<(), Box<dyn std::error::Error>> {
    let p0mm = 1.0 - p0uu;
    let p0um = 0.0;

    assert_eq!(p0mm + p0uu + p0um, 1.0);

    let mut results = Vec::new();

    for i in 0..n_boot {
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

        let res = Executor::new(problem.clone(), nm)
            .configure(|state| {
                state
                    // .param(vec![alpha, beta, weight, intercept])
                    .max_iters(1000)
            })
            .run()?;

        let best = res.state.best_param.unwrap();

        let alpha = best[0];
        let beta = best[1];

        let predicted_mm = (alpha * ((1.0 - alpha).powi(2) - (1.0 - beta).powi(2) - 1.0))
            / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));
        let predicted_um = (4.0 * alpha * beta * (alpha + beta - 2.0))
            / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));
        let predicted_uu = (beta * ((1.0 - beta).powi(2) - (1.0 - alpha).powi(2) - 1.0))
            / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0));

        results.push(Res {
            alpha,
            beta,
            weight: best[2],
            intercept: best[3],
            predicted_mm,
            predicted_um,
            predicted_uu,
        });

        println!("Progress: {}%", i / n_boot * 100);
    }
    Ok(())
}
