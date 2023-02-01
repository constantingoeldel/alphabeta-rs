pub type Pedigree = Array2<f64>;

use argmin::{
    core::{CostFunction, Executor},
    solver::neldermead::NelderMead,
};
use ndarray::Array2;
use rand::{distributions::Uniform, thread_rng, Rng};

use crate::divergence::divergence;

#[derive(Clone, Debug)]
struct Problem {
    pedigree: Pedigree,
    eqp_weight: f64,
    eqp: f64,
    p_mm: f64,
    p_um: f64,
    p_uu: f64,
}
#[derive(Clone, Debug)]
pub struct Params {
    pub alpha: f64,
    pub beta: f64,
    weight: f64,
    intercept: f64,
}

impl Params {
    pub fn new(max_divergence: f64) -> Self {
        let mut rng = thread_rng();
        let alpha = rng.sample(Uniform::new(-9.0, -2.0));
        let beta = rng.sample(Uniform::new(-9.0, -2.0));
        let weight = rng.sample(Uniform::new(0.0, 0.1));
        let intercept = rng.sample(Uniform::new(0.0, max_divergence));
        Params {
            alpha,
            beta,
            weight,
            intercept,
        }
    }
    fn to_vec(&self) -> Vec<f64> {
        vec![self.alpha, self.beta, self.weight, self.intercept]
    }
    fn from_vec(v: &[f64]) -> Self {
        Params {
            alpha: v[0],
            beta: v[1],
            weight: v[2],
            intercept: v[3],
        }
    }
}

impl CostFunction for Problem {
    type Output = f64;
    type Param = Vec<f64>;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let p = Params::from_vec(p);
        let divergence = divergence(
            &self.pedigree,
            self.p_mm,
            self.p_um,
            self.p_uu,
            p.alpha,
            p.beta,
            p.weight,
        );

        let mut square_sum = 0.0;

        for (div, ped) in divergence.dt1t2.iter().zip(self.pedigree.column(3)) {
            square_sum += (ped - p.intercept - div).powi(2)
                + self.eqp_weight * self.pedigree.len() as f64 * (div - self.eqp).powi(2);
        }

        println!("{square_sum:?}");
        Ok(square_sum)

        // Ok(p.iter().fold(0.0, |acc, x| acc + x.powi(2)))
    }
}

pub fn run(
    pedigree: Pedigree,
    p0uu: f64,
    eqp: f64,
    eqp_weight: f64,
    n_starts: i32,
) -> Result<(), Box<dyn std::error::Error>> {
    let p0mm = 1.0 - p0uu;
    let p0um = 0.0;
    let max_divergence = *pedigree
        .column(3)
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    assert_eq!(p0mm + p0uu + p0um, 1.0);

    struct Result {
        alpha: f64,
        beta: f64,
        weight: f64,
        intercept: f64,
        predicted_mm: f64,
        predicted_um: f64,
        predicted_uu: f64,
    }

    let mut results = Vec::new();

    // Optimization loop
    for i in 0..n_starts {
        // Draw random starting values

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
            Params::new(max_divergence).to_vec(),
            Params::new(max_divergence).to_vec(),
            Params::new(max_divergence).to_vec(),
            Params::new(max_divergence).to_vec(),
            Params::new(max_divergence).to_vec(),
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

        results.push(Result {
            alpha,
            beta,
            weight: best[2],
            intercept: best[3],
            predicted_mm,
            predicted_um,
            predicted_uu,
        });

        println!("Progress: {}%", i / n_starts * 100);
    }
    // Calculating the least squares error for all results and selecting the best one
    results.sort_by(|a, b| {
        let pedigree = pedigree.clone();
        let divergence_a = divergence(&pedigree, p0mm, p0um, p0uu, a.alpha, a.beta, a.weight);
        let divergence_b = divergence(&pedigree, p0mm, p0um, p0uu, b.alpha, b.beta, b.weight);

        let lse_a = divergence_a
            .dt1t2
            .iter()
            .zip(pedigree.column(3))
            .map(|(div, ped)| (ped - a.intercept - div).powi(2))
            .sum::<f64>();
        let lse_b = divergence_b
            .dt1t2
            .iter()
            .zip(pedigree.column(3))
            .map(|(div, ped)| (ped - b.intercept - div).powi(2))
            .sum::<f64>();
        println!("{lse_a} {lse_b}");
        lse_a.partial_cmp(&lse_b).unwrap() // Sort ascending
    });

    // Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
    // Caution: Calculating predicted divergence based on lowest LSQ model: check the biology!", "\n")

    let best = &results[0];
    let divergence = divergence(
        &pedigree,
        p0mm,
        p0um,
        p0uu,
        best.alpha,
        best.beta,
        best.weight,
    );

    let mut delta_t = Vec::new();
    for row in pedigree.rows() {
        delta_t.push(row[2] - row[3] - 2.0 * row[1])
    }

    let mut predicted_divergence = Vec::new();

    for (i, _row) in pedigree.rows().into_iter().enumerate() {
        predicted_divergence.push(
            best.intercept + divergence.dt1t2[i], /*+ delta_t[i] * best.alpha*/
        );
    }

    let mut residual = Vec::new();

    for (i, row) in pedigree.rows().into_iter().enumerate() {
        residual.push(row[3] - predicted_divergence[i]);
    }

    // Generating theoretical fit

    // Not needed for now

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::pedigree;

    use super::*;

    #[test]
    fn it_runs() {
        let pedigree =
            pedigree::read_pedigree_from_file("/home/cgoeldel/epigenomics/alphabeta/pedigree.txt");

        let result = run(pedigree, 0.75, 0.5, 0.7, 2);
        println!("{result:?}");
        assert_eq!(2 + 2, 4);
    }
}
