use std::{fmt::Display, fs::File, io::Write, path::Path};

use argmin::core::CostFunction;

use ndarray::Array1;
use rand::{distributions::Uniform, thread_rng, Rng};

use crate::divergence::{divergence, Divergence};
use pedigree::{DivergenceBetweenSamples, Pedigree};

#[derive(Clone, Debug)]
pub struct Problem {
    pub divergence: DivergenceBetweenSamples,
    pub eqp_weight: f64,
    pub eqp: f64,
    pub p_mm: f64,
    pub p_um: f64,
    pub p_uu: f64,
}
#[derive(Clone, Debug)]
pub struct Model {
    pub alpha: f64,
    pub beta: f64,
    pub weight: f64,
    pub intercept: f64,
}

pub type PredictedDivergence = Vec<f64>;
pub type Residuals = Vec<f64>;

impl Display for Model {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Model:\n\tAlpha: {}\n\tBeta: {}\n\tWeight: {}\n\tIntercept: {}",
            self.alpha, self.beta, self.weight, self.intercept
        )
    }
}

/// For testing purposes, we want to be able to create a model with known parameters.
impl Default for Model {
    fn default() -> Self {
        Model {
            alpha: 0.0001179555,
            beta: 0.0001180614,
            weight: 0.03693534,
            intercept: 0.003023981,
        }
    }
}

impl Model {
    pub fn new(max_divergence: f64) -> Self {
        let mut max = max_divergence;
        if max_divergence <= 0.0 {
            println!("Sample has a maximum divergence of zero! Check your data");
            max = 0.1;
        }

        let mut rng = thread_rng();
        let alpha = 10.0_f64.powf(rng.sample(Uniform::new(-9.0, -2.0)));
        let beta = 10.0_f64.powf(rng.sample(Uniform::new(-9.0, -2.0)));
        let weight = rng.sample(Uniform::new(0.0, 0.1));
        let intercept = rng.sample(Uniform::new(0.0, max));
        Model {
            alpha,
            beta,
            weight,
            intercept,
        }
    }
    /// Returns a new model with parameters that are randomly varied by up to 5% of their original value.
    ///
    /// I made sure to check that only positive, non-zero floats can be passed, but this is a nicer way to handle errors as the panic is not well-readable
    pub fn vary(&self) -> Self {
        let mut rng = thread_rng();
        const VARIANCE: f64 = 0.1; // 10% variance, somewhat arbitrarily chosen, but even 0.5 does not have a huge effect on the results.

        fn var(n: f64) -> Uniform<f64> {
            if n == 0.0 {
                println!("Warning: One of the parameters you passed was zero");
                return var(0.1);
            }

            let low = n - n.abs() * VARIANCE;
            let high = n + n.abs() * VARIANCE;

            if low >= high {
                println!("Wrong order: Low >= high ({low} >= {high})");
                Uniform::new(high, low)
            } else {
                Uniform::new(low, high)
            }
        }

        Model::from_vec(
            &self
                .to_vec()
                .iter()
                .map(|s| rng.sample(var(*s)))
                .collect::<Vec<f64>>(),
        )
    }

    pub fn to_vec(&self) -> Vec<f64> {
        vec![self.alpha, self.beta, self.weight, self.intercept]
    }
    pub fn from_vec(v: &[f64]) -> Self {
        Model {
            alpha: v[0],
            beta: v[1],
            weight: v[2],
            intercept: v[3],
        }
    }

    pub fn to_array(&self) -> Array1<f64> {
        Array1::from(self.to_vec())
    }

    pub fn est_mm(&self) -> f64 {
        (self.alpha * ((1.0 - self.alpha).powi(2) - (1.0 - self.beta).powi(2) - 1.0))
            / ((self.alpha + self.beta) * ((self.alpha + self.beta - 1.0).powi(2) - 2.0))
    }

    pub fn est_um(&self) -> f64 {
        (4.0 * self.alpha * self.beta * (self.alpha + self.beta - 2.0))
            / ((self.alpha + self.beta) * ((self.alpha + self.beta - 1.0).powi(2) - 2.0))
    }

    pub fn est_uu(&self) -> f64 {
        (self.beta * ((1.0 - self.beta).powi(2) - (1.0 - self.alpha).powi(2) - 1.0))
            / ((self.alpha + self.beta) * ((self.alpha + self.beta - 1.0).powi(2) - 2.0))
    }
    pub fn to_file(&self, path: &Path) -> std::io::Result<()> {
        println!("Writing model to file: {}", path.display());
        let mut file = File::create(path).unwrap();
        let content = format!(
            "Alpha {}\nBeta {}\nWeight {}\n Intercept {}\n",
            self.alpha, self.beta, self.weight, self.intercept
        );

        file.write_all(content.as_bytes())
    }
}

impl Default for Problem {
    fn default() -> Self {
        let divergence = Pedigree::divergence_from_file("./data/pedigree.txt");
        let p_uu = 0.75;
        let p_mm = 1.0 - p_uu;
        let p_um = 0.0;
        let eqp = 0.5;
        let eqp_weight = 0.7;
        Problem {
            divergence,
            p_mm,
            p_um,
            p_uu,
            eqp_weight,
            eqp,
        }
    }
}

impl CostFunction for Problem {
    type Output = f64;
    type Param = Vec<f64>;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let p = Model::from_vec(p);
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
                + self.eqp_weight
                    * self.pedigree.nrows() as f64
                    * (divergence.p_uu - self.eqp).powi(2);
        }

        Ok(square_sum)
    }
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use super::*;

    fn cost_function_tester<C: CostFunction>(c: C)
    where
        <C as argmin::core::CostFunction>::Param: From<Vec<f64>>,
        <C as argmin::core::CostFunction>::Output: std::cmp::PartialEq<f64> + Debug,
    {
        let param = Model::default().to_vec().into();
        let result = C::cost(&c, &param);
        let result = result.unwrap();
        assert_eq!(result, 0.0006700888539608879);
    }

    #[test]
    fn test_cost_function() {
        let p = Problem::default();
        cost_function_tester::<Problem>(p);
    }
}
