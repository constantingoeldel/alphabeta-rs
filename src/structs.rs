use std::{fmt::Display, ops::Deref};

use argmin::core::CostFunction;
use ndarray::{Array1, Array2, ArrayView};
use rand::{distributions::Uniform, thread_rng, Rng};

use crate::divergence::divergence;

#[derive(Clone, Debug)]
pub struct Pedigree(Array2<f64>);

#[derive(Clone, Debug)]
pub struct Problem {
    pub pedigree: Pedigree,
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

pub struct StandardDeviations {
    pub alpha: f64,
    pub beta: f64,
    pub alpha_beta: f64,
    pub weight: f64,
    pub intercept: f64,
    pub p_mm: f64,
    pub p_um: f64,
    pub p_uu: f64,
}

impl Display for Model {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Model:\t Alpha: {}\tBeta: {}\tWeight: {}\tIntercept: {}",
            self.alpha, self.beta, self.weight, self.intercept
        )
    }
}

impl Display for StandardDeviations {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Standard Deviations:\t Alpha: {}\tBeta: {}\tBeta/Alpha: {}\tWeight: {}\tIntercept: {}\tPr_mm: {}\tPr_um: {}\tPr_uu: {}",
            self.alpha,
            self.beta,
            self.alpha_beta,
            self.weight,
            self.intercept,
            self.p_mm,
            self.p_um,
            self.p_uu,
        )
    }
}

impl Pedigree {
    pub fn from_file(filename: &str) -> Self {
        let mut pedigree = Array2::<f64>::zeros((0, 4));
        let file = std::fs::read_to_string(filename).unwrap();
        file.split('\n').skip(1).for_each(|line| {
            let mut entries = line.split(' ');
            let row = [
                entries.next().unwrap().parse::<f64>().unwrap(),
                entries.next().unwrap().parse::<f64>().unwrap(),
                entries.next().unwrap().parse::<f64>().unwrap(),
                entries.next().unwrap().parse::<f64>().unwrap(),
            ];
            pedigree.push_row(ArrayView::from(&row)).unwrap();
        });
        Pedigree(pedigree)
    }
}

impl Deref for Pedigree {
    type Target = Array2<f64>;

    fn deref(&self) -> &Self::Target {
        &self.0
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
        let mut rng = thread_rng();
        let alpha = 10.0_f64.powf(rng.sample(Uniform::new(-9.0, -2.0)));
        let beta = 10.0_f64.powf(rng.sample(Uniform::new(-9.0, -2.0)));
        let weight = rng.sample(Uniform::new(0.0, 0.1));
        let intercept = rng.sample(Uniform::new(0.0, max_divergence));
        Model {
            alpha,
            beta,
            weight,
            intercept,
        }
    }
    /// Returns a new model with parameters that are randomly varied by up to 5% of their original value.
    pub fn vary(&self) -> Self {
        let mut rng = thread_rng();
        Model::from_vec(
            &self
                .to_vec()
                .iter()
                .map(|s| rng.sample(Uniform::new(s - s.abs() * 0.05, s + s.abs() * 0.05)))
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
}

impl Default for Problem {
    fn default() -> Self {
        let pedigree = Pedigree::from_file("./pedigree.txt");
        let p_uu = 0.75;
        let p_mm = 1.0 - p_uu;
        let p_um = 0.0;
        let eqp = 0.5;
        let eqp_weight = 0.7;
        Problem {
            pedigree,
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
        //   println!("{:?}", p);
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
        // println!("{:#?}", divergence);

        let mut square_sum = 0.0;

        for (div, ped) in divergence.dt1t2.iter().zip(self.pedigree.column(3)) {
            square_sum += (ped - p.intercept - div).powi(2)
                + self.eqp_weight
                    * self.pedigree.nrows() as f64
                    * (divergence.puuinf_est - self.eqp).powi(2);

            //    println!("{:?} {:?} {:?} {:?}", ped, p.intercept, div, square_sum)
        }

        println!("{square_sum:?}");
        Ok(square_sum)

        // Ok(p.iter().fold(0.0, |acc, x| acc + x.powi(2)))
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
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_cost_function() {
        let p = Problem::default();
        cost_function_tester::<Problem>(p);
    }
}
