use std::ops::Deref;

use argmin::core::CostFunction;
use ndarray::{Array2, ArrayView};
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

pub struct Res {
    pub alpha: f64,
    pub beta: f64,
    pub weight: f64,
    pub intercept: f64,
    pub predicted_mm: f64,
    pub predicted_um: f64,
    pub predicted_uu: f64,
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

impl Model {
    pub fn new(max_divergence: f64) -> Self {
        let mut rng = thread_rng();
        let alpha = rng.sample(Uniform::new(-9.0, -2.0));
        let beta = rng.sample(Uniform::new(-9.0, -2.0));
        let weight = rng.sample(Uniform::new(0.0, 0.1));
        let intercept = rng.sample(Uniform::new(0.0, max_divergence));
        Model {
            alpha,
            beta,
            weight,
            intercept,
        }
    }

    pub fn vary(&self) -> Self {
        let mut rng = thread_rng();
        Model {
            alpha: rng.sample(Uniform::new(self.alpha * 0.5, self.alpha * 1.5)),
            beta: rng.sample(Uniform::new(self.beta * 0.5, self.beta * 1.5)),
            weight: rng.sample(Uniform::new(self.weight * 0.5, self.weight * 1.5)),
            intercept: rng.sample(Uniform::new(self.intercept * 0.5, self.intercept * 1.5)),
        }
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
                + self.eqp_weight * self.pedigree.len() as f64 * (div - self.eqp).powi(2);
        }

        println!("{square_sum:?}");
        Ok(square_sum)

        // Ok(p.iter().fold(0.0, |acc, x| acc + x.powi(2)))
    }
}
