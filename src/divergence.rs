use std::ops::Mul;

use argmin_math::ArgminInv;
use ndarray::array;

use ndarray::Array2;

pub struct Divergence {
    pub dt1t2: Vec<f64>,
    pub puuinf_est: f64,
}
// https://numpy.org/doc/stable/reference/generated/numpy.linalg.matrix_power.html
fn matrix_power(matrix: &Array2<f64>, power: i8) -> Array2<f64> {
    if power < 0 {
        return matrix_power(&matrix.inv().unwrap(), power.abs());
    }

    if power == 0 {
        // return identity matrix
        return Array2::eye(matrix.nrows());
    }
    let mut result = matrix.clone();

    for _ in 1..power {
        result = result.dot(matrix);
    }
    result
}

pub fn divergence(
    pedigree: &Array2<f64>,
    p_mm: f64,
    _p_um: f64,
    p_uu: f64,
    alpha: f64,
    beta: f64,
    weight: f64,
) -> Divergence {
    // State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM ### Is the second field correct?
    let sv_gzero = array![p_uu, (weight) * p_mm, (1.0 - weight) * p_mm];

    // 	Defining the generation (or transition) matrix
    let genmatrix = genmatrix(alpha, beta);

    let mut dt1t2 = Vec::new();

    // 	Calculating theoretical divergence for every observed pair in 'pedigree.txt'
    for p in pedigree.rows() {
        // 			Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
        let svt0 = sv_gzero.t().dot(&matrix_power(&genmatrix, p[0] as i8));
        let svt1_mm = array![0.0, 0.0, 1.0]
            .t()
            .dot(&matrix_power(&genmatrix, (p[1] - p[0]) as i8));
        let svt2_mm = array![0.0, 0.0, 1.0]
            .t()
            .dot(&matrix_power(&genmatrix, (p[2] - p[0]) as i8));
        let svt1_um = array![0.0, 1.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, (p[1] - p[0]) as i8));
        let svt2_um = array![0.0, 1.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, (p[2] - p[0]) as i8));
        let svt1_uu = array![1.0, 0.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, (p[1] - p[0]) as i8));
        let svt2_uu = array![1.0, 0.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, (p[2] - p[0]) as i8));

        // Conditional divergences
        let dt1t2_mm = 0.5
            * (svt1_mm[0].mul(&svt2_mm[1])
                + svt1_mm[1].mul(&svt2_mm[0])
                + svt1_mm[1].mul(&svt2_mm[2])
                + svt1_mm[2].mul(&svt2_mm[1]))
            + (svt1_mm[0].mul(&svt2_mm[2]) + svt1_mm[2].mul(&svt2_mm[0]));

        let dt1t2_um = 0.5.mul(
            (svt1_um[0]).mul(&svt2_um[1])
                + svt1_um[1].mul(&svt2_um[0])
                + svt1_um[1].mul(&svt2_um[2])
                + svt1_um[2].mul(&svt2_um[1]),
        ) + (svt1_um[0].mul(&svt2_um[2]) + svt1_um[2].mul(&svt2_um[0]));

        let dt1t2_uu = 0.5.mul(
            (svt1_uu[0]).mul(&svt2_uu[1])
                + svt1_uu[1].mul(&svt2_uu[0])
                + svt1_uu[1].mul(&svt2_uu[2])
                + svt1_uu[2].mul(&svt2_uu[1]),
        ) + (svt1_uu[0].mul(&svt2_uu[2]) + svt1_uu[2].mul(&svt2_uu[0]));

        dt1t2.push(svt0[0] * (dt1t2_uu) + svt0[1] * (dt1t2_um) + svt0[2] * (dt1t2_mm));
    }

    // Pr(UU) at equilibrium given alpha and beta
    let puuinf_est = p_uu_est(alpha, beta);
    Divergence { dt1t2, puuinf_est }
}

pub fn p_uu_est(alpha: f64, beta: f64) -> f64 {
    (beta * ((1.0 - beta).powi(2) - (1.0 - alpha).powi(2) - 1.0))
        / ((alpha + beta) * ((alpha + beta - 1.0).powi(2) - 2.0))
}

pub fn genmatrix(alpha: f64, beta: f64) -> Array2<f64> {
    array![
        [
            (1.0 - alpha).powi(2),
            2.0 * (1.0 - alpha) * alpha,
            alpha.powi(2)
        ],
        [
            0.25 * (beta + 1.0 - alpha).powi(2),
            0.5 * (beta + 1.0 - alpha) * (alpha + 1.0 - beta),
            0.25 * (alpha + 1.0 - beta).powi(2)
        ],
        [
            beta.powi(2),
            2.0 * (1.0 - beta) * beta,
            (1.0 - beta).powi(2)
        ]
    ]
}

#[cfg(test)]
mod test {
    use crate::*;

    use super::*;

    macro_rules! assert_close_epsilon {
        ($x:expr, $y:expr, $d:expr) => {
            if !(($x - $y).abs() < $d) {
                panic!(
                    "assertion failed: `abs(left - right) < {}`, (left: `{}`, right: `{}`)",
                    $d, $x, $y
                );
            }
        };
    }

    macro_rules! assert_close {
        ($x:expr, $y:expr ) => {
            if !(($x - $y).abs() < 1e-6) {
                panic!(
                    "assertion failed: `abs(left - right) < {}`, (left: `{}`, right: `{}`)",
                    1e-6, $x, $y
                );
            }
        };
    }

    #[test]
    fn test_p_uu_est() {
        // Compare estimated steady state methylation to R
        assert_close!(p_uu_est(3.974271e-09, 1.519045e-07), 0.9745041);
        let _p = Model::new(1.0);
    }

    #[test]
    fn same_as_r() {
        let pedigree = Pedigree::from_file("/home/cgoeldel/epigenomics/alphabeta/pedigree.txt");
        let divergence = divergence(
            &pedigree,
            0.25,
            0.0,
            0.75,
            3.974271e-09,
            1.519045e-07,
            0.06892953,
        );

        let r = include_str!("../divergence.txt")
            .split('\n')
            .map(|s| s.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();

        assert_eq!(divergence.dt1t2.len(), r.len());
        for (i, r) in r.iter().enumerate() {
            assert!(divergence.dt1t2[i].is_normal());
            assert_close!(divergence.dt1t2[i], *r);
        }
    }
}
