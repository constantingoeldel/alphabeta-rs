use argmin_math::ArgminInv;

use ndarray::array;

use ndarray::Array2;

use crate::alphabeta::steady_state;
use crate::pedigree::Pedigree;

#[derive(Debug)]
pub struct Divergence {
    pub dt1t2: Vec<f64>,
    pub puuinf_est: f64,
}
// https://numpy.org/doc/stable/reference/generated/numpy.linalg.matrix_power.html
pub fn matrix_power(matrix: &Array2<f64>, power: i8) -> Array2<f64> {
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
    pedigree: &Pedigree,
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
    // 	Calculating theoretical divergence for every observed pair in 'data/data/pedigree.txt'
    for p in pedigree.rows() {
        let (t0, t1, t2) = (p[0] as i8, p[1] as i8, p[2] as i8);

        // 			Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
        let svt0 = sv_gzero.t().dot(&matrix_power(&genmatrix, t0));

        let t1t0 = matrix_power(&genmatrix, t1 - t0);
        let t2t0 = matrix_power(&genmatrix, t2 - t0);

        let svt1_mm = t1t0.row(2);
        let svt2_mm = t2t0.row(2);
        let svt1_um = t1t0.row(1);
        let svt2_um = t2t0.row(1);
        let svt1_uu = t1t0.row(0);
        let svt2_uu = t2t0.row(0);

        // Conditional divergences
        let dt1t2_mm = 0.5_f64
            * (svt1_mm[0] * svt2_mm[1]
                + svt1_mm[1] * svt2_mm[0]
                + svt1_mm[1] * svt2_mm[2]
                + svt1_mm[2] * svt2_mm[1])
            + (svt1_mm[0] * svt2_mm[2] + svt1_mm[2] * svt2_mm[0]);

        let dt1t2_um = 0.5_f64
            * ((svt1_um[0]) * svt2_um[1]
                + svt1_um[1] * svt2_um[0]
                + svt1_um[1] * svt2_um[2]
                + svt1_um[2] * svt2_um[1])
            + (svt1_um[0] * svt2_um[2] + svt1_um[2] * svt2_um[0]);

        let dt1t2_uu = 0.5_f64
            * ((svt1_uu[0]) * svt2_uu[1]
                + svt1_uu[1] * svt2_uu[0]
                + svt1_uu[1] * svt2_uu[2]
                + svt1_uu[2] * svt2_uu[1])
            + (svt1_uu[0] * svt2_uu[2] + svt1_uu[2] * svt2_uu[0]);

        dt1t2.push(svt0[0] * (dt1t2_uu) + svt0[1] * (dt1t2_um) + svt0[2] * (dt1t2_mm));
    }
    // Pr(UU) at equilibrium given alpha and beta
    let puuinf_est = steady_state(alpha, beta);
    Divergence { dt1t2, puuinf_est }
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
    use crate::{structs::Model, *};

    use super::*;

    #[test]
    fn test_p_uu_est() {
        // Compare estimated steady state methylation to R
        assert_close!(steady_state(3.974271e-09, 1.519045e-07), 0.9745041);
        let _p = Model::new(1.0);
    }

    #[test]
    fn test_matrix_power_is_identity_when_power_is_zero() {
        let m = array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let m_power = matrix_power(&m, 0);
        assert_eq!(
            m_power,
            array![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
    }
    #[test]
    fn same_as_r() {
        let pedigree = Pedigree::from_file("./data/pedigree.txt");
        let divergence = divergence(
            &pedigree,
            0.25,
            0.0,
            0.75,
            3.974271e-09,
            1.519045e-07,
            0.06892953,
        );

        let r = include_str!("../data/divergence.txt")
            .split('\n')
            .map(|s| s.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();

        assert_eq!(divergence.dt1t2.len(), r.len());
        for (i, r) in r.iter().enumerate() {
            assert!(divergence.dt1t2[i].is_normal());
            assert_close!(divergence.dt1t2[i], *r);
        }
    }

    #[test]
    fn dot_product_is_the_same_as_row() {
        let genmatrix = genmatrix(0.2, 0.5);
        let dot = array![0.0, 0.0, 1.0].t().dot(&matrix_power(&genmatrix, 2));

        let row = matrix_power(&genmatrix, 2).row(2).to_owned();
        assert_eq!(dot, row);
    }
    #[test]
    fn implementation_stays_same() {
        let genmatrix = genmatrix(0.2, 0.5);
        let t2 = 2;
        let t1 = 1;
        let t0 = 0;
        let svt1_mm = array![0.0, 0.0, 1.0]
            .t()
            .dot(&matrix_power(&genmatrix, t1 - t0));
        let svt2_mm = array![0.0, 0.0, 1.0]
            .t()
            .dot(&matrix_power(&genmatrix, t2 - t0));
        let svt1_um = array![0.0, 1.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, t1 - t0));
        let svt2_um = array![0.0, 1.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, t2 - t0));
        let svt1_uu = array![1.0, 0.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, t1 - t0));
        let svt2_uu = array![1.0, 0.0, 0.0]
            .t()
            .dot(&matrix_power(&genmatrix, t2 - t0));

        let svt1_mm_row = &matrix_power(&genmatrix, t1 - t0).row(2).to_owned();
        let svt2_mm_row = &matrix_power(&genmatrix, t2 - t0).row(2).to_owned();
        let svt1_um_row = &matrix_power(&genmatrix, t1 - t0).row(1).to_owned();
        let svt2_um_row = &matrix_power(&genmatrix, t2 - t0).row(1).to_owned();
        let svt1_uu_row = &matrix_power(&genmatrix, t1 - t0).row(0).to_owned();
        let svt2_uu_row = &matrix_power(&genmatrix, t2 - t0).row(0).to_owned();

        assert_eq!(svt1_mm, svt1_mm_row);
        assert_eq!(svt2_mm, svt2_mm_row);
        assert_eq!(svt1_um, svt1_um_row);
        assert_eq!(svt2_um, svt2_um_row);
        assert_eq!(svt1_uu, svt1_uu_row);
        assert_eq!(svt2_uu, svt2_uu_row);
    }
}
