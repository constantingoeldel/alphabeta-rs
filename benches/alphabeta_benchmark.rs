use alphabeta::{
    alphabeta::steady_state,
    divergence::{divergence, genmatrix, matrix_power},
    pedigree::Pedigree,
    structs::Model,
};
use argmin_math::ArgminMul;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ndarray::array;

fn criterion_benchmark(c: &mut Criterion) {
    let p_uu = black_box(0.5);
    let p_um = black_box(0.0);
    let p_mm = black_box(0.5);
    let model = black_box(Model::default());

    let pedigree = black_box(Pedigree::from_file("./data/pedigree.txt"));
    let sv_gzero = array![p_uu, (model.weight) * p_mm, (1.0 - model.weight) * p_mm];
    let p = pedigree.row(0);

    c.bench_function("genmatrix", |b| {
        b.iter(|| genmatrix(model.alpha, model.beta))
    });
    let genmatrix = black_box(genmatrix(model.alpha, model.beta));
    c.bench_function("matrix square", |b| b.iter(|| genmatrix.mul(&genmatrix)));
    c.bench_function("matrix cube", |b| {
        b.iter(|| genmatrix.mul(&genmatrix.mul(&genmatrix)))
    });

    c.bench_function("matrix power", |b| b.iter(|| matrix_power(&genmatrix, 3)));

    c.bench_function("divergence hot loop", |b| {
        b.iter(|| {
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
            svt0[0] * (dt1t2_uu) + svt0[1] * (dt1t2_um) + svt0[2] * (dt1t2_mm)
        })
    });

    c.bench_function("divergence", |b| {
        b.iter(|| {
            divergence(
                &pedigree,
                p_mm,
                p_um,
                p_uu,
                model.alpha,
                model.beta,
                model.weight,
            )
        })
    });

    c.bench_function("p_inf", |b| {
        b.iter(|| steady_state(black_box(model.alpha), black_box(model.beta)))
    });
    c.bench_function("integer conversion", |b| {
        b.iter(|| {
            let (t0, t1, t2) = (p[0] as i8, p[1] as i8, p[2] as i8);
            (t0, t1, t2)
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
