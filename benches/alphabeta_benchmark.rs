use alphabeta::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ndarray::{array, Array2};

fn criterion_benchmark(c: &mut Criterion) {
    let p_uu = black_box(0.5);
    let p_um = black_box(0.0);
    let p_mm = black_box(0.5);
    let model = black_box(Model::default());

    let pedigree = black_box(Pedigree::from_file("./pedigree.txt"));
    let sv_gzero = array![p_uu, (model.weight) * p_mm, (1.0 - model.weight) * p_mm];
    let p = pedigree.row(0);

    c.bench_function("genmatrix", |b| {
        b.iter(|| genmatrix(model.alpha, model.beta))
    });

    let genmatrix = genmatrix(model.alpha, model.beta);
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
    c.bench_function("matrix power", |b| {
        b.iter(|| {
            let mut matrix = Array2::<f64>::zeros((3, 3));
            matrix.fill(0.5);
            matrix_power(&matrix, 5)
        })
    });
    c.bench_function("state vectors", |b| {
        b.iter(|| {
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
            (svt0, svt1_mm, svt2_mm, svt1_um, svt2_um, svt1_uu, svt2_uu)
        })
    });

    c.bench_function("conditional divergences", |b| {
        b.iter(|| {
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

            //  dt1t2.push(svt0[0] * (dt1t2_uu) + svt0[1] * (dt1t2_um) + svt0[2] * (dt1t2_mm));

            (dt1t2_mm, dt1t2_um, dt1t2_uu)
        })
    });

    c.bench_function("p_inf", |b| {
        b.iter(|| p_uu_est(black_box(model.alpha), black_box(model.beta)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
