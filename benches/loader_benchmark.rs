use alphabeta::loader::Methylome;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn loader_benchmark(c: &mut Criterion) {
    c.bench_function("genmatrix", |b| {
        b.iter(|| {
            let meth = Methylome::from_nodelist_and_edgelist(
                "/mnt/fast/epigenomics/MA3_new_total_original_methylome/nodelist.txt",
                "/mnt/fast/epigenomics/MA3_new_total_original_methylome/edgelist.txt",
            )
            .unwrap();
            for node in meth.nodes() {
                match &node.methylome {
                    Some(m) => assert_eq!(m.len(), 410000),
                    None => panic!("No methylome found"),
                }
            }
        })
    });
}

criterion_group!(benches, loader_benchmark);
criterion_main!(benches);
