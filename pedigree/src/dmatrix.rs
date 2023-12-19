use petgraph::{adj::NodeIndex, algo::astar, prelude::UnGraph, visit::Dfs, Direction::Incoming};

use ndarray::{array, Array2, Axis};
use polars::prelude::*;

use crate::{DivergenceBetweenSamples, Node, Pedigree, Return};

#[derive(Debug)]
pub struct DMatrix(Array2<f64>);

impl DMatrix {
    pub fn from(pedigree: &Pedigree, posterior_max: f64) -> Return<Self> {
        println!();
        println!();
        println!("Calculating divergences between all pairs of nodes...");
        println!();
        let n_nodes = pedigree.df_count();
        let mut divergences = Array2::<f64>::zeros((n_nodes, n_nodes));

        let pb = progress_bars::Progress::new(
            "Calculating Divergences ",
            binomial_coefficient(n_nodes, 2),
        );

        // Go over all pairs of nodes, excluding self-pairs
        for (i, first) in pedigree.dfs().enumerate() {
            for (j, second) in pedigree.dfs().skip(i + 1).enumerate() {
                pb.inc(1);

                // if first.height() != second.height() {
                //     println!(
                //         "Lengths do not match, all bets are off: {} vs {}",
                //         first.height(),
                //         second.height()
                //     );
                //     divergences[[i, j]] = 0.0;
                //     continue;
                // }

                let mut divergence = 0;
                let mut compared_sites = 0;

                let categorical_to_numeric = |s: &Series| {
                    s.iter()
                        .map(|v| match v {
                            AnyValue::Utf8("U") => 0,
                            AnyValue::Utf8("I") => 1,
                            AnyValue::Utf8("M") => 2,
                            _ => panic!("Unknown status"),
                        })
                        .collect::<Series>()
                };

                let div = categorical_to_numeric(first.column("status")?)
                    - categorical_to_numeric(second.column("status")?);

                //     let df_1 = first
                //         .select([col("status"), col("posteriormax")])
                //         .collect()
                //         .unwrap();

                //     let df_2 = second
                //         .select([col("status"), col("posteriormax")])
                //         .collect()
                //         .unwrap();

                //    for (f,s) in df_1.column("status").unwrap().iter().zip(df_2.column("status").unwrap().iter())

                // Go over all sites in the first sample
                // IMPORTANT: It is assumed that the same sites are included in the datasets and that the sites are sorted by position
                // for (k, f) in first.sites.iter().enumerate() {
                //     let s = second
                //         .sites
                //         .get(k)
                //         .expect("Partner methylation site must exists");

                //     // This is the correct place for posterior_max filtering
                //     // Well - almost correct, because we should also filter out sites from other nodes at the same position if posterior_max is insufficient in just one sample
                //     if f.posteriormax < posterior_max || s.posteriormax < posterior_max {
                //         continue;
                //     }

                //     divergence += f.status_numeric().abs_diff(s.status_numeric());
                //     compared_sites += 1;
                // }

                // let divergence = divergence as f64 / (2.0 * compared_sites as f64);
                let sum = abs(&div).unwrap().sum().unwrap();
                dbg!(i, j, sum);
                divergences[[i, j]] = sum; // There is only one row in the dataframe
            }
        }
        pb.finish();
        Ok(DMatrix(divergences))
    }
    /// Convert graph of divergences to pedigree
    pub fn convert(&self, pedigree: &Pedigree) -> DivergenceBetweenSamples {
        // let e = edges
        //     .iter()
        //     .map(|e| {
        //         (
        //             e.from.id,
        //             e.to.id,
        //             e.from.generation.abs_diff(e.to.generation) as usize,
        //             // self.0.get((e.from.id, e.to.id)).unwrap(),
        //         )
        //     })
        //     .collect::<Vec<(usize, usize, usize)>>();

        // let graph = UnGraph::<usize, usize, usize>::from_edges(e);

        // Iterate over all pairs of nodes
        let mut divergence = Array2::<f64>::default((0, 4));
        for (i, source) in pedigree.nodes_with_sites().enumerate() {
            for (j, target) in pedigree.nodes_with_sites().skip(i + 1).enumerate() {
                // find latest common ancestor in directed acyclic graph
                fn lcm<'g>(pedigree: &'g Pedigree, n: &'g Node, m: &'g Node) -> Option<&'g Node> {
                    let parent = |n: &Node| -> Option<&Node> {
                        pedigree
                            .node_weight(pedigree.neighbors_directed(n.id.into(), Incoming).next()?)
                    };

                    match n.generation.cmp(&m.generation) {
                        std::cmp::Ordering::Less => lcm(pedigree, n, parent(m)?),
                        std::cmp::Ordering::Greater => lcm(pedigree, parent(n)?, m),
                        std::cmp::Ordering::Equal => {
                            if n.id == m.id {
                                Some(n)
                            } else {
                                lcm(pedigree, parent(n)?, parent(m)?)
                            }
                        }
                    }
                }

                let t0 = lcm(pedigree, source, target).unwrap().generation as f64;

                let t1 = source.generation as f64;
                let t2 = target.generation as f64;

                let div = self.0.get((i, j)).unwrap().to_owned();
                let generational_distance = t2 + t1 - 2.0 * t0;
                dbg!(t0, t1, t2, div, generational_distance);
                divergence
                    .push(Axis(0), array![t0, t1, t2, div].view())
                    .expect("Could not insert row into divergence list");
            }
        }

        divergence
    }
}

fn binomial_coefficient(n: usize, k: usize) -> usize {
    let mut result = 1;
    let k = k.min(n - k); // take advantage of symmetry
    for i in 0..k {
        result *= n - i;
        result /= i + 1;
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_binomial_coefficient() {
        assert_eq!(binomial_coefficient(5, 0), 1);
        assert_eq!(binomial_coefficient(5, 1), 5);
        assert_eq!(binomial_coefficient(5, 2), 10);
        assert_eq!(binomial_coefficient(5, 3), 10);
        assert_eq!(binomial_coefficient(5, 4), 5);
        assert_eq!(binomial_coefficient(5, 5), 1);
    }
}
