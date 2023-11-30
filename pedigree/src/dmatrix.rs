use petgraph::{algo::astar, prelude::UnGraph};

use ndarray::{array, Array2, Axis};

use crate::{Edge, Node, Pedigree};

#[derive(Debug)]
pub struct DMatrix(Array2<f64>);

impl DMatrix {
    pub fn from(nodes: &Vec<Node>, posterior_max: f64) -> Self {
        println!();
        println!();
        println!("Calculating divergences between all pairs of nodes...");
        println!();
        let mut divergences = Array2::<f64>::zeros((nodes.len(), nodes.len()));

        let pb = progress_bars::Progress::new(
            "Calculating Divergences ",
            binomial_coefficient(nodes.len(), 2),
        );

        // Go over all pairs of nodes, excluding self-pairs
        for (i, first) in nodes.iter().enumerate() {
            for (j, second) in nodes.iter().skip(i + 1).enumerate() {
                pb.inc(1);
                assert!(first.sites.as_ref().is_some());
                assert!(second.sites.as_ref().is_some());
                if first.sites.as_ref().unwrap().len() != second.sites.as_ref().unwrap().len() {
                    println!(
                        "Lengths do not match, all bets are off: {} vs {}",
                        first.sites.as_ref().unwrap().len(),
                        second.sites.as_ref().unwrap().len()
                    );
                    divergences[[i, j]] = 0.0;
                    continue;
                }
                assert_eq!(
                    first.sites.as_ref().unwrap().len(),
                    second.sites.as_ref().unwrap().len()
                );

                let mut divergence = 0;
                let mut compared_sites = 0;

                // Go over all sites in the first sample
                // IMPORTANT: It is assumed that the same sites are included in the datasets and that the sites are sorted by position
                for (k, f) in first.sites.as_ref().unwrap().iter().enumerate() {
                    let s = second
                        .sites
                        .as_ref()
                        .unwrap()
                        .get(k)
                        .expect("Partner methylation site must exists");

                    // TODO: remove, is duplicate
                    if f.posteriormax < posterior_max || s.posteriormax < posterior_max {
                        continue;
                    }

                    divergence += f.status_numeric().abs_diff(s.status_numeric());
                    compared_sites += 1;
                }

                let divergence = divergence as f64 / (2.0 * compared_sites as f64);
                divergences[[i, j]] = divergence;
            }
        }
        pb.finish();
        DMatrix(divergences)
    }
    /// Convert graph of divergences to pedigree
    pub fn convert(&self, nodes: &[Node], edges: &[Edge]) -> Pedigree {
        let e = edges
            .iter()
            .map(|e| {
                (
                    e.from.id,
                    e.to.id,
                    e.from.generation.abs_diff(e.to.generation) as usize,
                    // self.0.get((e.from.id, e.to.id)).unwrap(),
                )
            })
            .collect::<Vec<(usize, usize, usize)>>();

        let graph = UnGraph::<usize, usize, usize>::from_edges(e);

        let mut pedigree = Pedigree(Array2::<f64>::default((0, 4)));

        for (i, source) in nodes.iter().enumerate() {
            for (j, target) in nodes.iter().skip(i + 1).enumerate() {
                if source.id == target.id {
                    // Explicitly skip self-pairs
                    continue;
                }

                let path = astar(
                    &graph,
                    source.id.into(),
                    |finish| finish == target.id.into(),
                    |e| *e.weight(),
                    |_| 0,
                );

                match path {
                    None => continue,
                    Some(path) => {
                        let distance = path.0;
                        let visited = path.1;

                        let t0: f64 = visited
                            .iter()
                            .map(|n| {
                                edges
                                    .iter()
                                    .find_map(|e| {
                                        if e.from.id == n.index() {
                                            return Some(e.from.generation);
                                        }
                                        if e.to.id == n.index() {
                                            return Some(e.to.generation);
                                        }
                                        None
                                    })
                                    .unwrap()
                            })
                            .min()
                            .unwrap() as f64;

                        let t1 = source.generation as f64;
                        let t2 = target.generation as f64;

                        let div = self.0.get((i, j)).unwrap().to_owned();

                        assert_eq!(distance as f64, t1 - t0 + t2 - t0);
                        pedigree
                            .0
                            .push(Axis(0), array![t0, t1, t2, div].view())
                            .expect("Could not insert row into pedigree");
                    }
                }
            }
        }

        pedigree
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
