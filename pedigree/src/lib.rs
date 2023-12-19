mod dataframe;
mod dmatrix;
mod error;

use petgraph::{
    graph::DiGraph, Direction::Incoming,
};
use polars::prelude::*;
use std::{
    fs::{self, File},
    io::{BufRead, Write},
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
};

use dmatrix::DMatrix;
pub use error::Error;
pub type Return<T> = std::result::Result<T, Error>;



use ndarray::{Array2, ArrayView};
#[derive(Debug)]
pub struct Node {
    id: u32,
    name: String,
    generation: u32,
    // meth: bool,
    sites: Option<DataFrame>,
}

impl Node {
    // fn proportion_unmethylated(&self) -> Option<f64> {
    //     match &self.sites {
    //         None => None,
    //         Some(df) => {
    //             let len = df.height() as f64;
    //             let umeth_count = df
    //                 .clone()
    //                 .lazy()
    //                 .group_by("status")
    //                 .agg([col("*").count()])
    //                 // .filter(col("status").eq(lit("U")))
    //                 // .sum()
    //                 .collect()
    //                 .unwrap();

    //             Some(umeth_count["status"].max::<f64>().unwrap() / len)
    //         }
    //     }
    // }

    fn avg_meth_lvl(&self) -> Option<f64> {
        self.sites.as_ref().map(|df| {
            df.clone()
                .lazy()
                .select(&[col("rc.meth.lvl")])
                .mean()
                .collect()
                .unwrap()
                .get(0)
                .unwrap()
                .get(0)
                .unwrap()
                .try_extract::<f64>()
                .unwrap()
        })
    }
}

/// A pedigree describes the divergence bewteen any two samples in a population.
/// It's a matrix with four columns, containing the following information:
///
/// t0: The generation of the last common ancestor between two samples (Can be one of the samples if direct heritage).
///
/// t1: The generation of the first sample.
///
/// t2: The generation of the second sample.
///
/// d: The divergence between the two samples.
///
/// The length of the pedigree is the number of possible pairs of samples, for which methlyation data is available => n * (n - 1) / 2
pub type DivergenceBetweenSamples = Array2<f64>;
type UnmethylatedLevel = f64;

/// Directed Graph containing the pedigree
///
/// Edge weights are the generational time between two samples, "1" by default
#[derive(Debug)]
pub struct Pedigree(DiGraph<Node, u8>);

impl Pedigree {
    /// Read a pedigree from a file.
    ///
    /// The file must be a tab-separated file with four columns:
    ///
    /// `t0`: The generation of the last common ancestor between two samples (Can be one of the samples if direct heritage).
    ///
    /// `t1`: The generation of the first sample.
    ///
    /// `t2`: The generation of the second sample.
    ///
    ///` d`: The divergence between the two samples.
    ///
    /// The first line of the file is ignored.
    /// I chose not to return a result, as this function is meant to statically read a file and therefore it is preferable to panic if the file is not found or parsing errors occur.
    pub fn divergence_from_file(filename: &str) -> DivergenceBetweenSamples {
        let mut divergence = Array2::<f64>::zeros((0, 4));
        let file = std::fs::read_to_string(filename).unwrap();
        file.split('\n').skip(1).for_each(|line| {
            if line.is_empty() {
                return;
            }
            let mut entries = line.split(' ');
            let row = [
                entries.next().unwrap().parse::<f64>().unwrap(),
                entries.next().unwrap().parse::<f64>().unwrap(),
                entries.next().unwrap().parse::<f64>().unwrap(),
                entries.next().unwrap().parse::<f64>().unwrap(),
            ];
            divergence.push_row(ArrayView::from(&row)).unwrap();
        });
        divergence
    }

    pub fn divergence_to_file<P>(&self, path: P, posterior_max_filter: f64) -> std::io::Result<()>
    where
        P: AsRef<Path> + std::fmt::Debug,
    {
        println!("Writing divergence of pedigree to file: {:?}", path);
        let mut file = File::create(path)?;
        let mut content = String::new();
        content += "time0\ttime1\ttime2\tD.value\n";
        for row in self.divergence(posterior_max_filter).rows() {
            content.push_str(&format!("{}\t{}\t{}\t{}\n", row[0], row[1], row[2], row[3]));
        }
        file.write_all(content.as_bytes())
    }

    pub fn divergence(&self, posterior_max_filter: f64) -> DivergenceBetweenSamples {
        let divergence = DMatrix::from(self, posterior_max_filter);
        divergence.unwrap().convert(self)
    }

    pub fn avg_umeth_lvl(&self) -> UnmethylatedLevel {
        self.node_weights()
            .filter_map(|n| n.avg_meth_lvl())
            .map(|n| 1.0 - n)
            .sum::<f64>()
            / self.node_count() as f64
    }

    /// Iteratator over all nodes of the pedigree that contain a sites dataframe
    pub fn dfs(&self) -> impl Iterator<Item = &DataFrame> {
        self.node_weights()
            .filter_map(|n| n.sites.as_ref())
            .collect::<Vec<_>>()
            .into_iter()
    }

    pub fn nodes_with_sites(&self) -> impl Iterator<Item = &Node> {
        self.node_weights()
            .filter(|n| n.sites.is_some())
            .collect::<Vec<_>>()
            .into_iter()
    }

    /// Number of nodes in the pedigree that contain a sites dataframe
    pub fn df_count(&self) -> usize {
        self.node_weights()
            .fold(0, |acc, n| if n.sites.is_some() { acc + 1 } else { acc })
    }

    pub fn build<P>(
        nodelist: P,
        edgelist: P,
        _posterior_max_filter: f64,
        custom_column_names: Option<Vec<&'static str>>,
    ) -> Return<Self>
    where
        P: AsRef<Path>,
    {
        let mut p = DiGraph::new();
        let _nodes = dataframe::load(&nodelist, Some(&["filename", "node", "gen", "meth"]));
        let _edges = dataframe::load(&edgelist, Some(&["from", "to", "dt"]));

        let nodes = fs::read_to_string(&nodelist)?;
        let edges = fs::read_to_string(&edgelist)?;
        let pb = progress_bars::Progress::new("Loading Methylation Data", nodes.len());
        for (i, line) in nodes.split(['\n', '\r']).skip(1).enumerate() {
            pb.inc(1);
            if line.is_empty() {
                continue;
            }
            let mut entries = line.split([',', '\t', ' ']);
            let file = entries.next().ok_or(Error::NoFile(line.into()))?;
            let name = entries.next().ok_or(Error::NoName(line.into()))?;
            let generation = entries
                .next()
                .ok_or(Error::NoGeneration(line.into()))?
                .parse::<u32>()
                .map_err(|_| Error::GenerationUnparseable(line.into()))?;
            let meth = entries.next().ok_or(Error::NoMethylation(line.into()))?;
            let meth = meth == "Y";

            let node = Node {
                id: i as u32,
                name: name.into(),
                generation,

                sites: if meth {
                    let file = if file.starts_with('/') {
                        file.into()
                    } else {
                        let mut full_path = PathBuf::from(nodelist.as_ref().parent().unwrap());
                        full_path.push(file);
                        full_path
                    };

                    if !file.exists() {
                        return Err(Error::MethylationFile(file));
                    }

                    Some(dataframe::load(file, custom_column_names.as_deref()))
                } else {
                    None
                },
            };

            p.add_node(node);
        }
        pb.finish();

        for line in edges.split(['\n', '\r']).skip(1) {
            if line.is_empty() {
                continue;
            }

            let mut entries = line.split(['\t', ' ', ',']);
            let from = entries.next().ok_or(Error::NoFrom)?;
            let to = entries.next().ok_or(Error::NoTo)?;
            let dt = entries
                .next()
                .map_or(1, |t| t.parse().expect("Could not parse time difference"));

            let find_nx_by_name = |name: &str| {
                p.node_indices()
                    .find(|&nx| p.node_weight(nx).unwrap().name == name)
            };

            let from_nx = find_nx_by_name(from)
                .ok_or(Error::EdgelistContainingNonExistentNode(from.into()))?;
            let to_nx =
                find_nx_by_name(to).ok_or(Error::EdgelistContainingNonExistentNode(to.into()))?;

            p.add_edge(from_nx, to_nx, dt);
        }

        Ok(Pedigree(p))

        // let mut nodes: Vec<Node> = nodes.iter().filter(|n| n.meth).cloned().collect();
    }

    fn sites_per_node(&self) -> usize {
        self.node_weights()
            .find_map(|n| n.sites.as_ref())
            .map(|df| df.height())
            .unwrap_or(0)
    }

    fn get_predecessor_sites(&self, node: &Node) -> Option<&DataFrame> {
        let pred = self.neighbors_directed(node.id.into(), Incoming).next()?;

        match &self.node_weight(pred)?.sites {
            None => self.get_predecessor_sites(self.node_weight(pred)?),
            Some(df) => Some(df),
        }
    }

    /// Based on the inheritance pattern defined by the pedigree graph
    /// this method filters out all the sites that don't meet a stability requirement.
    ///
    /// The stability is quantified by the total number of switches between methylated and unmethylated
    /// states in the lineage of a site.
    ///
    /// Inclusive min, exclusive max
    pub fn filter_by_inheritance_stability(mut self, min: u32, max: u32) -> Return<Self> {
        let mut scores = Series::new("scores", vec![0; self.sites_per_node()]);
        for (_i, node) in self.node_weights().enumerate() {
            if let Some(df) = &node.sites {
                let status = df.column("status")?;

                let pred_sites = self.get_predecessor_sites(node);

                match pred_sites {
                    Some(pred_sites) => {
                        let pred_status = pred_sites.column("status")?;
                        let diff = status - pred_status;
                        scores = scores + diff;
                    }
                    None => continue,
                }
            }
        }

        let scores: BooleanChunked = scores
            .rechunk()
            .iter()
            .map(|s| s.ge(&AnyValue::UInt32(min)) && s.lt(&AnyValue::UInt32(max)))
            .collect_ca("scores");

        self.node_weights_mut()
            .for_each(|n| match n.sites.as_ref() {
                None => {}
                Some(df) => {
                    let df = df.clone().filter(&scores).unwrap();
                    n.sites = Some(df);
                }
            });

        Ok(self)
    }
}

impl Deref for Pedigree {
    type Target = DiGraph<Node, u8>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Pedigree {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    // #[test]
    // fn proportion_unmethylated() {
    //     let path = PathBuf::from("../data/methylome/G0.txt");

    //     let df = load(path, None);

    //     let node = Node {
    //         id: 0,
    //         name: "G0".into(),
    //         generation: 0,
    //         sites: Some(df),
    //     };
    //     assert!(node.proportion_unmethylated().is_some());
    //     assert_eq!(node.proportion_unmethylated(), Some(0.5));
    // }

    #[test]
    fn build_pedigree() {
        let nodelist = Path::new("../data/nodelist.txt");
        let edgelist = Path::new("../data/edgelist.txt");

        let div = Pedigree::build(nodelist, edgelist, 0.99, None)
            .unwrap()
            .divergence(0.99);

        assert_eq!(div.shape(), &[4 * 3 / 2, 4]);
    }

    // #[test]
    // fn wildtype_pedigree() {
    //     let nodelist = Path::new("./data/nodelist.txt");
    //     let edgelist = Path::new("./data/edgelist.txt");
    //     let comparison = Pedigree::from_file(
    //         "./data/desired_output/pedigree-pdata_epimutation_rate_estimation_window_gene_0.txt",
    //     );

    //     let pedigree = Pedigree::build(nodelist, edgelist, 0.99).expect("Could not build pedigree");

    //     // TODO: Enable
    //     // assert_close!(pedigree.1, 0.991008120326199);

    //     let pedigree = pedigree.0;

    //     // for (i, j) in pedigree.iter().zip(comparison.iter()) {
    //     // TODO: Enable
    //     // assert_close!(i, j);
    //     // }
    // }
}
