mod dmatrix;

use petgraph::{data::FromElements, graph::DiGraph, visit::IntoEdges};

use std::{
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
};

mod error;
use dmatrix::DMatrix;
pub use error::Error;
pub type Return<T> = std::result::Result<T, Error>;

use methylome::{MethylationSite, MethylationStatus};

use ndarray::{Array2, ArrayView};
#[derive(Debug)]
struct Node {
    id: usize,
    name: String,
    generation: u32,
    meth: bool,
    sites: Vec<MethylationSite>, // Might lead to memory issues if there are too many sites. But it makes lookup really fast.
}

impl Node {
    fn proportion_unmethylated(&self) -> f64 {
        let p_uu_count = self
            .sites
            .iter()
            .filter(|s| s.status == MethylationStatus::U)
            .count();

        p_uu_count as f64 / self.sites.len() as f64
    }

    fn avg_meth_lvl(&self) -> f64 {
        self.sites.iter().map(|s| s.meth_lvl).sum::<f64>() / self.sites.len() as f64
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
#[derive(Clone, Debug)]
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
        todo!();

        // let divergence = DMatrix::from(self.nodes(), posterior_max_filter);
        // divergence.convert(self.nodes(), self.edges())
    }

    pub fn avg_umeth_lvl(&self) -> UnmethylatedLevel {
        self.node_weights()
            .map(|n| 1.0 - n.rc_meth_lvl.unwrap())
            .sum::<f64>()
            / self.node_count() as f64
    }

    pub fn build<P>(nodelist: P, edgelist: P, posterior_max_filter: f64) -> Return<Self>
    where
        P: AsRef<Path>,
    {
        let mut p = DiGraph::new();
        let nodes = fs::read_to_string(nodelist)?;
        let edges = fs::read_to_string(edgelist)?;
        let pb = progress_bars::Progress::new("Loading Methylation Data", nodes.len());
        for (i, line) in nodes.split(['\n', '\r']).skip(1).enumerate() {
            pb.inc(1);
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

            if meth {
                let file = if file.starts_with("/") {
                    file.into()
                } else {
                    let mut file = PathBuf::from(nodelist.as_ref().parent().unwrap());
                    file.push(file);
                    file
                };

                if !file.exists() {
                    return Err(Error::MethylationFile(file));
                }
            }

            let node = Node {
                id: i,
                name: name.into(),
                generation,
                meth,

                sites: if meth {
                    let f = File::open(file).map_err(Error::NodeFile)?;
                    let reader = BufReader::new(f);
                    reader
                        .lines()
                        .filter_map(|l| {
                            MethylationSite::from_methylome_file_line(
                                &l.expect("Could not read methylation file line"),
                            )
                        })
                        .collect::<Vec<MethylationSite>>()
                } else {
                    Vec::new()
                },
            };

            p.add_node(node);
        }
        pb.finish();

        for line in edges.split(['\n', '\r']).skip(1) {
            let mut entries = line.split(['\t', ' ', ',']);
            let from = entries.next().ok_or(Error::NoFrom)?;
            let to = entries.next().ok_or(Error::NoTo)?;
            let dt = entries
                .next()
                .map_or(1, |t| t.parse().expect("Could not parse time difference"));

            let find_nx_by_name = |name: &str| {
                for nx in p.node_indices() {
                    if p.node_weight(nx).unwrap().name == name {
                        return Some(nx);
                    }
                }
                None
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
    #[test]
    fn build_pedigree() {
        let nodelist = Path::new("./data/nodelist.txt");
        let edgelist = Path::new("./data/edgelist.txt");

        let div = Pedigree::build(nodelist, edgelist, 0.99)
            .expect("Could not build pedigree")
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
