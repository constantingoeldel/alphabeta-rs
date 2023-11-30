mod dmatrix;

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
#[derive(Clone, Debug)]
struct Node {
    id: usize,
    file: PathBuf,
    name: String,
    generation: u32,
    meth: bool,
    proportion_unmethylated: Option<f64>,
    rc_meth_lvl: Option<f64>,
    sites: Option<Vec<MethylationSite>>, // Might lead to memory issues if there are too many sites. But it makes lookup really fast.
}
#[derive(Debug)]
struct Edge<'a> {
    from: &'a Node,
    to: &'a Node,
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
#[derive(Clone, Debug)]
pub struct Pedigree(Array2<f64>);
type UnmethylatedLevel = f64;

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
    pub fn from_file(filename: &str) -> Self {
        let mut pedigree = Array2::<f64>::zeros((0, 4));
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
            pedigree.push_row(ArrayView::from(&row)).unwrap();
        });
        Pedigree(pedigree)
    }

    pub fn to_file<P>(&self, path: P) -> std::io::Result<()>
    where
        P: AsRef<Path> + std::fmt::Debug,
    {
        println!("Writing pedigree to file: {:?}", path);
        let mut file = File::create(path)?;
        let mut content = String::new();
        content += "time0\ttime1\ttime2\tD.value\n";
        for row in self.rows() {
            content.push_str(&format!("{}\t{}\t{}\t{}\n", row[0], row[1], row[2], row[3]));
        }
        file.write_all(content.as_bytes())
    }
    pub fn build<P>(
        nodelist: P,
        edgelist: P,
        posterior_max_filter: f64,
    ) -> Return<(Self, UnmethylatedLevel)>
    where
        P: AsRef<Path>,
    {
        let nodes = fs::read_to_string(nodelist)?;
        let edges = fs::read_to_string(edgelist)?;
        let nodes: Vec<Node> = nodes
            .split(['\n', '\r'])
            .skip(1)
            .enumerate()
            .filter_map(|(i, line)| {
                let mut entries = line.split([',', '\t', ' ']);
                Some(Node {
                    id: i,
                    file: PathBuf::from(entries.next()?),
                    name: String::from(entries.next()?),
                    generation: entries.next()?.parse::<u32>().ok()?,
                    meth: entries.next()? == "Y",
                    proportion_unmethylated: None,
                    rc_meth_lvl: None,
                    sites: None,
                })
            })
            .collect();
        if nodes.is_empty() {
            return Err(Error::NodeParsing);
        }

        let edges: Vec<Edge> = edges
            .split(['\n', '\r'])
            .skip(1)
            .filter_map(|line| {
                let mut entries = line.split(['\t', ' ', ',']);
                let from = entries.next()?;
                let to = entries.next()?;
                Some(Edge {
                    from: nodes.iter().find(|n| n.name == from)?,
                    to: nodes.iter().find(|n| n.name == to)?,
                })
            })
            .collect();

        let mut nodes: Vec<Node> = nodes.iter().filter(|n| n.meth).cloned().collect();

        let pb = progress_bars::Progress::new("Loading Methylation Data", nodes.len());

        for node in nodes.iter_mut() {
            pb.inc(1);
            // let mut file = PathBuf::from(nodelist.parent().unwrap());
            // file.push(&node.file);
            let f = File::open(&node.file).map_err(Error::NodeFile)?;
            let reader = BufReader::new(f);

            let mut sites = Vec::new();

            for line in reader.lines() {
                let line = line.expect("Could not read line");
                let methylation = MethylationSite::from_methylome_file_line(&line, false);

                if methylation.is_none() {
                    continue;
                }

                let methylation = methylation.unwrap();
                sites.push(methylation);
            }

            let valid_sites: Vec<&MethylationSite> = sites
                .iter()
                .filter(|s| s.posteriormax >= posterior_max_filter)
                .collect();

            let p_uu_count = valid_sites
                .iter()
                .filter(|s| s.status == MethylationStatus::U)
                .count();

            let p_uu_share = p_uu_count as f64 / valid_sites.len() as f64;

            let avg_meth_lvl =
                valid_sites.iter().map(|s| s.meth_lvl).sum::<f64>() / valid_sites.len() as f64;

            node.proportion_unmethylated = Some(p_uu_share);
            node.rc_meth_lvl = Some(avg_meth_lvl);
            node.sites = Some(sites);
        }
        pb.finish();

        let tmp0uu_meth_lvl = nodes
            .iter()
            .map(|n| 1.0 - n.rc_meth_lvl.unwrap())
            .sum::<f64>()
            / nodes.len() as f64;

        let divergence = DMatrix::from(&nodes, posterior_max_filter);
        let pedigree = divergence.convert(&nodes, &edges);
        Ok((pedigree, tmp0uu_meth_lvl))
    }
}

impl Deref for Pedigree {
    type Target = Array2<f64>;

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

        let pedigree = Pedigree::build(nodelist, edgelist, 0.99).expect("Could not build pedigree");

        assert_eq!(pedigree.0.shape(), &[4 * 3 / 2, 4]);
        pedigree
            .0
            .to_file(Path::new("./data/pedigree_generated.txt"))
            .unwrap();
        // TODO: enable
        // assert_close!(pedigree.1, 0.4567024);
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
