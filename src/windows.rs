use std::{
    fmt::Display,
    fs::{self, File, OpenOptions},
    io::{self, BufRead, Write},
};

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use itertools::Itertools;

use crate::{
    arguments::Windows as Args,
    files::lines_from_file,
    genes::{Gene, GenesByStrand, Region},
    methylation_site::MethylationSite,
    *,
};

pub type Window = Vec<MethylationSite>;
#[derive(Debug, PartialEq)]
pub struct Windows {
    pub upstream: Vec<Window>,
    pub gene: Vec<Window>,
    pub downstream: Vec<Window>,
}

impl Windows {
    pub fn new(max_gene_length: u32, args: Args) -> Self {
        let gene_window_count = if args.absolute {
            max_gene_length / args.window_step
        } else {
            100
        };
        let up_down_window_count = if args.absolute {
            args.cutoff / args.window_step
        } else {
            100
        };
        Windows {
            upstream: vec![Vec::new(); up_down_window_count as usize],
            gene: vec![Vec::new(); gene_window_count as usize],
            downstream: vec![Vec::new(); up_down_window_count as usize],
        }
    }

    fn empty() -> Self {
        Windows {
            upstream: Vec::new(),
            downstream: Vec::new(),
            gene: Vec::new(),
        }
    }

    pub fn get(&self, region: Region) -> &Vec<Window> {
        match region {
            Region::Upstream => &self.upstream,
            Region::Gene => &self.gene,
            Region::Downstream => &self.downstream,
        }
    }
    pub fn get_mut<'a>(&'a mut self, location: &Region) -> &'a mut Vec<Window> {
        match location {
            Region::Upstream => &mut self.upstream,
            Region::Gene => &mut self.gene,
            Region::Downstream => &mut self.downstream,
        }
    }

    pub fn iter_upstream(
        &self,
    ) -> std::slice::Iter<'_, std::vec::Vec<methylation_site::MethylationSite>> {
        self.upstream.iter()
    }

    pub fn iter_gene(
        &self,
    ) -> std::slice::Iter<'_, std::vec::Vec<methylation_site::MethylationSite>> {
        self.gene.iter()
    }

    pub fn iter_downstream(
        &self,
    ) -> std::slice::Iter<'_, std::vec::Vec<methylation_site::MethylationSite>> {
        self.downstream.iter()
    }

    pub fn inverse(mut self) -> Self {
        self.upstream = self.downstream.iter().rev().map(|a| a.to_owned()).collect();
        self.gene = self.gene.iter().rev().map(|a| a.to_owned()).collect();
        self.downstream = self.upstream.iter().rev().map(|a| a.to_owned()).collect();
        self
    }

    pub fn steady_state_methylation(&self) -> Vec<f64> {
        let obs = |acc: f64, cur: &MethylationSite| {
            // if cur.count_total == 0 {
            //     acc
            // } else {
            //     acc + cur.count_methylated as f64 / cur.count_total as f64
            // }
            // match cur.status {
            //     MethylationStatus::U => acc,
            //     MethylationStatus::I => acc + 0.5,
            //     MethylationStatus::M => acc + 1.0,
            // }
            acc + cur.meth_lvl
        };

        let mut upstream: Vec<f64> = self
            .upstream
            .iter()
            .map(|w| w.iter().fold(0.0, obs) / w.len() as f64)
            .collect();
        let mut gene: Vec<f64> = self
            .gene
            .iter()
            .map(|w| w.iter().fold(0.0, obs) / w.len() as f64)
            .collect();
        let mut downstream: Vec<f64> = self
            .downstream
            .iter()
            .map(|w| w.iter().fold(0.0, obs) / w.len() as f64)
            .collect();

        gene.append(&mut downstream);
        upstream.append(&mut gene);
        upstream
    }

    pub fn print_steady_state_methylation(methylations: &[f64]) -> String {
        let mut output = String::new();
        for (_i, average) in methylations.iter().enumerate() {
            output.push_str(&format!("{average}\n"));
        }
        output
    }

    /// Print steady state methylations for all windows for all nodes.
    ///
    /// Useful when a single directory containes nodes of several samples.
    pub fn print_all_steady_state_methylations(
        nodes: Vec<String>,
        methylations: Vec<Vec<f64>>,
    ) -> String {
        let mut output = String::new();
        for (_i, (values, node)) in methylations.iter().zip(nodes.iter()).enumerate() {
            output.push_str(node);
            output.push(';');

            for value in values {
                output.push_str(&format!("{value};"));
            }
            output.push('\n');
        }
        output
    }

    pub fn distribution(&self) -> Vec<i32> {
        let mut u: Vec<i32> = self.upstream.iter().map(|w| w.len() as i32).collect();
        let mut g: Vec<i32> = self.gene.iter().map(|w| w.len() as i32).collect();
        let mut d: Vec<i32> = self.downstream.iter().map(|w| w.len() as i32).collect();
        g.append(&mut d);
        u.append(&mut g);
        u
    }

    pub fn print_distribution(distribution: &[i32]) -> String {
        // In CSV format
        let mut output = String::new();

        for (_i, count) in distribution.iter().enumerate() {
            output.push_str(&format!("{count}\n"));
        }

        output
    }

    /// Load an existing extraction into memory for further analysis.
    /// Depends on the standard structure as outputted from this program.
    ///
    /// Important: Methylome files must be in the
    fn load_existing(args: Args, file_name: &str) -> Result<Self> {
        let output_dir = fs::read_dir(&args.output_dir)?;

        let mut result = Windows::empty();
        for dir in output_dir {
            let dir = dir.unwrap();

            if !dir.file_type().unwrap().is_dir() {
                break;
            }

            let region = match dir.path() {
                p if p.ends_with("gene") => Region::Gene,
                p if p.ends_with("downstream") => Region::Downstream,
                p if p.ends_with("upstream") => Region::Upstream,
                _ => break,
            };

            let windows = fs::read_dir(dir.path())?;

            for window in windows {
                let window = window.unwrap();

                if !window.file_type().unwrap().is_dir() {
                    break;
                }

                for file in fs::read_dir(window.path())? {
                    let file = file.unwrap();

                    let name = files::file_name(&file.path());

                    if name != file_name {
                        continue;
                    }

                    let lines = lines_from_file(&file.path())?;

                    let mut sites = Vec::new();
                    let mut error_count = 0;
                    for line in lines {
                        let line = line.unwrap();
                        let m = MethylationSite::from_methylome_file_line(&line, args.invert);

                        match m {
                            Some(site) => sites.push(site),
                            None => error_count += 1,
                        }
                    }

                    if error_count > 5 || sites.len() < 5 {
                        dbg!(error_count, &sites.len());
                        continue;
                    } else {
                        result.get_mut(&region).push(sites);
                    }
                }
            }
        }

        Ok(result)
    }

    pub fn print_all_distributions(nodes: Vec<String>, distibutions: &[Vec<i32>]) -> String {
        let mut output = String::new();
        for (_i, (values, node)) in distibutions.iter().zip(nodes.iter()).enumerate() {
            output.push_str(node);
            output.push(';');

            for value in values {
                output.push_str(&format!("{value};"));
            }
            output.push('\n');
        }
        output
    }

    pub fn save(&self, args: Args, filename: String) -> Result<()> {
        for windows in vec![
            (&self.upstream, "upstream"),
            (&self.gene, "gene"),
            (&self.downstream, "downstream"),
        ]
        .iter()
        {
            for (window, cg_sites) in windows.0.iter().enumerate() {
                let output_file = args.output_dir.join(format!(
                    "{}/{}/{}",
                    windows.1,
                    window * args.window_step as usize,
                    filename
                ));
                let mut file = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&output_file)?;

                let metadata = file.metadata();
                if metadata.unwrap().len() == 0 {
                    // On first write to file, create header line
                    file.write_all("seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl\tcontext.trinucleotide\n".as_bytes())?;
                }
                file.write_all(cg_sites.iter().map(|e| &e.original).join("\n").as_bytes())?;
            }
        }
        Ok(())
    }

    pub fn extract(
        methylome_file: File,
        genome: Vec<GenesByStrand>,
        max_gene_length: u32,
        args: Args,
        file_name: String,
        bars: &MultiProgress,
    ) -> Result<Self> {
        // Check for results already present
        if !args.force {
            if let Ok(windows) = Self::load_existing(args.clone(), &file_name) {
                println!("Skipping extraction as parsable files are found in output directory");
                return Ok(windows);
            }
        }

        let mut last_gene: Option<&Gene> = None;

        let mut windows = Windows::new(max_gene_length, args.clone());

        // size estimation
        const BYTES_PER_LINE: u64 = 34113682 / 950045; // Taken from a random sample, used to estimate number of lines without actually counting
        let n_lines = methylome_file.metadata().unwrap().len() / BYTES_PER_LINE;

        // Progress bars
        let sty = ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-");

        let pb = bars.add(ProgressBar::new(n_lines));
        pb.set_style(sty);
        pb.set_message(file_name);

        let lines = io::BufReader::new(methylome_file).lines();

        // read file, skip header row
        for line_result in lines.skip(1) {
            pb.inc(1);
            if let Ok(line) = line_result {
                // If cg site could not be extracted from a file line, continue with the next line. Happens on header rows, for example.
                let Some(cg) = MethylationSite::from_methylome_file_line(&line, args.invert) else {continue;};

                if last_gene.is_none() || !cg.is_in_gene(last_gene.unwrap(), &args) {
                    last_gene = cg.find_gene(&genome, &args);
                }
                if let Some(gene) = last_gene {
                    cg.place_in_windows(gene, &mut windows, &args);
                    continue;
                }
            }
        }
        pb.finish();

        Ok(windows)
    }
}

impl Display for Windows {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Upstream: {:?}\n\nGene: {:?}\n\nDownstream: {:?}\n\n",
            self.upstream, self.gene, self.downstream
        )
    }
}

#[cfg(test)]
mod test {

    use crate::arguments::Windows as Args;
    #[test]
    fn new_absolute() {
        let args = Args {
            methylome: "/home/constantin/methylome/within_gbM_genes".into(),
            genome: "/home/constantin/methylome/gbM_gene_anotation_extract_Arabidopsis.bed".into(),
            window_size: 512,
            window_step: 256,

            output_dir: "/home/constantin/windows".into(),
            absolute: true,
            cutoff: 2048,
            ..Default::default()
        };
        let windows = super::Windows::new(4096, args);
        assert_eq!(windows.upstream.len(), 8);
        assert_eq!(windows.gene.len(), 16);
        assert_eq!(windows.downstream.len(), 8);
    }
    #[test]
    fn new_relative() {
        let args = Args {
            methylome: "/home/constantin/methylome/within_gbM_genes".into(),
            genome: "/home/constantin/methylome/gbM_gene_anotation_extract_Arabidopsis.bed".into(),
            window_size: 5,
            window_step: 1,

            output_dir: "/home/constantin/windows".into(),
            absolute: false,
            cutoff: 2048,
            ..Default::default()
        };
        let windows = super::Windows::new(4096, args);
        assert_eq!(windows.upstream.len(), 100);
        assert_eq!(windows.gene.len(), 100);
        assert_eq!(windows.downstream.len(), 100);
    }
}
