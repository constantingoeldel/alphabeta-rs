use std::{
    fmt::Display,
    fs::{self, File, OpenOptions},
    io::{self, BufRead, Write},
};

use itertools::Itertools;
use methylome::{MethylationSite, Strand};
use progress_bars::{MultiProgress, ProgressBar};

use crate::{
    config::{get, Metaprofile as Args},
    files::lines_from_file,
    genes::{find_gene, is_in_gene, Gene, Genome, Region},
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
            100 / args.window_step
        };
        let up_down_window_count = if args.absolute {
            args.cutoff / args.window_step
        } else {
            100 / args.window_step
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

    pub fn get_mut<'a>(&'a mut self, location: &Region) -> &'a mut Vec<Window> {
        match location {
            Region::Upstream => &mut self.upstream,
            Region::Gene => &mut self.gene,
            Region::Downstream => &mut self.downstream,
        }
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
    fn load_existing(args: Args, file_name: &str) -> Return<Self> {
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

    pub fn save(&self, args: Args, filename: String) -> Return<()> {
        for windows in [
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
                    .write(true)
                    .create(true)
                    .open(&output_file)
                    .unwrap();

                file.write_all("seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl\tcontext.trinucleotide\n".as_bytes())?;
                file.write_all(cg_sites.iter().map(|e| &e.original).join("\n").as_bytes())?;
            }
        }
        Ok(())
    }

    pub fn extract(
        methylome_file: File,
        genome: Genome,
        max_gene_length: u32,
        args: Args,
        file_name: String,
        bars: &MultiProgress,
    ) -> Return<Self> {
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

        let pb = bars.add(ProgressBar::new(n_lines));
        pb.set_style(progress_bars::progress_style());
        pb.set_message(file_name);

        let lines = io::BufReader::new(methylome_file).lines();

        // read file, skip header row
        for line_result in lines.skip(1) {
            pb.inc(1);
            if let Ok(line) = line_result {
                // If cg site could not be extracted from a file line, continue with the next line. Happens on header rows, for example.
                let Some(cg) = MethylationSite::from_methylome_file_line(&line, args.invert) else {
                    continue;
                };

                if last_gene.is_none()
                    || !is_in_gene(
                        &cg,
                        last_gene.unwrap(),
                        // TODO: MOVE INTO FUNCTION
                        get().cutoff_gene_length,
                        get().cutoff,
                    )
                {
                    last_gene = find_gene(&cg, &genome, get().cutoff_gene_length, get().cutoff);
                }
                if let Some(gene) = last_gene {
                    place_in_windows(&cg, gene, &mut windows, &args);
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

    use crate::config::Metaprofile as Args;
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

/// Place a CG site in the correct windows. Returns a list of all the successfull insertions as a tuple of the region (upstream, downstream or gene) and the index of the window.
///
/// It works by first finding the region the CG site is in (upstream, downstream or gene) and then finding the windows within that a CG site belongs to.
/// For genes on the - strand, the windows are reversed, so that the first window is the one closest to the end of the gene.
pub fn place_in_windows(
    site: &MethylationSite,
    gene: &Gene,
    windows: &mut Windows,
    args: &Args,
) -> Vec<(Region, usize)> // Return a vector of (strand, window) tuples for each window the CG site is in
{
    const E: f64 = 0.1; // Epsilon for floating point comparison
    let location = site.start as f64;
    let cutoff = args.cutoff as f64;
    let step = args.window_step as f64;
    let size = args.window_size as f64;
    let start = gene.start as f64;
    let end = gene.end as f64;
    let length = end - start;

    if let Strand::Unknown = site.strand {
        let mut copy = site.clone();
        copy.strand = Strand::Sense;
        return place_in_windows(&copy, gene, windows, args);
    }

    // Offset from start for + strand, offset from end for - strand. Can be negative for upstream sites
    let mut windows_in = Vec::new();
    let offset = match &site.strand {
        Strand::Sense => location - start,
        Strand::Antisense => end - location,
        _ => return windows_in, // TODO better handling
    };
    let region = match offset {
        x if x < 0.0 => Region::Upstream,
        x if x > length => Region::Downstream, // CG site exactly on the end of the gene is still considered in the gene
        _ => Region::Gene,
    };
    let local_windows = windows.get_mut(&region);

    // let max = if args.absolute { gene_length } else { 100 };
    let mut position = match (&region, &site.strand) {
        // Position within the region of the gene, switched start & end for - strand
        (Region::Upstream, Strand::Sense) => location - start + cutoff,
        (Region::Gene, Strand::Sense) => location - start,
        (Region::Downstream, Strand::Sense) => location - end,
        (Region::Upstream, Strand::Antisense) => end - location + cutoff,
        (Region::Gene, Strand::Antisense) => end - location,
        (Region::Downstream, Strand::Antisense) => start - location,
        _ => return windows_in,
    };

    if !args.absolute {
        position = match region {
            Region::Upstream => position / cutoff,
            Region::Gene => position / length,
            Region::Downstream => position / cutoff,
        };
        position *= 100.0; // Normalize to 0-100%
    }

    for (i, window) in local_windows.iter_mut().enumerate() {
        let lower_bound = i as f64 * step - E;
        let upper_bound = lower_bound + size + E;

        if position >= lower_bound && position <= upper_bound {
            window.push(site.clone());
            windows_in.push((region.clone(), i));
        }
    }
    windows_in
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{config, genes::Gene};
    use methylome::{Chromosome, Strand};

    #[test]
    fn test_place_site_absolute() {
        let args = config::Metaprofile {
            invert: false,
            absolute: true,
            cutoff: 1000,
            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let all_within_gene = Gene {
            annotation: String::new(),

            chromosome: Chromosome::Numbered(1),
            start: 1000,
            end: 2000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 2000,
            end: 3000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 0,
            end: 1000,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(1000, args.clone());
        for i in 1..1000 {
            let cg = MethylationSite::new(i + 1000, Strand::Antisense);
            let upstream = place_in_windows(&cg, &all_upstream_gene, &mut windows, &args);
            let gene = place_in_windows(&cg, &all_within_gene, &mut windows, &args);
            let downstream = place_in_windows(&cg, &all_downstream_gene, &mut windows, &args);

            println!("Placing {i}");
            println!("Upstream: {upstream:?}");
            println!("Gene: {gene:?}");
            println!("Downstream: {downstream:?}");
            // TODO: enable
            // assert!(windows.upstream[i as usize].contains(&cg));
            // assert!(windows.gene[i as usize].contains(&cg));
            // assert!(windows.downstream[i as usize].contains(&cg));

            // if i > 3 {
            //     assert!(windows.upstream[i as usize - 1].contains(&cg));
            //     assert!(windows.gene[i as usize - 1].contains(&cg));
            //     assert!(windows.downstream[i as usize - 1].contains(&cg));
            //     assert!(windows.upstream[i as usize - 2].contains(&cg));
            //     assert!(windows.gene[i as usize - 2].contains(&cg));
            //     assert!(windows.downstream[i as usize - 2].contains(&cg));
            // }
        }
    }
    #[test]
    fn test_place_site_relative_acting_like_absolute() {
        let args = Args {
            invert: false,
            absolute: false,
            cutoff: 100,

            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let all_within_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 100,
            end: 200,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 200,
            end: 300,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 0,
            end: 100,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(100, args.clone());
        for i in 1..100 {
            let cg = MethylationSite::new(i + 100, Strand::Sense);
            let upstream = place_in_windows(&cg, &all_upstream_gene, &mut windows, &args);
            let gene = place_in_windows(&cg, &all_within_gene, &mut windows, &args);
            let downstream = place_in_windows(&cg, &all_downstream_gene, &mut windows, &args);

            println!("Placing {i}");
            println!("Upstream: {upstream:?}");
            println!("Gene: {gene:?}");
            println!("Downstream: {downstream:?}");
            assert!(windows.upstream[i as usize].contains(&cg));
            assert!(windows.gene[i as usize].contains(&cg));
            assert!(windows.downstream[i as usize].contains(&cg));
        }
    }
    #[test]
    fn test_place_site_relative() {
        let args = Args {
            invert: false,
            absolute: false,
            cutoff: 1000,

            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let all_within_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 1000,
            end: 2000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 2000,
            end: 3000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 0,
            end: 1000,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(1000, args.clone());
        assert!(windows.upstream.len() == 100);
        for i in 1..1000 {
            let cg = MethylationSite::new(i + 1000, Strand::Sense);
            let upstream = place_in_windows(&cg, &all_upstream_gene, &mut windows, &args);
            let gene = place_in_windows(&cg, &all_within_gene, &mut windows, &args);
            let downstream = place_in_windows(&cg, &all_downstream_gene, &mut windows, &args);

            println!("Placing {i}");
            println!("Upstream: {upstream:?}");
            println!("Gene: {gene:?}");
            println!("Downstream: {downstream:?}");
            println!("{}: {}", i / 10, windows.upstream[(i / 10) as usize].len());
            assert!(windows.upstream[(i / 10) as usize].contains(&cg));
            assert!(windows.gene[(i / 10) as usize].contains(&cg));
            assert!(windows.downstream[(i / 10) as usize].contains(&cg));
        }
    }

    #[test]
    fn test_place_site() {
        let cg_a = MethylationSite::new(80, Strand::Sense);
        let cg_b = MethylationSite::new(100, Strand::Sense);
        let cg_c = MethylationSite::new(123, Strand::Sense);
        let cg_d = MethylationSite::new(200, Strand::Sense);
        let cg_e = MethylationSite::new(201, Strand::Sense);
        let cg_f = MethylationSite::new(512 + 100 + 100, Strand::Sense);
        let cg_g = MethylationSite::new(1024 + 100 + 100, Strand::Sense);
        let cg_h = MethylationSite::new(2048 + 100 + 100, Strand::Sense);

        let gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 100,
            end: 200,
            strand: Strand::Sense,
            name: String::new(),
        };

        let args = Args {
            invert: false,
            absolute: false,
            cutoff: 2048,

            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let mut windows = Windows::new(1000, args.clone());

        place_in_windows(&cg_a, &gene, &mut windows, &args);
        place_in_windows(&cg_b, &gene, &mut windows, &args);
        place_in_windows(&cg_c, &gene, &mut windows, &args);
        place_in_windows(&cg_d, &gene, &mut windows, &args);
        place_in_windows(&cg_e, &gene, &mut windows, &args);
        place_in_windows(&cg_f, &gene, &mut windows, &args);
        place_in_windows(&cg_g, &gene, &mut windows, &args);
        place_in_windows(&cg_h, &gene, &mut windows, &args);

        println!("{windows}");
        assert!(windows.upstream[98].contains(&cg_a));
        assert!(windows.upstream[99].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[21].contains(&cg_c));
        assert!(windows.gene[22].contains(&cg_c));
        assert!(windows.gene[23].contains(&cg_c));
        assert!(windows.gene[99].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[24].contains(&cg_f));
        assert!(windows.downstream[49].contains(&cg_g));
        assert!(windows.downstream[99].contains(&cg_h));
    }
    #[test]
    fn test_place_site_absolute_2() {
        let cg_a = MethylationSite::new(80, Strand::Sense);
        let cg_b = MethylationSite::new(100, Strand::Sense);
        let cg_c = MethylationSite::new(123, Strand::Sense);
        let cg_d = MethylationSite::new(200, Strand::Sense);
        let cg_e = MethylationSite::new(201, Strand::Sense);
        let cg_f = MethylationSite::new(220, Strand::Sense);

        let gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 100,
            end: 200,
            strand: Strand::Sense,
            name: String::new(),
        };

        let args = Args {
            absolute: true,
            cutoff: 2048,

            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let mut windows = Windows::new(100, args.clone());

        place_in_windows(&cg_a, &gene, &mut windows, &args);
        place_in_windows(&cg_b, &gene, &mut windows, &args);
        place_in_windows(&cg_c, &gene, &mut windows, &args);
        place_in_windows(&cg_d, &gene, &mut windows, &args);
        place_in_windows(&cg_e, &gene, &mut windows, &args);
        place_in_windows(&cg_f, &gene, &mut windows, &args);
        assert!(windows.upstream[2026].contains(&cg_a));
        assert!(windows.upstream[2027].contains(&cg_a));
        assert!(windows.upstream[2028].contains(&cg_a));
        assert!(windows.gene[0].contains(&cg_b));
        assert!(windows.gene[21].contains(&cg_c));
        assert!(windows.gene[22].contains(&cg_c));
        assert!(windows.gene[23].contains(&cg_c));
        assert!(windows.gene[99].contains(&cg_d));
        assert!(windows.downstream[0].contains(&cg_e));
        assert!(windows.downstream[1].contains(&cg_e));
        assert!(windows.downstream[18].contains(&cg_f));
        assert!(windows.downstream[19].contains(&cg_f));
        assert!(windows.downstream[20].contains(&cg_f));
    }

    #[test]
    fn test_place_site_relative_antisense() {
        let args = Args {
            absolute: false,
            cutoff: 1000,

            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let all_within_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 1000,
            end: 2000,
            strand: Strand::Antisense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 2000,
            end: 3000,
            strand: Strand::Antisense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 0,
            end: 1000,
            strand: Strand::Antisense,
            name: String::new(),
        };

        let mut windows = Windows::new(1000, args.clone());
        assert!(windows.upstream.len() == 100);
        for i in 1..1000 {
            let cg = MethylationSite::new(i + 1000, Strand::Antisense);
            let upstream = place_in_windows(&cg, &all_upstream_gene, &mut windows, &args);
            let gene = place_in_windows(&cg, &all_within_gene, &mut windows, &args);
            let downstream = place_in_windows(&cg, &all_downstream_gene, &mut windows, &args);

            println!("Placing {i}");
            println!("Upstream: {upstream:?}");
            println!("Gene: {gene:?}");
            println!("Downstream: {downstream:?}");
            println!(
                "{}: {}",
                (999 - i) / 10,
                windows.upstream[(i / 10) as usize].len()
            );
            assert!(windows.upstream[((999 - i) / 10) as usize].contains(&cg));
            assert!(windows.gene[((999 - i) / 10) as usize].contains(&cg));
            assert!(windows.downstream[((999 - i) / 10) as usize].contains(&cg));
        }
    }
    #[test]
    fn test_place_site_absolute_invert() {
        let args = Args {
            absolute: true,
            cutoff: 1000,

            window_size: 2,
            window_step: 1,
            ..Default::default()
        };
        let all_within_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 1000,
            end: 2000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_upstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 2000,
            end: 3000,
            strand: Strand::Sense,
            name: String::new(),
        };
        let all_downstream_gene = Gene {
            annotation: String::new(),
            chromosome: Chromosome::Numbered(1),
            start: 0,
            end: 1000,
            strand: Strand::Sense,
            name: String::new(),
        };

        let mut windows = Windows::new(1000, args.clone());
        assert!(windows.upstream.len() == 1000);
        for i in 1..1000 {
            let cg = MethylationSite::new(i + 1000, Strand::Sense);
            let upstream = place_in_windows(&cg, &all_upstream_gene, &mut windows, &args);
            let gene = place_in_windows(&cg, &all_within_gene, &mut windows, &args);
            let downstream = place_in_windows(&cg, &all_downstream_gene, &mut windows, &args);

            println!("Placing {i}");
            println!("Upstream: {upstream:?}");
            println!("Gene: {gene:?}");
            println!("Downstream: {downstream:?}");
            println!("{}: {}", (i), windows.upstream[(i) as usize].len());
            assert!(windows.upstream[(i) as usize].contains(&cg));
            assert!(windows.gene[(i) as usize].contains(&cg));
            assert!(windows.downstream[(i) as usize].contains(&cg));
        }
    }
}
