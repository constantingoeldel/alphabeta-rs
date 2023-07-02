use std::fmt::Display;

use crate::{genes::Genome, Error};
use arguments::Windows as Args;
use itertools::Itertools;

#[macro_export]
macro_rules! print_dev {
    ($($rest:tt)*) => {
        #[cfg(debug_assertions)]
        std::println!($($rest)*)
    }
}

use crate::{
    genes::{Gene, Region, Strand},
    windows::Windows,
    *,
};

/// Basic Data Type for a Methylation Site
///
/// Includes conversion from a .bed file line in several formats. If your program fails,
/// check if the format of your .bed file is supported.
///
/// If not, simply add your own conversion function and submit a pull request.
///
/// Most times, only one position is mentioned, but in some cases, a range is given.
/// If only one location, we assign `start` to that and `end` to `start + 1`, as a cg site is two nucleotides long.
///
#[derive(Clone, PartialEq, Debug)]
pub struct MethylationSite {
    pub chromosome: Chromosome,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub original: String,
    pub context: String,
    pub count_methylated: u32,
    pub count_total: u32,
    pub posteriormax: f64,
    pub status: MethylationStatus,
    pub meth_lvl: f64,
    pub context_trinucleotide: String,
}

// TODO: Is there a better name for the regular chromosomes than "Numbered"?
#[derive(Clone, PartialEq, Debug, Eq, Hash)]
pub enum Chromosome {
    Numbered(u8),
    Mitochondrial,
    Chloroplast,
}

impl TryFrom<&str> for Chromosome {
    type Error = Error;
    fn try_from(s: &str) -> Result<Self, Error> {
        let s = s.trim_start_matches("chr");
        match s {
            "M" => Ok(Chromosome::Mitochondrial),
            "C" => Ok(Chromosome::Chloroplast),
            _ => match s.parse::<u8>() {
                Ok(n) => Ok(Chromosome::Numbered(n)),
                Err(_) => Err(Error::Chromosome(s.to_string())),
            },
        }
    }
}

impl Default for MethylationSite {
    fn default() -> Self {
        MethylationSite {
            chromosome: Chromosome::Numbered(1),
            start: 1,
            end: 2,
            strand: Strand::Unknown,
            original: String::new(),
            context: String::new(),
            count_methylated: 1,
            count_total: 1,
            posteriormax: 0.999,
            status: MethylationStatus::M,
            meth_lvl: 0.1,
            context_trinucleotide: String::new(),
        }
    }
}

/// The three different kinds of methylations statusses that are distinguished in the AlphaBeta paper.
///
/// U: unmethylated on both alleles
/// M: methylated on both alleles
/// I: methylated on one allele, unmethylated on the other (intermediate)
#[derive(Clone, PartialEq, Debug)]
pub enum MethylationStatus {
    U,
    M,
    I,
}

impl From<char> for MethylationStatus {
    fn from(c: char) -> Self {
        match c {
            'M' => MethylationStatus::M,
            'I' => MethylationStatus::I,
            'U' => MethylationStatus::U,
            _ => {
                println!(
                    "Warning: Encountered invalid methylation status: {c}. Parsed as Unmethylated"
                );
                MethylationStatus::U
            }
        }
    }
}

impl MethylationSite {
    pub fn new(location: u32, strand: Strand) -> Self {
        MethylationSite {
            start: location,
            end: location + 1,
            strand,
            ..Default::default()
        }
    }
    /// Convert status to numeric value following the AlphaBeta paper.
    /// u/u => 0
    /// u/m => 1
    /// m/m => 2
    pub fn status_numeric(&self) -> u32 {
        match self.status {
            MethylationStatus::M => 2,
            MethylationStatus::U => 0,
            MethylationStatus::I => 1,
        }
    }

    /// Create a new CG site from a line of a methylation file.
    /// Only yields a CG site if the line is formatted correctly and is a CG site.
    /// If invalid, an error is returned.
    ///
    /// One pitfall of this implementation is the `collect tuple` call, which only yields a `Some` value if the line has exactly 9 or 10 tab-separated fields.
    /// Some files provide the trinucleotide, others don't, so there are two versions.
    ///
    /// If you need to parse a different format, simply add a new function and submit a pull request.
    pub fn from_methylome_file_line(s: &str, invert_strand: bool) -> Option<Self> {
        let first_format = |s: &str| {
            s.split('\t')
                .collect_tuple()
                .filter(|(_, _, _, context, _, _, _, _, _)| context == &"CG")
                .map(
                    |(
                        chromosome,
                        location,
                        strand,
                        context,
                        count_methylated,
                        count_total,
                        posteriormax,
                        status,
                        meth_lvl,
                    )| {
                        Ok(MethylationSite {
                            chromosome: chromosome.try_into()?,
                            start: location.parse::<u32>()?,
                            end: location.parse::<u32>()? + 1,
                            strand: if (strand == "+") ^ invert_strand {
                                Strand::Sense
                            } else {
                                Strand::Antisense
                            },
                            original: s.to_owned(),
                            context: String::from(context),
                            count_methylated: count_methylated.parse::<u32>()?,
                            count_total: count_total.parse::<u32>()?,
                            posteriormax: posteriormax.parse::<f64>()?,
                            status: status
                                .chars()
                                .next()
                                .ok_or(Error::Simple("Status could not included in file"))?
                                .into(),
                            meth_lvl: meth_lvl.parse::<f64>()?,
                            context_trinucleotide: String::from("XXX"),
                        })
                    },
                )
                .unwrap_or(Err(Error::MethlyationSiteFormat))
        };

        let second_format = |s: &str| {
            s.split('\t')
                .collect_tuple()
                .filter(|(_, _, _, context, _, _, _, _, _, _)| context == &"CG")
                .map(
                    |(
                        chromosome,
                        location,
                        strand,
                        context,
                        count_methylated,
                        count_total,
                        posteriormax,
                        status,
                        meth_lvl,
                        trinucleotide,
                    )| {
                        Ok(MethylationSite {
                            chromosome: chromosome.try_into()?,
                            start: location.parse::<u32>()?,
                            end: location.parse::<u32>()? + 1,
                            strand: if (strand == "+") ^ invert_strand {
                                Strand::Sense
                            } else {
                                Strand::Antisense
                            },
                            original: s.to_owned(),
                            context: String::from(context),
                            count_methylated: count_methylated.parse::<u32>()?,
                            count_total: count_total.parse::<u32>()?,
                            posteriormax: posteriormax.parse::<f64>()?,
                            status: status
                                .chars()
                                .next()
                                .ok_or(Error::Simple("Status could not be parsed"))?
                                .into(),
                            meth_lvl: meth_lvl.parse::<f64>()?,
                            context_trinucleotide: String::from(trinucleotide),
                        })
                    },
                )
                .unwrap_or(Err(Error::MethlyationSiteFormat))
        };

        let third_format = |s: &str| {
            s.split('\t')
                .collect_tuple()
                .filter(|(_, _, _, context, _, _, _, _, _, _, _)| context == &"CG")
                .map(
                    |(
                        chromosome,
                        first_nucleotide,
                        second_nucleotide,
                        context,
                        _, // No idea what that is
                        strand,
                        count_methylated,
                        count_total,
                        posteriormax,
                        status,
                        meth_lvl,
                    )| {
                        Ok(MethylationSite {
                            chromosome: chromosome.try_into()?,
                            start: first_nucleotide.parse::<u32>()?,
                            end: second_nucleotide.parse::<u32>()?,
                            strand: if (strand == "+") ^ invert_strand {
                                Strand::Sense
                            } else {
                                Strand::Antisense
                            },
                            original: s.to_owned(),
                            context: String::from(context),
                            count_methylated: count_methylated.parse::<u32>()?,
                            count_total: count_total.parse::<u32>()?,
                            posteriormax: posteriormax.parse::<f64>()?,
                            status: status
                                .chars()
                                .next()
                                .ok_or(Error::Simple("Status could not be parsed"))?
                                .into(),
                            meth_lvl: meth_lvl.parse::<f64>()?,
                            context_trinucleotide: String::from("XXX"),
                        })
                    },
                )
                .unwrap_or(Err(Error::MethlyationSiteFormat))
        };

        // Chromatin State files
        // 1	131800	132400	E10
        let chromatin_state = |s: &str| {
            s.split(['\t', ' '])
                .collect_tuple()
                .map(|(chromosome, start, end, state)| {
                    Ok(MethylationSite {
                        chromosome: chromosome.try_into()?,
                        strand: Strand::Unknown,
                        start: start.parse::<u32>()?,
                        end: end.parse::<u32>()?,
                        context: String::from(state),
                        context_trinucleotide: String::new(),
                        count_methylated: 0,
                        count_total: 0,
                        meth_lvl: 0.0,
                        original: s.to_owned(),
                        posteriormax: 0.0,
                        status: MethylationStatus::U,
                    })
                })
                .unwrap_or(Err(Error::MethlyationSiteFormat))
        };

        // #bedGraph section chr1:1-3480
        // chr1	1	4	1
        let bigwig_format = |s: &str| {
            s.split(['\t', ' '])
                .collect_tuple()
                .map(|(chromosome, start, end, _length)| {
                    Ok(MethylationSite {
                        chromosome: chromosome.try_into()?,
                        strand: Strand::Unknown,
                        start: start.parse::<u32>()?,
                        end: end.parse::<u32>()?,
                        context: String::from("Modification"),
                        context_trinucleotide: String::new(),
                        count_methylated: 0,
                        count_total: 0,
                        meth_lvl: 0.0,
                        original: s.to_owned(),
                        posteriormax: 0.0,
                        status: MethylationStatus::U,
                    })
                })
                .unwrap_or(Err(Error::MethlyationSiteFormat))
        };
        // Heterogeneity Score Files
        // Requested by Patrick Wolf
        // chr | start | end | score | _ | _
        let heterogenity_score_files = |s: &str| {
            s.split(['\t', ' '])
                .collect_tuple()
                .map(|(chromosome, start, end, _score)| {
                    dbg!(chromosome, start, end);
                    Ok(MethylationSite {
                        chromosome: chromosome.try_into()?,
                        start: start.parse()?,
                        end: end.parse()?,
                        ..Default::default()
                    })
                })
                .unwrap_or(Err(Error::MethlyationSiteFormat))
        };

        let results: Result<MethylationSite, crate::error::Error> = first_format(s)
            .or_else(|_| second_format(s))
            .or_else(|_| third_format(s))
            .or_else(|_| chromatin_state(s))
            .or_else(|_| bigwig_format(s))
            .or_else(|_| heterogenity_score_files(s));

        match results {
            Ok(methylation_site) => Some(methylation_site),
            Err(e) => {
                // Log lines that could not be parsed
                // Ignore header
                if !s.contains("seqnames") {
                    print_dev!("Non-fatal error when parsing methylation site: {e}\nLine: {s}");
                }
                None
            }
        }
    }

    /// Checks weather a given CG site belongs to a specific gene. The cutoff is the number of bases upstream and downstream of the gene to consider the CG site in the gene. For example, a cutoff of 1000 would consider a CG site 1000 bases upstream of the gene to be in the gene.
    /// To strictly check weather a CG site is within the gene region, pass a cutoff of 0.
    ///
    /// Passing a negative cutoff is possible but leads to undefined behaviour if used together with ``find_gene``.
    pub fn is_in_gene(&self, gene: &Gene, args: &Args) -> bool {
        let cutoff = if args.cutoff_gene_length {
            gene.end - gene.start
        } else {
            args.cutoff
        };
        self.chromosome == gene.chromosome
            && gene.start <= self.start + cutoff
            && self.end <= gene.end + cutoff
            && self.strand == gene.strand
    }

    /// Find the gene within a genome that a CG site belongs to. Due to binary search, searching is O(log n) where n is the number of genes in the genome.
    /// Therefore, this method is efficient to use on large genomes.
    ///
    /// The lifetime of the genome is longer than the lifetime of the CG site.
    /// GG sites exist only while a single methylation file is being processed but the genome is loaded once and exists for the entire program
    pub fn find_gene<'long>(&self, genome: &'long Genome, args: &Args) -> Option<&'long Gene> {
        // This can fail sensibly if there are chromosomes in the methylome files that are not in the annotation files
        let chromosome = genome.get(&self.chromosome)?;

        let strand: &Vec<Gene> = match self.strand {
            Strand::Sense => &chromosome.sense, // This is a performance hit. Is there a better way to do this?
            Strand::Antisense => &chromosome.antisense,
            Strand::Unknown => &chromosome.combined,
        };

        let search = |gene: &Gene| {
            if args.cutoff_gene_length {
                gene.end + (gene.end - gene.start)
            } else {
                gene.end + args.cutoff
            }
        };

        let first_matching_gene_index = strand
            .binary_search_by_key(&self.start, search)
            .unwrap_or_else(|x| x); // Collapse exact match on gene end and closest previous match into one, as both are valid

        if strand.len() < first_matching_gene_index + 1 {
            return None;
        }

        let gene = &strand[first_matching_gene_index];

        if self.is_in_gene(gene, args) {
            return Some(gene);
        }

        None
    }
    /// Place a CG site in the correct windows. Returns a list of all the successfull insertions as a tuple of the region (upstream, downstream or gene) and the index of the window.
    ///
    /// It works by first finding the region the CG site is in (upstream, downstream or gene) and then finding the windows within that a CG site belongs to.
    /// For genes on the - strand, the windows are reversed, so that the first window is the one closest to the end of the gene.
    pub fn place_in_windows(
        &self,
        gene: &Gene,
        windows: &mut Windows,
        args: &Args,
    ) -> Vec<(Region, usize)> // Return a vector of (strand, window) tuples for each window the CG site is in
    {
        const E: f64 = 0.1; // Epsilon for floating point comparison
        let location = self.start as f64;
        let cutoff = args.cutoff as f64;
        let step = args.window_step as f64;
        let size = args.window_size as f64;
        let start = gene.start as f64;
        let end = gene.end as f64;
        let length = end - start;

        if let Strand::Unknown = self.strand {
            let mut copy = self.clone();
            copy.strand = Strand::Sense;
            return copy.place_in_windows(gene, windows, args);
        }

        // Offset from start for + strand, offset from end for - strand. Can be negative for upstream sites
        let mut windows_in = Vec::new();
        let offset = match &self.strand {
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
        let mut position = match (&region, &self.strand) {
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
                window.push(self.clone());
                windows_in.push((region.clone(), i));
            }
        }
        windows_in
    }
}

impl Display for MethylationSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CG site is located on the {} strand of chromosome {} at bp {}",
            self.strand, self.chromosome, self.start
        )
    }
}

#[cfg(test)]
mod tests {

    use std::fs::read_to_string;

    use super::MethylationSite;
    use crate::{
        arguments::Windows as Args,
        genes::{Gene, Strand},
        methylation_site::Chromosome,
        windows::Windows,
    };

    #[test]
    fn test_instantiate_from_methylome_file_line() {
        let line = "1	23151	+	CG	0	8	0.9999	U	0.0025";
        let cg = MethylationSite::from_methylome_file_line(line, false).unwrap();
        assert_eq!(cg.chromosome, Chromosome::Numbered(1));
    }

    #[test]
    fn test_bigwig_format() {
        let line = "chr1	7	11	3";
        let cg = MethylationSite::from_methylome_file_line(line, true).unwrap();
        assert_eq!(cg.chromosome, Chromosome::Numbered(1));
        assert_eq!(cg.start, 7);
        assert_eq!(cg.end, 11);
    }

    #[test]
    fn test_instantiate_from_methylome_file_line_invalid_line() {
        let line = "1	23151	+	CG	0	8	0.9999	";
        let cg = MethylationSite::from_methylome_file_line(line, false);
        assert!(cg.is_none());
    }

    #[test]
    fn test_cmt3_line_not_cg() {
        let line = "1	25600	+	CHH	0	94	0.9999	U	0.0043	CAT";
        let site = MethylationSite::from_methylome_file_line(line, false);
        assert!(site.is_none());
    }

    #[test]
    fn test_instantiate_from_methylome_file_line_invalid_chromosome() {
        let line = "X	23151	+	CG	0	8	0.9999	U	0.0025";
        let cg = MethylationSite::from_methylome_file_line(line, false);
        assert!(cg.is_none());
    }

    #[test]
    fn test_heterogenity_file_format() {
        let file = read_to_string("data/heterogenity_score_files.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s, false))
            .collect();
        assert_eq!(sites.get(0).unwrap().start, 809);
        assert_eq!(sites.len(), 509825);
    }

    #[test]
    fn test_r_example_methylome_data() {
        let file = read_to_string("data/methylome/G0.txt").unwrap();
        let sites: Vec<_> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s, false))
            .collect();
        assert_eq!(sites.len(), 500); // Only 500 of the > 2800 sites are CG sites

        let file = read_to_string("data/methylome/G1_2.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s, false))
            .collect();
        assert_eq!(sites.len(), 500);

        let file = read_to_string("data/methylome/G4_2.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s, false))
            .collect();
        assert_eq!(sites.len(), 500);

        let file = read_to_string("data/methylome/G4_8.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s, false))
            .collect();
        assert_eq!(sites.len(), 500);
    }

    #[test]
    fn test_extract_gene() {}

    #[test]
    fn test_place_site_absolute() {
        let args = Args {
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
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

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
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

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
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

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

        cg_a.place_in_windows(&gene, &mut windows, &args);
        cg_b.place_in_windows(&gene, &mut windows, &args);
        cg_c.place_in_windows(&gene, &mut windows, &args);
        cg_d.place_in_windows(&gene, &mut windows, &args);
        cg_e.place_in_windows(&gene, &mut windows, &args);
        cg_f.place_in_windows(&gene, &mut windows, &args);
        cg_g.place_in_windows(&gene, &mut windows, &args);
        cg_h.place_in_windows(&gene, &mut windows, &args);

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

        cg_a.place_in_windows(&gene, &mut windows, &args);
        cg_b.place_in_windows(&gene, &mut windows, &args);
        cg_c.place_in_windows(&gene, &mut windows, &args);
        cg_d.place_in_windows(&gene, &mut windows, &args);
        cg_e.place_in_windows(&gene, &mut windows, &args);
        cg_f.place_in_windows(&gene, &mut windows, &args);
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
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

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
            let upstream = cg.place_in_windows(&all_upstream_gene, &mut windows, &args);
            let gene = cg.place_in_windows(&all_within_gene, &mut windows, &args);
            let downstream = cg.place_in_windows(&all_downstream_gene, &mut windows, &args);

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
