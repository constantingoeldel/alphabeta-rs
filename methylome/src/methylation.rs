use std::fmt::Display;

use crate::{chromosome::Chromosome, strand::Strand, Error, Return};
use itertools::Itertools;

#[macro_export]
macro_rules! print_dev {
    ($($rest:tt)*) => {
        #[cfg(debug_assertions)]
        std::println!($($rest)*)
    }
}

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
    pub fn from_methylome_file_line(s: &str) -> Option<Self> {
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
                            strand: if strand == "+" {
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
                            strand: if strand == "+" {
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
                            strand: if strand == "+" {
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

        let results: Return<MethylationSite> = first_format(s)
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
    use crate::methylation::Chromosome;

    #[test]
    fn test_instantiate_from_methylome_file_line() {
        let line = "1	23151	+	CG	0	8	0.9999	U	0.0025";
        let cg = MethylationSite::from_methylome_file_line(line).unwrap();
        assert_eq!(cg.chromosome, Chromosome::Numbered(1));
    }

    #[test]
    fn test_bigwig_format() {
        let line = "chr1	7	11	3";
        let cg = MethylationSite::from_methylome_file_line(line).unwrap();
        assert_eq!(cg.chromosome, Chromosome::Numbered(1));
        assert_eq!(cg.start, 7);
        assert_eq!(cg.end, 11);
    }

    #[test]
    fn test_instantiate_from_methylome_file_line_invalid_line() {
        let line = "1	23151	+	CG	0	8	0.9999	";
        let cg = MethylationSite::from_methylome_file_line(line);
        assert!(cg.is_none());
    }

    #[test]
    fn test_cmt3_line_not_cg() {
        let line = "1	25600	+	CHH	0	94	0.9999	U	0.0043	CAT";
        let site = MethylationSite::from_methylome_file_line(line);
        assert!(site.is_none());
    }

    #[test]
    fn test_instantiate_from_methylome_file_line_invalid_chromosome() {
        let line = "X	23151	+	CG	0	8	0.9999	U	0.0025";
        let cg = MethylationSite::from_methylome_file_line(line);
        assert!(cg.is_none());
    }

    #[test]
    fn test_heterogenity_file_format() {
        let file = read_to_string("data/heterogenity_score_files.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s))
            .collect();
        assert_eq!(sites.get(0).unwrap().start, 809);
        assert_eq!(sites.len(), 509825);
    }

    #[test]
    fn test_r_example_methylome_data() {
        let file = read_to_string("data/methylome/G0.txt").unwrap();
        let sites: Vec<_> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s))
            .collect();
        assert_eq!(sites.len(), 500); // Only 500 of the > 2800 sites are CG sites

        let file = read_to_string("data/methylome/G1_2.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s))
            .collect();
        assert_eq!(sites.len(), 500);

        let file = read_to_string("data/methylome/G4_2.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s))
            .collect();
        assert_eq!(sites.len(), 500);

        let file = read_to_string("data/methylome/G4_8.txt").unwrap();
        let sites: Vec<MethylationSite> = file
            .split('\n')
            .filter_map(|s| MethylationSite::from_methylome_file_line(s))
            .collect();
        assert_eq!(sites.len(), 500);
    }

    #[test]
    fn test_extract_gene() {}
}
