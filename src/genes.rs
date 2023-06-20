use std::{collections::HashMap, fmt::Display, ops::Deref};

use itertools::Itertools;

use crate::{
    error::{self, Error},
    methylation_site::Chromosome,
    print_dev,
};

pub type Result<T> = std::result::Result<T, error::Error>;

#[derive(Clone, Debug)]
pub enum Strand {
    Sense,
    Antisense,
    Unknown,
}

// TODO Can we somehow get rid of this? Cloning Genome is expensive
// Clone is needed for try_for_each_with of the rayon lib for processing the genes
// Although at the cloning moment, it is still empty, so not expensive
#[derive(Clone)]
pub struct Genome(HashMap<Chromosome, GenesByStrand>);

impl Deref for Genome {
    type Target = HashMap<Chromosome, GenesByStrand>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

// impl DerefMut for Genome {
//     fn deref_mut(&mut self) -> &mut Self::Target {
//         &mut self.0
//     }
// }

impl Genome {
    pub fn new(chromosomes: Vec<&Chromosome>) -> Self {
        let mut genome = HashMap::new();
        for chromosome in chromosomes {
            genome.insert(chromosome.to_owned(), GenesByStrand::new());
        }
        Genome(genome)
    }

    pub fn insert_gene(&mut self, chromsome: Chromosome, gene: Gene) {
        self.0.get_mut(&chromsome).unwrap().insert(gene);
    }

    pub fn sort(&mut self) {
        for genes_by_strand in self.0.values_mut() {
            genes_by_strand.sort();
        }
    }
}

impl<'a> TryFrom<&'a str> for Strand {
    type Error = Error;
    fn try_from(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Self::Sense),
            "-" => Ok(Self::Antisense),
            "*" => Ok(Self::Unknown),
            _ => Err(Error::Simple("Invalid strand")),
        }
    }
}

impl Strand {
    pub fn correct_format(s: &str) -> bool {
        matches!(s, "+" | "-" | "*")
    }

    pub fn invert(self, invert: bool) -> Self {
        if !invert {
            return self;
        }
        match self {
            Self::Sense => Self::Antisense,
            Self::Antisense => Self::Sense,
            Self::Unknown => Self::Unknown,
        }
    }
}

impl PartialEq for Strand {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Sense, Self::Antisense) => false,
            (Self::Antisense, Self::Sense) => false,
            _ => true, // If one of them is unknown, they are equal. WARNING: This leads to unexpected behaviour if you want to check for one specific strand! in that case use `if let Strand::X = self.strand {}
        }
    }
}

impl Eq for Strand {}

#[derive(Debug, Clone)]
pub enum Region {
    Upstream,
    Gene,
    Downstream,
}

impl Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Upstream => write!(f, "upstream"),
            Self::Gene => write!(f, "gene"),
            Self::Downstream => write!(f, "downstream"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gene {
    pub chromosome: Chromosome,
    pub start: u32,
    pub end: u32,
    pub name: String,
    pub annotation: String,
    pub strand: Strand,
}

#[derive(Clone)]
pub struct GenesByStrand {
    pub sense: Vec<Gene>,
    pub antisense: Vec<Gene>,
    pub combined: Vec<Gene>,
}

impl Default for GenesByStrand {
    fn default() -> Self {
        GenesByStrand::new()
    }
}

impl GenesByStrand {
    pub fn new() -> Self {
        GenesByStrand {
            sense: Vec::new(),
            antisense: Vec::new(),
            combined: Vec::new(),
        }
    }

    pub fn insert(&mut self, gene: Gene) {
        self.combined.push(gene.clone());
        match gene.strand {
            Strand::Sense => self.sense.push(gene),
            Strand::Antisense => self.antisense.push(gene),
            _ => (),
        }
    }

    pub fn sort(&mut self) {
        self.sense.sort_by(|a, b| a.start.cmp(&b.start));
        self.antisense.sort_by(|a, b| a.start.cmp(&b.start));
        self.combined.sort_by(|a, b| a.start.cmp(&b.start));
    }
}

impl Gene {
    pub fn from_annotation_file_line(s: &str, invert_strand: bool) -> Option<Self> {
        let first_format = |s: &str| {
            s.split([' ', '\t'])
                .collect_tuple()
                .filter(|(_, _, _, _, _, strand)| Strand::correct_format(strand))
                .map(|(chromosome, start, end, name, annotation, strand)| {
                    let strand: Strand = strand.try_into()?;
                    Ok(Gene {
                        chromosome: chromosome.try_into()?,
                        start: start.parse::<u32>()?,
                        end: end.parse::<u32>()?,
                        name: String::from(name),
                        annotation: String::from(annotation),
                        strand: strand.invert(invert_strand),
                    })
                })
        };

        // seqnames	start	end	width	strand	gbM.id
        // 1	6362000	6362200	201	*	AT1G18480
        let second_format = |s: &str| {
            s.split([' ', '\t'])
                .collect_tuple()
                .filter(|(_, _, _, _, strand, _)| Strand::correct_format(strand))
                .map(|(chromosome, start, end, _width, strand, name)| {
                    let strand: Strand = strand.try_into()?;
                    Ok(Gene {
                        chromosome: chromosome.try_into()?,
                        start: start.parse::<u32>()?,
                        end: end.parse::<u32>()?,
                        name: String::from(name),
                        annotation: String::new(),
                        strand: strand.invert(invert_strand),
                    })
                })
        };

        let results: Option<Result<Gene>> = first_format(s).or_else(|| second_format(s)).or(None);

        match results {
            Some(Ok(methylation_site)) => Some(methylation_site),
            Some(Err(e)) => {
                print_dev!("Non-fatal error when parsing gene: {e}\nLine: {s}");
                None
            }
            None => {
                print_dev!("Could not find a suitable parser for line: {s}");
                None
            }
        }
    }
}

impl Display for Gene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Gene {} is located on the {} strand of chromosome {} ranging from bp {} to bp {} ",
            self.name, self.strand, self.chromosome, self.start, self.end
        )
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Sense => write!(f, "+"),
            Strand::Antisense => write!(f, "-"),
            Strand::Unknown => write!(f, "*"),
        }
    }
}

impl Display for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Chromosome::Mitochondrial => write!(f, "M"),
            Chromosome::Chloroplast => write!(f, "C"),
            Chromosome::Numbered(n) => write!(f, "{n}"),
        }
    }
}
