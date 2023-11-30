use std::{collections::HashMap, fmt::Display, ops::Deref};

use itertools::Itertools;
use methylome::{Chromosome, MethylationSite, Strand};

use crate::{Error, Return};

use crate::print_dev;

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

        let results: Option<Return<Gene>> = first_format(s).or_else(|| second_format(s)).or(None);

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

/// Checks weather a given CG site belongs to a specific gene. The cutoff is the number of bases upstream and downstream of the gene to consider the CG site in the gene. For example, a cutoff of 1000 would consider a CG site 1000 bases upstream of the gene to be in the gene.
/// To strictly check weather a CG site is within the gene region, pass a cutoff of 0.
///
/// Passing a negative cutoff is possible but leads to undefined behaviour if used together with ``find_gene``.
pub fn is_in_gene(
    site: &MethylationSite,
    gene: &Gene,
    cutoff_gene_length: bool,
    cutoff: u32,
) -> bool {
    let cutoff = if cutoff_gene_length {
        gene.end - gene.start
    } else {
        cutoff
    };
    site.chromosome == gene.chromosome
        && gene.start <= site.start + cutoff
        && site.end <= gene.end + cutoff
        && site.strand == gene.strand
}

/// Find the gene within a genome that a CG site belongs to. Due to binary search, searching is O(log n) where n is the number of genes in the genome.
/// Therefore, this method is efficient to use on large genomes.
///
/// The lifetime of the genome is longer than the lifetime of the CG site.
/// GG sites exist only while a single methylation file is being processed but the genome is loaded once and exists for the entire program
pub fn find_gene<'long>(
    site: &MethylationSite,
    genome: &'long Genome,
    cutoff_gene_length: bool,
    cutoff: u32,
) -> Option<&'long Gene> {
    // This can fail sensibly if there are chromosomes in the methylome files that are not in the annotation files
    let chromosome = genome.get(&site.chromosome)?;

    let strand: &Vec<Gene> = match site.strand {
        Strand::Sense => &chromosome.sense, // This is a performance hit. Is there a better way to do this?
        Strand::Antisense => &chromosome.antisense,
        Strand::Unknown => &chromosome.combined,
    };

    let search = |gene: &Gene| {
        if cutoff_gene_length {
            gene.end + (gene.end - gene.start)
        } else {
            gene.end + cutoff
        }
    };

    let first_matching_gene_index = strand
        .binary_search_by_key(&site.start, search)
        .unwrap_or_else(|x| x); // Collapse exact match on gene end and closest previous match into one, as both are valid

    if strand.len() < first_matching_gene_index + 1 {
        return None;
    }

    let gene = &strand[first_matching_gene_index];

    if is_in_gene(site, gene, cutoff_gene_length, cutoff) {
        return Some(gene);
    }

    None
}
