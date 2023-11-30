use std::fmt::Display;

use crate::Error;

// TODO: Is there a better name for the regular chromosomes than "Numbered"?
#[derive(Clone, PartialEq, Debug, Eq, Hash)]
pub enum Chromosome {
    Numbered(u8),
    Mitochondrial,
    Chloroplast,
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
