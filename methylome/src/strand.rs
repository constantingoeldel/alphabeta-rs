use std::fmt::Display;

use crate::{Error, Return};

#[derive(Clone, Debug)]
pub enum Strand {
    Sense,
    Antisense,
    Unknown,
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

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Sense => write!(f, "+"),
            Strand::Antisense => write!(f, "-"),
            Strand::Unknown => write!(f, "*"),
        }
    }
}

impl<'a> TryFrom<&'a str> for Strand {
    type Error = Error;
    fn try_from(s: &str) -> Return<Self> {
        match s {
            "+" => Ok(Self::Sense),
            "-" => Ok(Self::Antisense),
            "*" => Ok(Self::Unknown),
            _ => Err(Error::Simple("Invalid strand")),
        }
    }
}
