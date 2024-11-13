use std::str::FromStr;

use eyre::bail;

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, PartialEq, Eq)]
pub enum Misassembly {
    MISJOIN,
    GAP,
    COLLAPSE,
    ERROR,
    HET,
}

impl FromStr for Misassembly {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "MISJOIN" => Misassembly::MISJOIN,
            "GAP" => Misassembly::GAP,
            "COLLAPSE" | "COLLAPSE_VAR" => Misassembly::COLLAPSE,
            "ERROR" => Misassembly::ERROR,
            "HET" => Misassembly::HET,
            _ => bail!("Invalid misassembly type. ({s})"),
        })
    }
}
