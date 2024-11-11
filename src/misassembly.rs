use std::str::FromStr;

use coitrees::{GenericInterval, IntervalNode};
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

pub fn get_misassembly_from_itv(
    itv: Option<IntervalNode<Option<String>, usize>>,
) -> Option<Misassembly> {
    itv.as_ref()
        .and_then(|m| {
            m.metadata()
                .as_ref()
                .and_then(|m| Misassembly::from_str(m).ok())
        })
        // Filter HETs that aren't serious misassemblies.
        .filter(|m| *m != Misassembly::HET)
}
