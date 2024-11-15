use eyre::bail;
use itertools::Itertools;

#[derive(Debug, Clone)]
pub enum CigarOp {
    Match,
    Mismatch,
    Insertion,
    Deletion,
}

impl TryFrom<char> for CigarOp {
    type Error = eyre::Error;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            '=' => Ok(CigarOp::Match),
            'X' => Ok(CigarOp::Mismatch),
            'I' => Ok(CigarOp::Insertion),
            'D' => Ok(CigarOp::Deletion),
            _ => bail!("Invalid cigar operation ({value})."),
        }
    }
}

pub fn get_cg_ops(cg: &str) -> eyre::Result<Vec<(u32, CigarOp)>> {
    let mut all_cg_ops = vec![];
    let mut bp_elem = None;
    let mut op_elem = None;
    for (is_digit, mut cg_ops) in &cg.chars().chunk_by(|c| c.is_ascii_digit()) {
        if is_digit {
            bp_elem = Some(cg_ops.join("").parse::<u32>()?);
        } else if let Some(cg_op) = cg_ops.next() {
            op_elem = Some(CigarOp::try_from(cg_op)?)
        } else {
            unreachable!()
        }
        if bp_elem.is_some() && op_elem.is_some() {
            all_cg_ops.push((bp_elem.unwrap(), op_elem.unwrap()));
            bp_elem = None;
            op_elem = None;
        }
    }
    Ok(all_cg_ops)
}
