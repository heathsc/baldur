use regex::{Match, Regex};

use crate::reference::Reference;

pub struct Region {
    tid: usize,
    start: usize,
    stop: usize,
    ctg_size: usize,
}

impl Region {
    pub fn from_str(reg_str: &str, reference: &Reference) -> anyhow::Result<Self> {
        let err = |s| Err(anyhow!("Could not parse region string '{}'", s,));

        let parse_x = |s: Match| parse_usize_with_commas(s.as_str());

        let reg = Regex::new(r#"^([^:]+):?([0-9,]+)?-?([0-9,]+)?"#).unwrap();
        if let Some(cap) = reg.captures(reg_str) {
            match (cap.get(1), cap.get(2), cap.get(3)) {
                (Some(c), None, None) => Self::new(c.as_str(), None, None, reference),
                (Some(c), Some(p), None) => Self::new(c.as_str(), parse_x(p), None, reference),
                (Some(c), Some(p), Some(q)) => {
                    Self::new(c.as_str(), parse_x(p), parse_x(q), reference)
                }
                _ => err(reg_str),
            }
        } else {
            err(reg_str)
        }
    }

    pub fn new(
        chrom: &str,
        start: Option<usize>,
        stop: Option<usize>,
        reference: &Reference,
    ) -> anyhow::Result<Self> {
        let ctg = reference
            .name2contig(chrom)
            .ok_or_else(|| anyhow!("Contig {} not present in input file", chrom))?;
        let start = start.unwrap_or(1).max(1) - 1;
        let ctg_size = ctg.size();
        assert!(ctg_size > 0, "Zero sized contig!");
        let stop = stop.unwrap_or(ctg_size).max(1).min(ctg_size) - 1;
        debug!("Chromosome region: {}:{}-{}", ctg.name(), start, stop);
        if stop >= start {
            Ok(Region {
                tid: ctg.tid(),
                start,
                stop,
                ctg_size,
            })
        } else {
            Err(anyhow!("Invalid range - stop < start".to_string(),))
        }
    }

    pub fn tid(&self) -> usize {
        self.tid
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn len(&self) -> usize {
        self.stop + 1 - self.start
    }

    pub fn ctg_size(&self) -> usize {
        self.ctg_size
    }
}


fn parse_usize_with_commas(s: &str) -> Option<usize> {
    s.replace(',', "").parse::<usize>().ok()
}
