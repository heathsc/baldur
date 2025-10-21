use crate::{cli::Config, context::N_CTXT, depth::Depth, reference::RefPos};

use super::deletions::Deletions;

type Qhist = [[usize; 4]; 64];

pub struct ProcWork<'a> {
    pub ref_seq: &'a [RefPos],
    pub depth: Depth,
    pub qual_hist: Qhist,
    pub ctxt_hist: [Qhist; N_CTXT],
    pub dels: Option<Deletions<'a>>,
}

impl<'a> ProcWork<'a> {
    pub fn new(cfg: &'a Config) -> anyhow::Result<Self> {
        let reg = cfg.region();
        let ref_seq = cfg.reference().contig(reg.tid()).unwrap().seq();

        let dels = if cfg.output_deletions() {
            Some(Deletions::new(
                cfg.region().ctg_size(),
                cfg.large_deletion_limit(),
                cfg.adjust(),
            ))
        } else {
            None
        };

        Ok(Self {
            depth: Depth::new(reg.len(), !cfg.no_call()),
            qual_hist: [[0; 4]; 64],
            ctxt_hist: [[[0; 4]; 64]; N_CTXT],
            ref_seq,
            dels,
        })
    }
}
