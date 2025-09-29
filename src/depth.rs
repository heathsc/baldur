use std::{collections::HashMap, fmt};

use crate::{alleles::InsAlleles, detailed_align::DetailedAlign, model::N_QUAL};

// BDHU^J are low quality ACGT-I, bdhu&j are low quality acgt_i
// S = skip
// - = del
// I = insertion
pub const BASES: &str = "ACGT-ISNBDHU^Jacgt_isnbdhu&j";
pub const BASES_LEN: usize = BASES.len();
pub const HIDX: usize = BASES_LEN >> 1;
pub const BASES_LIM: usize = 5; // Index of last 'base' (I)
pub const BASES_INS: usize = 5;
pub const BASES_SKIP: usize = 6; // Index of skip
pub const BASES_LOWQ: usize = 8; // Index of start of low qual bases

#[derive(Clone, Default)]
pub(crate) struct DepthCounts {
    pub(crate) cts: Vec<[usize; 2]>,
    pub(crate) qcts: Vec<[usize; N_QUAL]>,
}

impl DepthCounts {
    pub(crate) fn new(n_alls: usize) -> Self {
        Self {
            cts: vec![[0; 2]; n_alls],
            qcts: vec![[0; N_QUAL]; n_alls],
        }
    }
}

impl fmt::Display for DepthCounts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let cts = &self.cts;
        let s3 = cts.iter().fold(0, |s, x| s + x[0] + x[1]);
        let l = cts.len();
        let s2 = if l > 3 {
            let x = &cts[l - 3];
            s3 - x[0] - x[1]
        } else {
            s3
        };
        let s1 = if l > 1 {
            let x = &cts[l - 1];
            s2 - x[0] - x[1]
        } else {
            s2
        };
        write!(f, "{}\t{}\t{}", s3, s2, s1)?;
        // for x in cts.iter() {
        //    write!(f, "\t{}\t{}", x[0], x[1])?;
        //}
        Ok(())
    }
}

pub(crate) struct Depth {
    hash: HashMap<u8, usize>,
    pub(crate) counts: Vec<DepthCounts>,
    pub(crate) dalign: Option<Vec<DetailedAlign>>,
    pub(crate) ins_hash: HashMap<usize, InsAlleles>,
}

pub fn base_to_ix_strand(ix: usize) -> (usize, usize) {
    const HM1: usize = HIDX - 1;
    const HLIM: usize = HIDX + BASES_LIM;
    const HSKIP: usize = HIDX + BASES_SKIP;
    match ix {
        0..=BASES_LIM => (ix, 0),      // Forward A,C,G,T,Del,I
        BASES_SKIP => (4, 0),          // Forward skip counts as Del
        BASES_LOWQ..=HM1 => (6, 0),    // Forward low quality bases
        HIDX..=HLIM => (ix - HIDX, 1), // Reverse A,C,G,T,Del
        HSKIP => (4, 1),               // Reverse skip
        _ => (6, 1),                   // Reverse low quality bases
    }
}

impl Depth {
    pub(crate) fn new(size: usize, detail: bool) -> Self {
        let hash: HashMap<u8, usize> = BASES.bytes().enumerate().map(|(i, c)| (c, i)).collect();
        let counts = vec![DepthCounts::new(BASES_LIM + 2); size];
        let dalign = if detail {
            Some(Vec::with_capacity(1024))
        } else {
            None
        };
        let ins_hash = HashMap::new();
        Self {
            hash,
            counts,
            dalign,
            ins_hash,
        }
    }

    pub(crate) fn add_obs_vec(
        &mut self,
        obs: &[u8],
        qual: &[u8],
        ins_hash: HashMap<usize, (usize, u8)>,
        id: &str,
    ) {
        assert_eq!(obs.len(), self.counts.len());
        if let Some(d) = self.dalign.as_mut()
            && let Some(dall) = DetailedAlign::from_obs_vec(obs, qual, ins_hash, &self.hash, id)
        {
            d.push(dall)
        }
        let itr = obs.iter().zip(qual.iter().map(|x| *x as usize));
        for ((o, q), ct) in itr.zip(self.counts.iter_mut()) {
            if let Some(&ix) = self.hash.get(o) {
                let (k, s) = base_to_ix_strand(ix);
                ct.cts[k][s] += 1;
                ct.qcts[k][q] += 1;
            }
        }
    }
}
