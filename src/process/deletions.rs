use std::{
    cmp::Ordering,
    collections::BTreeMap,
    fmt::{self, Formatter},
};

use rstar::RTree;

mod del_prob;
mod del_set;
mod rtree;

use del_prob::LikeContrib;
use del_set::{DelSet, ReadDels};

use crate::process::deletions::del_prob::calc_del_prob;

#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
pub enum DelType {
    Del,
    Split,
    Wildtype
}

impl fmt::Display for DelType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Del => "Deletion",
                Self::Split => "Split",
                Self::Wildtype => "Wildtype",
            }
        )
    }
}

/// A deletion found within a read.  The sign of `size` indicates the
/// the strand.  Note that because we are dealing with circular genomes
/// we can not infer the strand from the start and end positions alone
#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
pub struct Deletion {
    coords: DelCoords,
    size: isize,
    dtype: DelType,
    ix: u64,
    ct: [usize; 2], // Observatons of this deletion on the top and bottom strands
}

impl fmt::Display for Deletion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.coords.start,
            self.coords.end,
            self.size,
            self.dtype,
            self.ix,
            self.ct[0],
            self.ct[1]
        )
    }
}

#[derive(Debug)]
pub struct ReadExtent {
    coords: DelCoords,
    target_size: isize,
    reversed: bool,
    obs_deletion: Option<u64>,
}

impl ReadExtent {
    fn new(
        start: isize,
        end: isize,
        reversed: bool,
        obs_deletion: Option<u64>,
        target_size: isize,
    ) -> Self {
        assert!(start < end);
        Self {
            coords: DelCoords::new(start, end),
            reversed,
            obs_deletion,
            target_size,
        }
    }

    #[inline]
    fn coords(&self) -> (isize, isize) {
        self.coords.coords()
    }

    fn physical_start(&self) -> isize {
        let s = if self.reversed {
            self.coords.end
        } else {
            self.coords.start
        };

        s % self.target_size
    }

    #[inline]
    fn coords2(&self) -> Option<(isize, isize)> {
        if self.coords.end > self.target_size {
            Some((
                self.coords.start - self.target_size,
                self.coords.end - self.target_size,
            ))
        } else {
            None
        }
    }

    #[inline]
    fn obs_deletion(&self) -> Option<u64> {
        self.obs_deletion
    }
}

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq)]
pub struct DelCoords {
    start: isize,
    end: isize,
}

impl DelCoords {
    #[inline]
    pub fn new(start: isize, end: isize) -> Self {
        Self { start, end }
    }

    #[inline]
    pub fn coords(&self) -> (isize, isize) {
        (self.start, self.end)
    }
}

impl PartialOrd for DelCoords {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DelCoords {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Equal => self.end.cmp(&other.end),
            c => c,
        }
    }
}

/// We store deletions at least as large as min_size in a hash table so we can get summaries of
/// large deletions across all reads
pub struct Deletions {
    del_hash: Option<BTreeMap<DelCoords, Deletion>>,
    read_extents: Option<Vec<ReadExtent>>,
    min_size: isize,
    target_size: usize,
    adjust: isize,
}

impl Deletions {
    pub fn new(target_size: usize, min_size: usize, adjust: usize) -> Self {
        Self {
            del_hash: Some(BTreeMap::new()),
            read_extents: Some(Vec::new()),
            min_size: min_size as isize,
            adjust: adjust as isize,
            target_size,
        }
    }

    pub fn add_del(
        &mut self,
        start: usize,
        end: usize,
        reversed: bool,
        dtype: DelType,
    ) -> Option<u64> {
        let s = start as isize - self.adjust;
        let e = end as isize - self.adjust;

        let size = e - s + if s <= e { 1 } else { -1 };

        let adj = |x: isize| x % (self.target_size as isize);

        if size.abs() >= self.min_size {
            let (start, end) = if reversed {
                (adj(e), adj(s))
            } else {
                (adj(s), adj(e))
            };

            let dc = if start < end {
                DelCoords::new(start, end)
            } else {
                DelCoords::new(start, end + self.target_size as isize)
            };

            let dh = self.del_hash.as_mut().unwrap();
            let l: u64 = dh.len().try_into().unwrap();
            let del = dh.entry(dc).or_insert_with(|| Deletion {
                coords: dc,
                size,
                dtype,
                ix: l,
                ct: [0; 2],
            });

            if reversed {
                del.ct[1] += 1
            } else {
                del.ct[0] += 1
            }
            Some(del.ix)
        } else {
            None
        }
    }

    pub fn add_read_extent(
        &mut self,
        start: usize,
        end: usize,
        reversed: bool,
        del_idx: Option<u64>,
    ) {
        assert!(start < end);
        // Adjust coodinates so that start is in first genome copy if necessary
        let s = start % self.target_size;
        // Apply same adjustment to end coordinat
        let e = end - (start - s);

        self.read_extents.as_mut().unwrap().push(ReadExtent::new(
            s as isize - self.adjust,
            e as isize - self.adjust,
            reversed,
            del_idx,
            self.target_size as isize,
        ));
    }

    pub fn target_size(&self) -> usize {
        self.target_size
    }

    pub fn len(&self) -> usize {
        self.del_hash
            .as_ref()
            .map(|dh| dh.len())
            .unwrap_or_default()
    }

}

pub struct DeletionWork {
    del_tree: RTree<Deletion>,
    del_freq: Vec<f64>,
    target_size: usize,
}

impl DeletionWork {
    pub fn new(mut dels: Deletions) -> anyhow::Result<Option<Self>> {
        let n_dels = dels.len();
        if n_dels == 0 {
            Ok(None)
        } else {
            let target_size = dels.target_size();
            
            // Get likelihood contributions from reads
            let (lc, del_tree) = dels
                .get_like_data()
                .expect("Problem getting deletion contributions");

            // Estimate deletion frequencies
            let del_freq = lc.est_freq()?;
            Ok(Some(Self { del_tree, del_freq, target_size }))
        }
    }
    
    pub fn write_deletions(&self, prefix: &str) -> anyhow::Result<()> {
        rtree::write_deletions(&self.del_tree, prefix, &self.del_freq)
    }

    pub fn get_del_prob(&mut self, x: usize) -> f64 {
        calc_del_prob(x, &self.del_tree, &self.del_freq, self.target_size)
    }
}
