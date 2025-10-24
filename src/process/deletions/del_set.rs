use std::{
    collections::{HashMap, hash_map},
    num::NonZeroUsize,
};

/// Store set of deletions indexed by an integer from 0..n
/// where n is the total number of deletions
#[derive(Debug, Hash, PartialEq, Eq)]
pub struct DelSet {
    m: Box<[u64]>,
    n_dels: usize,
}

impl DelSet {
    pub(super) fn new(n_dels: usize, n: usize) -> Self {
        let n_u64 = calc_n_u64(n_dels);
        let m = vec![0u64; n_u64 * n].into_boxed_slice();
        Self { m, n_dels }
    }
    
    pub(super) fn all_dset(&mut self) -> &mut [u64] {
        &mut self.m
    }
    
    pub(super) fn add_index(&mut self, i: usize, k: usize) {
        let n_u64 = calc_n_u64(self.n_dels);
        let (j, x) = get_coords(i);
        self.m[j + n_u64 * k] |= x
    }
    
    #[inline]
    pub fn dset(&self, k: usize) -> &[u64] {
        let n_u64 = calc_n_u64(self.n_dels);
        &self.m[k * n_u64 .. (k + 1) * n_u64]
    }

    pub fn mask(n: usize) -> Self {
        assert!(n > 0);
        let k = calc_n_u64(n);
        let mut ds = Self::new(n, 1);
        for m in ds.m[..k - 1].iter_mut() {
            *m = u64::MAX
        }
        ds.m[k - 1] = u64::MAX >> (63 - ((n - 1) & 63));
        ds
    }
}

#[inline]
pub(super) fn calc_n_u64(n: usize) -> usize {
    (n + 63) >> 6
}

#[inline]
pub(super) fn get_coords(i: usize) -> (usize, u64) {
    (i >> 6, 1u64 << (i & 63))
}

/// Summary of a read w.r.t. deletions
/// covered is a bitmap of which deletions are covered by the read
/// while obs_index has the index of the observed deletion (if present)
#[derive(Debug, Hash, PartialEq, Eq)]
pub struct ReadDelSet {
    dset: DelSet,
    obs_index: Option<NonZeroUsize>,
}

impl ReadDelSet {
    #[inline]
    pub fn obs_index(&self) -> Option<usize> {
        self.obs_index.map(|i| i.get() - 1)
    }

    pub fn new(dset: DelSet, i: Option<u64>) -> Self {
        let obs_index = i.and_then(|x| NonZeroUsize::new(x as usize + 1));
        Self { dset, obs_index }
    }

    #[inline]
    pub fn dset(&self, k: usize) -> &[u64] {
        self.dset.dset(k)
    }
}

pub struct ReadDels {
    hash: HashMap<ReadDelSet, u64>,
    n_dels: usize,
}

impl ReadDels {
    pub(super) fn new(n_dels: usize) -> Self {
        Self {
            hash: HashMap::new(),
            n_dels,
        }
    }

    pub(super) fn add(&mut self, covered: DelSet, obs_index: Option<u64>) {
        let ds = ReadDelSet::new(covered, obs_index);
        *self.hash.entry(ds).or_default() += 1;
    }

    pub(super) fn iter(&self) -> hash_map::Iter<'_, ReadDelSet, u64> {
        self.hash.iter()
    }

    pub(super) fn del_mask(&self) -> DelSet {
        DelSet::mask(self.n_dels)
    }
    
    pub(super) fn n_dels(&self) -> usize {
        self.n_dels
    }
}
