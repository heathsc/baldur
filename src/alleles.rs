use std::{
    cmp::Ordering,
    collections::HashMap,
    fmt,
    ops::{Deref, DerefMut},
    rc::Rc,
};

use log::Level::Trace;

use crate::vcf::VcfRes;
use crate::{
    cli::{Config, ThresholdType},
    depth::HIDX,
    depth::*,
    detailed_align::DetailedAlign,
    model::N_QUAL,
};

// Inserted allele counts
pub(crate) struct InsAllele {
    pub(crate) allele: Rc<[u8]>,
    pub(crate) cts: [usize; 2],
    pub(crate) qcts: [usize; N_QUAL],
}

impl InsAllele {
    pub(crate) fn new(allele: Rc<[u8]>) -> Self {
        Self {
            allele,
            cts: [0; 2],
            qcts: [0; N_QUAL],
        }
    }
}

#[derive(Default)]
pub(crate) struct InsAlleles {
    pub(crate) hash: HashMap<Rc<[u8]>, usize>,
    pub(crate) alleles: Vec<InsAllele>,
}

impl InsAlleles {
    pub(crate) fn get_or_insert(&mut self, all: &[u8]) -> usize {
        if let Some(ix) = self.hash.get(all) {
            *ix
        } else {
            let ix = self.alleles.len();
            let all = Rc::from(all.to_vec().into_boxed_slice());
            self.hash.insert(Rc::clone(&all), ix);
            self.alleles.push(InsAllele::new(all));
            ix
        }
    }
    pub(crate) fn add_obs(&mut self, ix: usize, q: u8, reverse: bool) {
        let k = if reverse { 1 } else { 0 };
        let all = &mut self.alleles[ix];
        all.cts[k] += 1;
        all.qcts[q as usize] += 1;
    }
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, Debug)]
pub enum Trunc {
    No,
    Left(usize),
    Right,
    Both(usize),
}

impl Default for Trunc {
    fn default() -> Self {
        Self::No
    }
}

#[derive(Clone, Hash, PartialEq, Eq)]
pub struct AllDesc {
    v: Vec<u8>,
    ix: Option<usize>,
    trunc: Trunc,
}

impl AllDesc {
    pub fn new(sz: usize) -> Self {
        Self {
            v: Vec::with_capacity(sz),
            ix: None,
            trunc: Trunc::No,
        }
    }

    pub fn make(v: Vec<u8>, ix: Option<usize>, trunc: Trunc) -> Self {
        Self { v, ix, trunc }
    }

    pub fn make_rev(&self) -> Self {
        Self {
            v: self.v.iter().map(|&c| c + HIDX as u8).collect(),
            ix: self.ix,
            trunc: self.trunc,
        }
    }

    pub fn push(&mut self, c: u8) {
        self.v.push(c)
    }

    pub fn clear(&mut self) {
        self.v.clear()
    }

    pub fn set_trunc(&mut self, trunc: Trunc) {
        self.trunc = trunc
    }

    pub fn start(&self) -> usize {
        match self.trunc {
            Trunc::Left(x) | Trunc::Both(x) => x,
            _ => 0,
        }
    }
}

impl Deref for AllDesc {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.v
    }
}

impl DerefMut for AllDesc {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.v
    }
}

impl fmt::Display for AllDesc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(f, "{}\t", self.start())?;
            for i in self.iter().map(|&x| x as usize) {
                write!(f, "{}", BASES.as_bytes()[i] as char)?
            }
        } else {
            for i in self.iter().map(|&x| x as usize) {
                let i1 = i % HIDX;
                if i1 != 4 && i1 != BASES_SKIP {
                    write!(f, "{}", BASES.as_bytes()[i] as char)?
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn get_allele_counts(
    adesc: &[AllDesc],
    x1: usize,
    x: usize,
    ref_len: usize,
    dalign: &[DetailedAlign],
    ins_hash: &HashMap<usize, InsAlleles>,
) -> (Vec<[usize; 2]>, Vec<[usize; N_QUAL]>, Vec<bool>) {
    let ref_all_len = adesc[0].len();
    let indel_flag: Vec<_> = adesc
        .iter()
        .map(|all| all.len() != ref_all_len || all.contains(&4))
        .collect();

    let mut all_hash = HashMap::with_capacity(2 * adesc.len());
    for (ix, ds) in adesc.iter().enumerate() {
        let ds1 = ds.clone();
        // Forward allele
        all_hash.insert(ds1, (ix, 0));
        let ds2 = ds.make_rev();
        // Reverse allele
        all_hash.insert(ds2, (ix, 1));
    }

    let mut obs_all = AllDesc::new(x + 1 - x1);
    let n_alls = adesc.len();
    let mut cts: Vec<[usize; 2]> = vec![[0; 2]; n_alls];
    let mut qcts: Vec<[usize; N_QUAL]> = vec![[0; N_QUAL]; n_alls];
    for dr in dalign.iter() {
        let q = dr.extract_allele(x1, x, ref_len, &mut obs_all, ins_hash);
        if obs_all.is_empty() {
            obs_all.push(if dr.reverse() { 4 + HIDX as u8 } else { 4 });
        }
        if q < N_QUAL as u8
            && let Some((ix, str)) = all_hash.get(&obs_all)
        {
            cts[*ix][*str] += 1;
            qcts[*ix][q as usize] += 1;
        }
    }
    (cts, qcts, indel_flag)
}

/// Find start and length of longest continuous deletion in allele
fn longest_del(all: &AllDesc) -> Option<(usize, usize)> {
    let mut del = None;
    let mut best: Option<(usize, usize)> = None;
    let mut state = false;
    for (i, &c) in all.iter().enumerate() {
        if state {
            let (j, k) = del.take().unwrap();
            if c == 4 {
                del = Some((j, k + 1))
            } else {
                state = false;
                if let Some((_, k1)) = best.as_ref() {
                    if k > *k1 {
                        best = Some((j, k))
                    }
                } else {
                    best = Some((j, k))
                }
            }
        } else if c == 4 {
            state = true;
            del = Some((i, 1))
        }
    }
    if let Some((j, k)) = del.take() {
        if let Some((_, k1)) = best.as_ref() {
            if k > *k1 {
                best = Some((j, k))
            }
        } else {
            best = Some((j, k))
        }
    }
    best
}

#[derive(Default, Clone, Debug)]
struct DelAllele {
    x: usize,
    len: usize,
    v: Vec<usize>,
    ct: [usize; 2],
}

#[derive(Default)]
struct DelCluster {
    alleles: Vec<DelAllele>,
    mean_x: f64,
    mean_len: f64,
    ss_x: f64,
    ss_len: f64,
    n: [usize; 2],
}

impl fmt::Display for DelCluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.2}\t{:.2}\t{}",
            self.mean_x,
            self.mean_len,
            self.n[0] + self.n[1]
        )
    }
}

fn update_mean_ss(mean: &mut f64, ss: &mut f64, n: f64, mean1: f64, ss1: f64, n1: f64) {
    let delta = mean1 - *mean;
    *mean += n1 * delta / n;
    *ss += ss1 + n1 * delta * (mean1 - *mean);
}

impl DelCluster {
    fn total(&self) -> usize {
        self.n[0] + self.n[1]
    }

    // Add DelAllele to a cluster
    fn add_allele(&mut self, del: DelAllele) -> &mut Self {
        self.n[0] += del.ct[0];
        self.n[1] += del.ct[1];
        let n = (del.ct[0] + del.ct[1]) as f64;
        let nn = self.total() as f64;
        update_mean_ss(&mut self.mean_x, &mut self.ss_x, nn, del.x as f64, 0.0, n);
        update_mean_ss(
            &mut self.mean_len,
            &mut self.ss_len,
            nn,
            del.len as f64,
            0.0,
            n,
        );
        self.alleles.push(del);
        self
    }

    // Merge two clusters
    fn merge(&mut self, other: Self) {
        self.n[0] += other.n[0];
        self.n[1] += other.n[1];
        let n1 = other.total() as f64;
        let n = self.total() as f64;
        update_mean_ss(
            &mut self.mean_x,
            &mut self.ss_x,
            n,
            other.mean_x,
            other.ss_x,
            n1,
        );
        update_mean_ss(
            &mut self.mean_len,
            &mut self.ss_len,
            n,
            other.mean_len,
            other.ss_len,
            n1,
        );
        self.alleles.extend_from_slice(&other.alleles);
    }

    fn dist(&self, other: &Self) -> f64 {
        // Find out how much the SS of the combined cluster would be increased over the sum of the SS of the component clusters
        let nn = (self.total() + other.total()) as f64;
        let mut mn = self.mean_x;
        let mut ss_x = self.ss_x;
        update_mean_ss(
            &mut mn,
            &mut ss_x,
            nn,
            other.mean_x,
            other.ss_x,
            other.total() as f64,
        );
        mn = self.mean_len;
        let mut ss_len = self.ss_len;
        update_mean_ss(
            &mut mn,
            &mut ss_len,
            nn,
            other.mean_len,
            other.ss_len,
            other.total() as f64,
        );
        ((ss_x + ss_len) / nn - (self.ss_x + self.ss_len) / (self.total() as f64)
            + (other.ss_x + other.ss_len) / (other.total() as f64))
            .sqrt()
    }

    fn get_median_and_percentiles<K, F>(&mut self, pcntg: f64, field: F) -> Option<(K, K, K)>
    where
        F: Fn(&DelAllele) -> K,
        K: Ord + Copy,
    {
        assert!((0.0..=1.0).contains(&pcntg));
        let n = self.total();
        if n == 0 {
            None
        } else if n == 1 {
            let x = field(&self.alleles[0]);
            Some((x, x, x))
        } else {
            self.alleles.sort_unstable_by_key(&field);
            // Find entries corresponding to median and percentiles
            let mut wk = [
                (n >> 1, None),
                (((n as f64) * pcntg) as usize, None),
                (((n as f64) * (1.0 - pcntg)) as usize, None),
            ];
            let mut ct = 0;
            for a in self.alleles.iter() {
                ct += a.v.len();
                let mut finished = true;
                for (x, y) in wk.iter_mut() {
                    if y.is_none() {
                        if ct >= *x {
                            *y = Some(field(a))
                        } else {
                            finished = false
                        }
                    }
                }
                if finished {
                    break;
                }
            }
            Some((
                wk[0].1.take().unwrap(),
                wk[1].1.take().unwrap(),
                wk[2].1.take().unwrap(),
            ))
        }
    }
}

pub struct LargeDeletion {
    pub start: usize,
    pub length: usize,
    pub start_ci: [usize; 2],
    pub length_ci: [usize; 2],
    pub reads: Vec<usize>,
    pub counts: [[usize; 2]; 2],
}

impl LargeDeletion {
    fn make(mut cl: DelCluster, total: [usize; 2]) -> Self {
        // We use the median and 2.5%, 97.5% percentiles independently for the start
        // and length (is there a more logical way to do this?)

        let (start, sx1, sx2) = cl
            .get_median_and_percentiles(0.025, |a| a.x)
            .expect("Empty cluster!");
        let (length, lx1, lx2) = cl
            .get_median_and_percentiles(0.025, |a| a.len)
            .expect("Empty cluster!");

        // Copy indices of all reads involved in deletion to reads vec
        let mut reads = Vec::with_capacity(cl.n[0] + cl.n[1]);
        for a in cl.alleles.iter_mut() {
            reads.extend_from_slice(&a.v)
        }

        // Counts
        assert!(cl.n[0] <= total[0] && cl.n[1] <= total[1]);
        let counts = [[total[0] - cl.n[0], total[1] - cl.n[1]], [cl.n[0], cl.n[1]]];

        LargeDeletion {
            start,
            length,
            start_ci: [sx1, sx2],
            length_ci: [lx1, lx2],
            reads,
            counts,
        }
    }

    pub fn end(&self) -> usize {
        self.start + self.length - 1
    }

    pub fn pos_ci(&self, x: usize) -> isize {
        assert!(x < 2);
        let pos = self.start as isize;
        (self.start_ci[x] as isize) - pos
    }

    pub fn len_ci(&self, x: usize) -> isize {
        assert!(x < 2);
        let len = self.length as isize;
        (self.length_ci[x] as isize) - len
    }

    pub fn n(&self) -> usize {
        self.counts[0][0] + self.counts[0][1] + self.counts[1][0] + self.counts[1][1]
    }

    pub fn fq(&self) -> f64 {
        ((self.counts[1][0] + self.counts[1][1]) as f64) / (self.n() as f64)
    }
}

impl fmt::Display for LargeDeletion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} ({} - {}) {} ({} - {})",
            self.start,
            self.start_ci[0],
            self.start_ci[1],
            self.length,
            self.length_ci[0],
            self.length_ci[1]
        )
    }
}

/// Cluster large deletion alleles by start and end points
fn cluster_del_alleles(dv: &mut Vec<DelCluster>) {
    trace!("Clustering deletion alleles - starting with {}", dv.len());
    while dv.len() > 1 {
        let mut min = (1, 0, dv[0].dist(&dv[1]), dv[0].total() + dv[1].total());
        for (ix, cl) in dv[2..].iter().enumerate() {
            for (ix1, cl1) in dv[..ix + 2].iter().enumerate() {
                let d = cl.dist(cl1);
                let n = cl.total() + cl1.total();
                if d < min.2 || (d - min.2 < 1.0e-8 && n > min.3) {
                    min = (ix + 2, ix1, d, n)
                }
            }
        }
        if min.2 > 256.0 {
            break;
        }
        trace!(
            "Adding cluster {} to {}: distance = {}",
            &dv[min.0], &dv[min.1], min.2
        );
        // Note that swap_remove swaps the selected position with the last element in the list
        // It is important that the second cluster is not the last element otherwise it's index
        // could change!  Normally min.0 > min.1 unless the code above is changed, but we will
        // check in case
        assert!(min.0 > min.1);
        let cl = dv.swap_remove(min.0);
        dv[min.1].merge(cl);
    }
}

pub(crate) fn get_large_deletions(
    res: &[VcfRes],
    j: usize,
    k: usize,
    ref_len: usize,
    dalign: &[DetailedAlign],
    cfg: &Config,
) -> Vec<LargeDeletion> {
    let x = res[j].x;
    let y = res[k].x;
    let sz = if y >= x {
        y + 1 - x
    } else {
        y + ref_len + 1 - x
    };

    let mut obs_all = AllDesc::new(sz);
    let mut cts: Vec<[usize; 2]> = vec![[0; 2]];
    
    #[allow(clippy::type_complexity)]
    let mut ahash: HashMap<(usize, usize), (usize, [usize; 2], Vec<usize>)> = HashMap::new();
    
    let hidx = HIDX as u8;
    let mut n_alls = 1;
    for (read_ix, dr) in dalign.iter().enumerate() {
        if dr
            .extract_deletion_allele(x, y, ref_len, &mut obs_all)
            .is_some()
        {
            let start = obs_all.start();
            let str = if obs_all[0] >= hidx { 1 } else { 0 };
            if str == 1 {
                for c in obs_all.iter_mut() {
                    *c %= hidx;
                }
            }
            // Only take alleles that start with a non-deletion for the first base
            if obs_all[0].is_multiple_of(HIDX as u8) {
                // Find start position and length of longest continuous deletion in allele
                let ix = if let Some((mut del_start, length)) = longest_del(&obs_all) {
                    if length >= cfg.small_deletion_limit() {
                        del_start += start;
                        if let Some((i, ct, v)) = ahash.get_mut(&(del_start, length)) {
                            v.push(read_ix);
                            ct[str] += 1;
                            *i
                        } else {
                            let reads: Vec<usize> = vec![read_ix];
                            let mut ct = [0; 2];
                            ct[str] += 1;
                            ahash.insert((del_start, length), (n_alls, ct, reads));
                            n_alls += 1;
                            cts.push([0; 2]);
                            n_alls - 1
                        }
                    } else {
                        0
                    }
                } else {
                    0
                };
                cts[ix][str] += 1;
                if ix == 0 {
                    trace!("Non-deleted allele: {}", dr.read_id);
                }
            }
        }
    }
    let mut dv = Vec::with_capacity(n_alls);
    for ((start, len), (_, ct, v)) in ahash.drain() {
        let mut clust = DelCluster::default();
        let del = DelAllele {
            x: start + x,
            len,
            ct,
            v,
        };
        clust.add_allele(del);
        dv.push(clust);
    }
    // Cluster together alleles with similar start and end points
    cluster_del_alleles(&mut dv);
    // Remove low frequency clusters
    let total = cts
        .iter()
        .fold([0; 2], |s, ct| [s[0] + ct[0], s[1] + ct[1]]);
    let mut dv_old = dv;
    let mut tmp = Vec::new();
    for cl in dv_old.drain(..) {
        let f = (cl.total() as f64) / ((total[0] + total[1]) as f64);
        if f >= cfg.indel_threshold(ThresholdType::Hard) {
            tmp.push(cl);
        } else {
            trace!("Removing cluster {}: low frequency", cl);
        }
    }
    let total = tmp
        .iter()
        .fold(cts[0], |s, cl| [s[0] + cl.n[0], s[1] + cl.n[1]]);
    let mut dv: Vec<_> = tmp
        .drain(..)
        .map(|cl| LargeDeletion::make(cl, total))
        .collect();
    dv.sort_unstable_by(|d1, d2| match d1.start.cmp(&d2.start) {
        Ordering::Equal => d1.length.cmp(&d2.length),
        c => c,
    });

    if log_enabled!(Trace) {
        trace!(
            "Non deletion allele: {} {}",
            cts[0][0] + cts[0][1],
            ((cts[0][0] + cts[0][1]) as f64) / ((total[0] + total[1]) as f64)
        );
        for ld in dv.iter() {
            trace!("Cluster: {} {}", ld, ld.fq());
            for dv in ld.reads.iter().map(|&i| &dalign[i]) {
                trace!("\t{}", dv.read_id);
            }
        }
    }

    dv
}
