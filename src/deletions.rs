use std::{
    collections::{HashMap, hash_map},
    fmt::{self, Formatter},
    fs::File,
    io::{BufWriter, Write},
};

#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
pub enum DelType {
    Del,
    Split,
}

impl fmt::Display for DelType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Del => "Deletion",
                Self::Split => "Split",
            }
        )
    }
}

/// A deletion found within a read.  The sign of `size` indicates the
/// the strand.  Note that because we are dealing with circular genomes
/// we can not infer the strand from the start and end positions alone
#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
pub struct Deletion {
    start: usize,
    end: usize,
    size: isize,
    dtype: DelType,
    ix: usize,
    ct: [usize; 2],
}

impl fmt::Display for Deletion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.start, self.end, self.size, self.dtype, self.ix, self.ct[0], self.ct[1]
        )
    }
}

impl Deletion {
    pub fn start_end_pos(&self) -> (usize, usize) {
        (self.start, self.end)
    }
}

/// We store deletions at least as large as min_size in a hash table so we can get summaries of
/// large deletions across all reads
pub struct Deletions {
    del_hash: HashMap<(usize, usize), Deletion>,
    read_extents: Vec<[isize; 2]>,
    min_size: isize,
    target_size: usize,
    adjust: isize,
}

impl Deletions {
    pub fn new(target_size: usize, min_size: usize, adjust: usize) -> Self {
        Self {
            del_hash: HashMap::new(),
            read_extents: Vec::new(),
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
    ) -> Option<usize> {
        let s = start as isize - self.adjust;
        let e = end as isize - self.adjust;

        let size = e - s + if s <= e { 1 } else { -1 };

        let adj = |x: isize| (x % (self.target_size as isize)) as usize;

        if size.abs() >= self.min_size {
            let (start, end) = if reversed {
                (adj(e), adj(s))
            } else {
                (adj(s), adj(e))
            };

            let l = self.del_hash.len();
            let del = self
                .del_hash
                .entry((start, end))
                .or_insert_with(|| Deletion {
                    start,
                    end,
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

    pub fn add_read_extent(&mut self, start: usize, end: usize) {
        assert!(end > start);
        self.read_extents
            .push([start as isize - self.adjust, end as isize - self.adjust]);
    }

    pub fn iter<'a>(&'a self) -> hash_map::Values<'a, (usize, usize), Deletion> {
        self.del_hash.values()
    }

    pub fn read_extents(&self) -> &[[isize; 2]] {
        &self.read_extents
    }

    pub fn target_size(&self) -> usize {
        self.target_size
    }

    pub fn len(&self) -> usize {
        self.del_hash.len()
    }

    pub fn write_deletions(&self, prefix: &str) -> anyhow::Result<()> {
        let file_name = format!("{}_del.txt", prefix);
        let mut wrt = BufWriter::new(File::create(&file_name)?);

        for (_, d) in self.del_hash.iter() {
            writeln!(wrt, "{d}")?;
        }

        Ok(())
    }
}

pub struct DeletionWork {
    dels: Deletions,
    mask: Box<[u64]>,
    reads: Box<[[isize; 2]]>,
    del_cts: Box<[(usize, usize)]>,
    n_u64: usize,
}

const EXTEND_DEL: usize = 30;

impl DeletionWork {
    pub fn new(dels: Deletions) -> Self {
        let re = dels.read_extents();
        let n_dels = dels.len();
        let n_u64 = (n_dels + 63) >> 6;
        let mut mask = vec![0u64; re.len() * n_u64];
        let target_size = dels.target_size();

        for d in dels.iter() {
            let (s, e) = d.start_end_pos();
            let overlap_origin = e < s;
            let start = s.saturating_sub(EXTEND_DEL);
            let end = if overlap_origin {
                e + target_size + EXTEND_DEL
            } else {
                e + EXTEND_DEL
            };
            
            let i = d.ix;
            
            let (j, mk) = (i >> 6, 1u64 << (i & 63));

            for (x, m) in re.iter().zip(mask.chunks_mut(n_u64)) {
                // Is del covered by read?
                let tst = |a, b| (a as isize) >= x[0] && (b as isize) <= x[1];

                if tst(start, end)
                    || (!overlap_origin && tst(start + target_size, end + target_size))
                {
                    m[j] |= mk
                }
            }
        }
        let mut mask1 = Vec::new();
        let mut reads = Vec::new();
        if n_u64 > 0 {
            for (r, m) in re.iter().zip(mask.chunks(n_u64)) {
                if m.iter().any(|b| *b != 0) {
                    mask1.extend_from_slice(m);
                    reads.push(*r)
                }
            }
        }
        let del_cts: Box<[(usize, usize)]> = vec![(0,0); n_dels].into_boxed_slice();
        Self {
            dels,
            mask: mask1.into_boxed_slice(),
            reads: reads.into_boxed_slice(),
            del_cts,
            n_u64,
        }
    }

    pub fn get_del_prob(&mut self, x: usize) -> f64 {
        let n_u64 = self.n_u64;

        let mut msk = vec![0u64; n_u64];
        let mut msk1 = Vec::with_capacity(n_u64);

        let target_size = self.dels.target_size();

        let del_cts = &mut self.del_cts;
        let re = &self.reads;
        let mask = &self.mask;

        del_cts.fill((0,0));
        for d in self.dels.iter() {
            let (start, end) = d.start_end_pos();
            let z = if start <= end {
                (start..=end).contains(&x)
            } else {
                x >= start || x <= end
            };

            let i = d.ix;
            if z {
                let (j, mk) = (i >> 6, 1u64 << (i & 63));
                msk[j] |= mk;
                del_cts[i] = (d.ct[0] + d.ct[1], 0)
            }
        }
        if n_u64 > 0 {
            for (r, m) in re.iter().zip(mask.chunks(n_u64)) {
                msk1.clear();
                for (m1, m2) in m.iter().zip(msk.iter()) {
                    msk1.push(m1 & m2)
                }
                if msk1.iter().any(|y| *y != 0)
                    && ((r[0]..r[1]).contains(&(x as isize))
                        || (r[0]..r[1]).contains(&((x + target_size) as isize)))
                {
                    for (j, m) in msk1.iter().enumerate() {
                        let mut x = *m;
                        let mut i = 0;
                        while x != 0 {
                            if (x & 1) == 1 {
                                del_cts[j * 64 + i].1 += 1;
                            }
                            i += 1;
                            x >>= 1
                        }
                    }
                }
            }
        }

        let mut z = 1.0;
        for d in self.dels.iter() {
            let (a, b) = del_cts[d.ix];
            if a > b {
                panic!("Odd missed read for del at x: {x} - {a} {b}  {d} {:?}", d.ct);
            } else if b > 0 {
                let p = a as f64 / b as f64;
                z *= 1.0 - p;
            }
        }
        1.0 - z
    }
}
