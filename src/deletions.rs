use std::{
    collections::{hash_map, HashMap},
    fmt::{self, Formatter},
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
}

impl fmt::Display for Deletion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.start,
            self.end,
            if self.size.is_negative() { '-' } else { '+' },
            self.size,
            self.dtype,
        )
    }
}
/*
impl Deletion {
    pub fn size(&self) -> usize {
        self.size.abs() as usize
    }

    pub fn reversed(&self) -> bool {
        self.size.is_negative()
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }
}*/

/// We store deletions at least as large as min_size in a hash table so we can get summaries of
/// large deletions across all reads
pub struct Deletions {
    del_hash: HashMap<Deletion, usize>,
    min_size: isize,
    target_size: usize,
}

impl Deletions {
    pub fn new(target_size: usize, min_size: usize) -> Self {
        Self {
            del_hash: HashMap::new(),
            min_size: min_size as isize,
            target_size,
        }
    }

    pub fn add_del(&mut self, start: usize, end: usize, reversed: bool, dtype: DelType) {
        let start = start % self.target_size;
        let end = end % self.target_size;
        let size = if reversed {
            let mut s = (start as isize + 1) - (end as isize);
            if s < 1 {
                s += self.target_size as isize
            }
            -s
        } else {
            let mut s = (end as isize + 1) - (start as isize);
            if s < 1 {
                s += self.target_size as isize
            }
            s
        };
        assert!(size.abs() > 0);
        if size.abs() >= self.min_size {
            let del = Deletion {
                start,
                end,
                size,
                dtype,
            };
            let e = self.del_hash.entry(del).or_insert(0);
            *e += 1;
        }
    }

    pub fn iter(&self) -> hash_map::Iter<Deletion, usize> {
        self.del_hash.iter()
    }
}
