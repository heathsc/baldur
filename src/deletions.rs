use std::{
    collections::{HashMap, hash_map},
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
    reversed: bool,
    dtype: DelType,
}

impl fmt::Display for Deletion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.start,
            self.end,
            if self.reversed { '-' } else { '+' },
            self.size,
            self.dtype,
        )
    }
}

impl Deletion {
    pub fn start_end_pos(&self) -> (usize, usize) {
        // Set coordinates to forward strand and add 10bp padding
        if self.reversed {
            (self.end, self.start)
        } else {
            (self.start, self.end)
        }
    }
}
/// We store deletions at least as large as min_size in a hash table so we can get summaries of
/// large deletions across all reads
pub struct Deletions {
    del_hash: HashMap<Deletion, usize>,
    read_extents: Vec<[usize; 2]>,
    min_size: isize,
    target_size: usize,
}

impl Deletions {
    pub fn new(target_size: usize, min_size: usize) -> Self {
        Self {
            del_hash: HashMap::new(),
            read_extents: Vec::new(),
            min_size: min_size as isize,
            target_size,
        }
    }

    pub fn add_del(&mut self, start: usize, end: usize, reversed: bool, dtype: DelType) {
        let size = end as isize - start as isize;
        let start = start % self.target_size;
        let end = end % self.target_size;
        if size.abs() >= self.min_size {
            let del = Deletion {
                start,
                end,
                size,
                dtype,
                reversed,
            };
            let ct = self.del_hash.entry(del).or_insert(0);
            *ct += 1;
        }
    }

    pub fn add_read_extent(&mut self, start: usize, end: usize) {
        assert!(end > start);
        self.read_extents.push([start, end]);
    }

    pub fn iter<'a>(&'a self) -> hash_map::Iter<'a, Deletion, usize> {
        self.del_hash.iter()
    }

    pub fn read_extents(&self) -> &[[usize; 2]] {
        &self.read_extents
    }

    pub fn target_size(&self) -> usize {
        self.target_size
    }

    pub fn len(&self) -> usize {
        self.del_hash.len()
    }
}
