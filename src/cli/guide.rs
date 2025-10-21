use std::{collections::HashMap, io::BufRead, rc::Rc};

use anyhow::Context;
use rstar::AABB;

use crate::reference::Reference;

#[derive(Debug)]
pub struct Guide {
    name: Rc<str>,
    tid: usize,
    cut_site: isize,
    reverse_strand: bool,
}

impl Guide {
    pub fn tid(&self) -> usize {
        self.tid
    }

    pub fn cut_site(&self) -> isize {
        self.cut_site
    }
    
    pub fn get_bb(&self, rev: bool, ctg_size: isize) -> (AABB<[isize; 2]>, AABB<[isize; 2]>) {
        let cs = self.cut_site;

        let (x1, y1) = match (self.reverse_strand, rev) {
            (false, false) => (cs - 17, cs + 50),
            (false, true) => (cs - 50, cs + 5),
            (true, false) => (cs - 5, cs + 50),
            (true, true) => (cs - 50, cs + 17),
        };

        let (x2, y2) = if y1 < ctg_size {
            (x1 + ctg_size, y1 + ctg_size)
        } else {
            (x1 - ctg_size, y1 - ctg_size)
        };

        (
            AABB::from_corners([x1, 0], [y1, 1]),
            AABB::from_corners([x2, 0], [y2, 1]),
        )
    }

    pub fn from_str(s: &str, rf: &Reference) -> anyhow::Result<Self> {
        let fd: Vec<_> = s.trim().split('\t').take(4).collect();
        if fd.len() != 4 {
            return Err(anyhow!("Illegal short line: {s} (exepcted 4 columns)"));
        }
        let tid = rf
            .name2contig(fd[0])
            .map(|c| c.tid())
            .ok_or(anyhow!("Unknown contig {}", fd[0]))?;
        let cut_site = fd[1]
            .parse::<isize>()
            .with_context(|| format!("Illegal cut site: {}", fd[1]))?;
        let reverse_strand = match fd[2] {
            "+" | "" => false,
            "-" => true,
            _ => return Err(anyhow!("Illegal strand: {}", fd[2])),
        };
        let name: Rc<str> = Rc::from(fd[3]);
        Ok(Self {
            name,
            tid,
            cut_site,
            reverse_strand,
        })
    }
}

#[derive(Default)]
pub struct Guides {
    hash: HashMap<Rc<str>, Guide>,
}

impl Guides {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, guide: Guide) {
        let s = guide.name.clone();
        let s1 = s.clone();
        if self.hash.insert(s, guide).is_some() {
            warn!("Duplicate entry for guide {s1}")
        }
    }

    pub fn get(&self, name: &str) -> Option<&Guide> {
        self.hash.get(name)
    }

    pub fn from_reader<R: BufRead>(mut rdr: R, rf: &Reference) -> anyhow::Result<Self> {
        let mut guides = Self::new();
        let mut s = String::new();
        while rdr
            .read_line(&mut s)
            .with_context(|| "Error reading from guides file")?
            != 0
        {
            let g = Guide::from_str(&s, rf)?;
            guides.add(g);
            s.clear();
        }

        debug!("{} guide definitions read in", guides.hash.len());
        Ok(guides)
    }
}
