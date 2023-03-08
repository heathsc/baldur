use r_htslib::SamHeader;
use std::{collections::HashMap, io::BufRead, rc::Rc};

const BASE_TRANS: [u8; 256] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[derive(Copy, Clone, Debug)]
pub struct RefPos {
    base: u8,
    hpoly: u8,
    dinuc: u8,
}

impl RefPos {
    pub fn new(base: u8) -> Self {
        assert!(base <= 4);
        Self {
            base,
            hpoly: 0,
            dinuc: 0,
        }
    }
    pub fn base(&self) -> u8 {
        self.base
    }
    pub fn hpoly(&self) -> u8 {
        self.hpoly
    }
}

pub struct Contig {
    tid: usize,
    seq: Vec<RefPos>,
    name: Rc<str>,
}

impl Contig {
    pub fn tid(&self) -> usize {
        self.tid
    }

    pub fn size(&self) -> usize {
        self.seq.len()
    }

    pub fn seq(&self) -> &[RefPos] {
        &self.seq
    }

    pub fn name(&self) -> &str {
        &self.name
    }
}

struct FastaReader<R: BufRead> {
    rdr: R,
    buffer: String,
    line: usize,
}

impl<R: BufRead> FastaReader<R> {
    fn new(rdr: R) -> FastaReader<R> {
        let buffer = String::new();
        Self {
            rdr,
            buffer,
            line: 0,
        }
    }

    fn next_record(&mut self, hdr: &SamHeader) -> anyhow::Result<Option<Contig>> {
        if self.buffer.is_empty() && self.rdr.read_line(&mut self.buffer)? == 0 {
            return Ok(None);
        }
        self.line += 1;
        loop {
            let name = self
                .buffer
                .strip_prefix('>')
                .map(|s| Rc::from(s.trim_end()))
                .ok_or_else(|| anyhow!("Expected '>' at start of line {}", self.line))?;
            trace!("Found contig {}", name);
            let ctg = hdr.name2tid(&name as &str);
            trace!("SAM tid = {:?}", ctg);
            let mut seq = Vec::new();
            loop {
                self.buffer.clear();
                if self.rdr.read_line(&mut self.buffer)? == 0 {
                    break;
                }
                self.line += 1;
                if self.buffer.starts_with('>') {
                    break;
                }
                if ctg.is_some() {
                    for c in self
                        .buffer
                        .trim_end()
                        .as_bytes()
                        .iter()
                        .map(|x| BASE_TRANS[*x as usize])
                    {
                        if c > 0 {
                            seq.push(RefPos::new(c - 1))
                        } else {
                            return Err(anyhow!("Illegal base at line {}", self.line));
                        }
                    }
                }
            }
            if let Some(tid) = ctg {
                // Add homopolymer information
                add_homopolymer_info(&mut seq);

                // Add dinucleotide repeat information
                // Starting from 0...
                add_dinuc_repeat_info(&mut seq);
                // ...and starting from 1
                add_dinuc_repeat_info(&mut seq[1..]);

                return Ok(Some(Contig { name, tid, seq }));
            } else if self.buffer.is_empty() {
                return Ok(None);
            }
        }
    }
}

fn add_homopolymer_info(seq: &mut [RefPos]) {
    struct State {
        cbase: u8,
        count: usize,
    }
    let mut state: Option<State> = None;
    let mut hpoly_list = Vec::with_capacity(128);

    for (i, c) in seq.iter().enumerate() {
        state = match state {
            Some(s) => {
                if c.base == s.cbase {
                    Some(State {
                        cbase: s.cbase,
                        count: s.count + 1,
                    })
                } else {
                    if s.count >= 2 {
                        hpoly_list.push((i - s.count, s.count))
                    }
                    Some(State {
                        cbase: c.base,
                        count: 1,
                    })
                }
            }
            None => Some(State {
                cbase: c.base,
                count: 1,
            }),
        }
    }
    for (i, ct) in hpoly_list.drain(..) {
        let ctm = (ct - 1).min(15);
        for (k, c) in seq[i..i + ct].iter_mut().enumerate() {
            c.hpoly = ((k.min(15) << 4) | ctm) as u8;
        }
    }
}

fn add_dinuc_repeat_info(seq: &mut [RefPos]) {
    struct State {
        dinuc: [u8; 2],
        count: usize,
    }
    let mut state: Option<State> = None;
    let mut dinuc_list = Vec::with_capacity(128);
    for (i, c) in seq
        .chunks_exact(2)
        .map(|v| [v[0].base, v[1].base])
        .enumerate()
    {
        state = match state {
            Some(s) => {
                if c == s.dinuc {
                    Some(State {
                        dinuc: s.dinuc,
                        count: s.count + 1,
                    })
                } else {
                    if s.count >= 2 {
                        dinuc_list.push(((i - s.count) << 1, s.count))
                    }
                    if c[0] == c[1] {
                        None
                    } else {
                        Some(State { dinuc: c, count: 1 })
                    }
                }
            }
            None => {
                if c[0] == c[1] {
                    None
                } else {
                    Some(State { dinuc: c, count: 1 })
                }
            }
        }
    }
    for (i, ct) in dinuc_list.drain(..) {
        let ctm = (ct - 1).min(15);
        for (k, c) in seq[i..i + (ct << 1)].chunks_exact_mut(2).enumerate() {
            c[0].dinuc = ((k.min(15) << 4) | ctm) as u8;
        }
    }
}

#[derive(Default)]
pub struct Reference {
    contigs: HashMap<usize, Contig>,
    name2tid: HashMap<Rc<str>, usize>,
}

impl Reference {
    pub fn from_reader<R: BufRead>(rdr: R, hdr: &SamHeader) -> anyhow::Result<Self> {
        let mut contigs = HashMap::new();
        let mut name2tid = HashMap::new();
        let mut fasta_rdr = FastaReader::new(rdr);
        while let Some(ctg) = fasta_rdr.next_record(hdr)? {
            debug!("Read in contig {}", hdr.tid2name(ctg.tid));
            name2tid.insert(ctg.name.clone(), ctg.tid);
            contigs.insert(ctg.tid, ctg);
        }
        Ok(Reference { contigs, name2tid })
    }

    pub fn n_contigs(&self) -> usize {
        self.contigs.len()
    }

    pub fn contig(&self, tid: usize) -> Option<&Contig> {
        self.contigs.get(&tid)
    }

    pub fn name2contig(&self, s: &str) -> Option<&Contig> {
        self.name2tid.get(s).map(|x| self.contigs.get(x).unwrap())
    }

    pub fn contigs(&self) -> &HashMap<usize, Contig> {
        &self.contigs
    }
}
