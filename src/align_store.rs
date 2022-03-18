use crate::{
   cli::{Config, Region},
   reference::RefPos,
   context::*,
   process::ProcWork,
};

pub(crate) struct AlignStore<'a> {
   seq: Vec<u8>,
   qual: Vec<u8>,
   start: usize,  // start on reference
   pos: isize,    // current position on reference
   adjust: isize, // amount to adjust positions
   target_size: usize,
   cfg: &'a Config,
   changed: bool,
}

impl <'a> AlignStore<'a> {
   pub(crate) fn new(reg: &Region, cfg: &'a Config) -> Self {
      let seq: Vec<u8> = vec![b' '; reg.len()];
      let qual: Vec<u8> = vec![0; reg.len()];
      let pos = -1;
      let adjust = cfg.adjust() as isize;
      let target_size = cfg.region().ctg_size();
      Self {
         seq,
         qual,
         start: reg.start(),
         pos,
         adjust,
         target_size,
         cfg,
         changed: false,
      }
   }

   pub(crate) fn set_pos(&mut self, pos: usize) {
      self.pos = (pos as isize) - self.adjust;
      if self.pos < 1 {
         self.pos += self.target_size as isize
      }
   }

   pub(crate) fn advance(&mut self, n: usize) {
      self.pos += n as isize
   }

   pub(crate) fn in_region(&self) -> Option<usize> {
      let pos = self.pos % (self.target_size as isize);
      if pos >= self.start as isize && pos < (self.start + self.seq.len()) as isize {
         Some((pos as usize) - self.start)
      } else { None }
   }

   pub(crate) fn add(&mut self, src: &[(u8, u8, Context5)], ref_seq: &[RefPos], pw: &mut ProcWork) -> [usize; 2] {

      let is_del = |x: u8| x == b'-' || x == b'_' || x == b'^' || x == b'&';

      let mut ct = [0; 2];
      for (i, (base, c, ctxt)) in src.iter().enumerate() {
         if let Some(x) = self.in_region() {
            self.seq[x] = *base;
            self.qual[x] = *c >> 2;
            self.changed = true;
            let rpos = &ref_seq[x];
            let rb = rpos.base();
            let bs = *c & 3;
            if rb < 4 {
               let mm = if rb != bs { 1 } else { 0 };
               if !is_del(*base) {
                  ct[mm] += 1
               }
               if ctxt.context3().is_some() && self.cfg.rs(x).is_none() {
                  let q = (*c >> 2) as usize;
                  let ct = ctxt.context3().unwrap() as usize;
                  if !is_del(*base) {
                     pw.qual_hist[q][mm] += 1;
                     pw.ctxt_hist[ct][q][mm] += 1;
                     // We don't count the last base of a series as a non-deletion as we can't tell if there will be a deletion on the next base
                     if i < src.len() - 1 {
                        pw.qual_hist[q][2] += 1;
                        pw.ctxt_hist[ct][q][2] += 1;
                     }
                  } else {
                     pw.qual_hist[q][3] += 1;
                     pw.ctxt_hist[ct][q][3] += 1;
                  }
               }
            }
         }
         self.pos += 1;
      }
      ct
   }

   pub(crate) fn fill_to(&mut self, fill: u8, qual: u8, x: usize) {
      assert!(self.pos >= 0);
      let x = (x as isize) - self.adjust;
      loop {
         let pos = self.pos % (self.target_size as isize);
         if pos == x {
            break;
         }
         if self.in_region().is_some() {
            let tp = &mut self.seq[(pos as usize) - self.start];
            if *tp != b' ' {
               break;
            }
            *tp = fill;
            self.qual[(pos as usize) - self.start] = qual;
            self.changed = true;
         }
         self.pos += 1;
      }
   }

   pub(crate) fn seq(&self) -> &[u8] {
      &self.seq
   }

   pub(crate) fn qual(&self) -> &[u8] {
      &self.qual
   }

   pub(crate) fn changed(&self) -> bool {
      self.changed
   }
}

