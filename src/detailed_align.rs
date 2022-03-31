use std::collections::HashMap;

use crate::{
   model::N_QUAL,
   depth::*,
   alleles::{AllDesc, InsAlleles, Trunc},
};

pub(crate) struct DetailedAlign {
   pub(crate) seq_qual: Box<[(u8, u8)]>,
   pub(crate) ins_hash: HashMap<usize, (usize, u8)>,
   pub(crate) start: usize,
   pub(crate) read_id: String,
}

impl DetailedAlign {
   pub(crate) fn reverse(&self) -> bool {
      self.seq_qual[0].0 >= HIDX as u8
   }

   pub(crate) fn extract_allele(&self, x: usize, y: usize, ref_len: usize, v: &mut AllDesc, global_ins_hash: &HashMap<usize, InsAlleles>) -> u8 {
      assert!(x < ref_len && y < ref_len);

      let mut x1 = x;
      let x = if x < self.start { x + ref_len } else { x };
      let y = if y < x { y + ref_len } else { y };
      assert!(x <= y);
      let end_x = self.start + self.seq_qual.len();
      v.clear();
      let mut q = N_QUAL as u8;
      // Only get complete alleles
      if x >= self.start && y < end_x {
         for (c, qual) in self.seq_qual[x - self.start ..= y - self.start].iter( ) {
            let c1 = *c % (HIDX as u8);
            if c1 == BASES_INS as u8 {
               // Handle insertions
               let strand_corr = if self.reverse() { HIDX as u8 } else { 0 };
               // Get known insertion alleles for this position
               let ins_alls = global_ins_hash.get(&x1).expect("Missing insertion!");
               // Get allele number for this read
               let (ix, q1) = self.ins_hash.get(&x1).expect("Missing insertion index for read");
               // Get allele
               let all = &ins_alls.alleles[*ix];
               for c1 in all.allele.iter() {
                  v.push(*c1 + strand_corr)
               }
               q = q.min(*q1);
            } else {
               v.push(if c1 == BASES_SKIP as u8 { *c - (BASES_SKIP as u8 - 4) } else { *c });
               q = q.min(*qual);
            }
            x1 += 1;
         }
      }
      q
   }

   pub(crate) fn process_allele<F>(&self, x: usize, y: usize, mut f: F)
   where
      F: FnMut(usize, u8, u8),
   {
      assert!(x <= y);

      let mut x1 = x;
      let end_x = self.start + self.seq_qual.len();
      // Only get complete alleles
      if x >= self.start && y < end_x {
         for (c, qual) in self.seq_qual[x - self.start ..= y - self.start].iter( ) {
            f(x1, *c, *qual);
            x1 += 1;
         }
      }
   }
/*
   pub(crate) fn extract_filtered_allele<'a>(&self, x: usize, y: usize, ref_len: usize, v: &mut AllDesc, global_ins_hash: &HashMap<usize, InsAlleles>,
                                             vres: impl IntoIterator<Item = &'a VcfRes>) -> Option<u8> {
      assert!(x < ref_len && y < ref_len);

      let mut x1 = x;
      let x = if x < self.start { x + ref_len } else { x };
      let y = if y < x { y + ref_len } else { y };
      assert!(x <= y);
      let end_x = self.start + self.seq_qual.len();
      v.clear();
      let mut qsum: usize = 0;
      let mut qn: usize = 0;
      // Only get complete alleles
      if x >= self.start && y < end_x {
         let strand_corr = if self.reverse() { HIDX as u8 } else { 0 };
         for ((c, qual), vr) in self.seq_qual[x - self.start ..= y - self.start].iter( ).zip(vres) {
            let c1 = *c % (HIDX as u8);
            if c1 == BASES_INS as u8 {
               // Handle insertions
               // Get known insertion alleles for this position
               let ins_alls = global_ins_hash.get(&x1).expect("Missing insertion!");
               // Get allele number for this read
               let (ix, q1) = self.ins_hash.get(&x1).expect("Missing insertion index for read");
               let allele_good = if let Some(adv) = vr.adesc.as_ref() {
                  adv.iter().any(|ad| if let Some(k) = ad.ix() { k == *ix } else { false })
               } else { false };
               if allele_good {
                  // Get allele
                  let all = &ins_alls.alleles[*ix];
                  for c1 in all.allele.iter() {
                     v.push(*c1 + strand_corr)
                  }
               } else {
                  v.push(5 + strand_corr)
               }
               qsum += *q1 as usize;
               qn += 1;
            } else {
               let (c1, c) = if c1 == BASES_SKIP as u8{
                  (4, *c - (BASES_SKIP as u8 - 4))
               } else {
                  (c1, *c)
               };
               let c = if let Some(ad) = vr.adesc.as_ref() {
                  if vr.alleles.iter().any(|a| ad[a.ix][0] == c1) {
                     c
                  } else { 5 + strand_corr }
               } else if vr.alleles.iter().any(|a| a.ix == c1.into()) {
                  c
               } else {
                  5 + strand_corr
               };
               v.push(c);
               qsum += *qual as usize;
               qn += 1;
            }
            x1 += 1;
         }
      }
      if qn > 0 { Some( (qsum / qn) as u8) } else { None }
   } */

   pub(crate) fn extract_deletion_allele<'a>(&self, x: usize, y: usize, ref_len: usize, v: &mut AllDesc) -> Option<u8> {
      assert!(x < ref_len && y < ref_len);

      let x = if x < self.start { x + ref_len } else { x };
      let y = if y < x { y + ref_len } else { y };
      assert!(x <= y);
      let end_x = self.start + self.seq_qual.len();
      v.clear();
      let mut qsum: usize = 0;
      let mut qn: usize = 0;
      // Check if we overlap requested region
      if end_x > x && self.start <= y {
         // Handle truncated alleles
         let (a, b) = match (x >= self.start, y < end_x) {
            (true, true) => {
               // Allele is not truncated
               v.set_trunc(Trunc::No);
               (x - self.start, y + 1 - self.start)
            },
            (true, false) => {
               // Right side is truncated
               v.set_trunc(Trunc::Right);
               (x - self.start, self.seq_qual.len())
            },
            (false, true) => {
               // Left side is truncated
               v.set_trunc(Trunc::Left(self.start - x));
               (0, y + 1 - self.start)
            },
            (false, false) => {
               // Both sides are truncated
               v.set_trunc(Trunc::Both(self.start - x));
               (0, self.seq_qual.len())
            },
         };
         let strand_corr = if self.reverse() { HIDX as u8 } else { 0 };
         for (c, qual) in self.seq_qual[a .. b].iter() {
            let c1 = *c % (HIDX as u8);
            let c = if c1 == 4 || c1 == BASES_SKIP as u8 || c1 == (4 + BASES_LOWQ as u8) {
               4 + strand_corr
            } else {
               strand_corr
            };
            v.push(c);
            qsum += *qual as usize;
            qn += 1;
         }
      }
      if qn > 0 { Some( (qsum / qn) as u8) } else { None }
   }

   fn from_vec_slice(seq: &[u8], qual: &[u8], ins_hash: HashMap<usize, (usize, u8)>, start: usize, hash: &HashMap<u8, usize>, id: &str) -> Self {
      assert_eq!(seq.len(), qual.len());
      let mut seq_qual = Vec::with_capacity(seq.len());
      for(s, &q) in seq.iter().zip(qual.iter())
         .map(|(o, x)| (hash.get(o).map(|y| *y as u8)
                           .unwrap_or(BASES_LEN as u8), x))
      {
         seq_qual.push((s, q))
      }
      Self{ seq_qual: seq_qual.into_boxed_slice(), ins_hash, start, read_id: id.to_string()}
   }

   // Find block of observed read
   pub(crate) fn from_obs_vec(seq: &[u8], qual: &[u8], ins_hash: HashMap<usize, (usize, u8)>, hash: &HashMap<u8, usize>, id: &str) -> Option<Self> {
      if seq[0] == b' ' {
         // First base if unobserved.  Find first observed base
         if let Some((x, _)) = seq.iter().enumerate()
            .find(|(_, &c)| c != b' ') {
            // Find first base from end that is observed.  We can safely call unwrap()
            // because we know at least 1 base is observed
            let y = seq.iter().enumerate().rev()
               .find(|(_, &c)| c != b' ').map(|(i, _)| i)
               .unwrap();
            // Clone and make boxed slice of selected section
            Some(Self::from_vec_slice(&seq[x..=y],&qual[x..=y], ins_hash, x, hash, id))
         } else {
            // No observed base, so return None
            None
         }
      } else if let Some((y, _)) = seq.iter().enumerate()
         .find(|(_, &c)| c == b' ') {
         // First base is observed, and y is index of first unobserved base
         // Find next observed base (if this exists)
         if let Some(x) = seq[y + 1..].iter().enumerate()
            .find(|(_, &c)| c != b' ')
            .map(|(i, _)| i + y + 1) {
            // Copy observed section
            let l = seq[x..].len() + y;
            let mut seq_qual = Vec::with_capacity(l);
            for(s, &q) in seq[x..].iter().zip(qual[x..].iter())
               .map(|(o, x)| (hash.get(o).map(|y| *y as u8)
                                 .unwrap_or(BASES_LEN as u8), x)) {
               seq_qual.push((s, q))
            }
            for(s, &q) in seq[..y].iter().zip(qual[..y].iter())
               .map(|(o, x)| (hash.get(o).map(|y| *y as u8)
                                 .unwrap_or(BASES_LEN as u8), x)) {
               seq_qual.push((s, q))
            }
            Some(Self{seq_qual: seq_qual.into_boxed_slice(), ins_hash, start: x, read_id: id.to_string()})
         } else {
            // No following observed base, so just copy this first block
            Some(Self::from_vec_slice(&seq[..y],&qual[..y], ins_hash, 0, hash, id))
         }
      } else {
         // All bases are observed
         Some(Self::from_vec_slice(seq, qual, ins_hash, 0, hash, id))
      }
   }
}
