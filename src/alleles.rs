use std::{
   collections::HashMap,
   fmt,
   rc::Rc,
   ops::{Deref, DerefMut},
};

use crate::{
   model::N_QUAL,
   depth::*,
   detailed_align::DetailedAlign,
};

// Inserted allele counts
pub(crate) struct InsAllele {
   pub(crate) allele: Rc<[u8]>,
   pub(crate) cts: [usize; 2],
   pub(crate) qcts: [usize; N_QUAL],
}

impl InsAllele {
   pub(crate) fn new(allele: Rc<[u8]>) -> Self {
      Self{allele, cts: [0; 2], qcts: [0; N_QUAL] }
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

#[derive(Hash, PartialEq, Eq)]
pub struct AllDesc(pub Vec<u8>);

impl Deref for AllDesc {
   type Target = [u8];

   fn deref(&self) -> &Self::Target {
      &self.0
   }
}

impl DerefMut for AllDesc {
   fn deref_mut(&mut self) -> &mut Self::Target {
      &mut self.0
   }
}

impl fmt::Display for AllDesc {
   fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      for &i in self.0.iter() {
         write!(f, "{}", BASES.as_bytes()[i as usize] as char)?;
      }
      Ok(())
   }
}

pub(crate) fn get_allele_counts(adesc: &[AllDesc], x1: usize, x: usize, ref_len: usize,
                     dalign: &[DetailedAlign], ins_hash: &HashMap<usize, InsAlleles>)
                     -> (Vec<[usize; 2]>, Vec<[usize; N_QUAL]>, Vec<bool>) {

   let ref_all_len = adesc[0].len();
   let indel_flag: Vec<_> = adesc.iter().map(|all| all.len() != ref_all_len).collect();

   let mut all_hash = HashMap::with_capacity(2 * adesc.len());
   for (ix, ds) in adesc.iter().enumerate() {
      let ds1 = AllDesc(ds.to_vec());
      // Forward allele
      all_hash.insert(ds1, (ix, 0));
      let ds2 = AllDesc(ds.iter().map(|&c| c + HIDX as u8).collect());
      // Reverse allele
      all_hash.insert(ds2, (ix, 1));
   }

   let mut obs_all = AllDesc(Vec::with_capacity(x + 1 - x1));
   let n_alls = adesc.len();
   let mut cts: Vec<[usize; 2]> = vec!([0; 2]; n_alls);
   let mut qcts: Vec<[usize; N_QUAL]> = vec!([0; N_QUAL]; n_alls);
   for dr in dalign.iter() {
      let q = dr.extract_allele(x1, x, ref_len, &mut obs_all, ins_hash);
      if obs_all.is_empty() {
         obs_all.0.push(if dr.reverse() { 4 + HIDX as u8 } else { 4 });
      }
      if q < N_QUAL as u8 {
         if let Some((ix, str)) = all_hash.get(&obs_all) {
            cts[*ix][*str] += 1;
            qcts[*ix][q as usize] += 1;
         }
      }
   }
   (cts, qcts, indel_flag)
}

