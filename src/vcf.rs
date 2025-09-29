use std::{
   fmt::{self, Write as fmtWrite},
   io::{self, BufWriter, Write},
};

use log::Level::Trace;

use compress_io::{
   compress::{CompressIo, Writer},
   compress_type::CompressType,
};

use r_htslib::*;

use crate::{
   model::{freq_mle, AlleleRes, N_QUAL, MAX_PHRED, ModelRes},
   mann_whitney::mann_whitney,
   cli::{Config, ThresholdType},
   fisher::FisherTest,
   reference::RefPos,
   alleles::{AllDesc, LargeDeletion},
   stat_funcs::pnormc,
};

const FLT_FS: u32 = 1;
const FLT_QUAL_BIAS: u32 = 2;
const FLT_BLACKLIST: u32 = 4;
const FLT_Q30: u32 = 8;
const FLT_LOW_FREQ: u32 = 16;
const FLT_HOMO_POLY: u32 = 32;

const FLT_STR: [(&str, &str); 6] = [
   ("strand_bias", "Alleles unevenely distributed across strands", ),
   ("qual_bias", "Minor allele has lower quality values"),
   ("blacklist", "Position on black list"),
   ("q30", "Quality < 30"),
   ("low_freq", "Low heteroplasmy frequency"),
   ("homopolymer", "Indel starting in homopolymer region"),
];

pub(crate) struct VcfCalc<'a, 'b, 'c> {
   pub(crate) ftest: &'a FisherTest,
   pub(crate) seq_len: usize,
   pub(crate) ref_seq: &'b [RefPos],
   pub(crate) homopolymer_limit: u8,
   pub(crate) all_desc: Vec<AllDesc>, // AllDesc for 'standard' 5 alleles (A, C, G, T, Del)
   pub(crate) cfg: &'c Config,
}

pub(crate) struct Allele {
   pub(crate) res: AlleleRes,
   pub(crate) ix: usize,
}

pub(crate) struct VcfRes {
   pub(crate) alleles: Vec<Allele>,
   pub(crate) adesc: Option<Vec<AllDesc>>,
   pub(crate) x: usize,
   pub(crate) phred: u8,
}

// Additional allele specific results
#[derive(Default, Copy, Clone)]
struct ExtraRes {
   flt: u32,
   avg_qual: f64,
   fisher_strand: f64,
   wilcox: f64,
}

impl<'a, 'b, 'c> VcfCalc<'a, 'b, 'c> {
   pub fn get_allele_freqs(&self, x: usize, cts: &[[usize; 2]], qcts: &[[usize; N_QUAL]]) -> VcfRes {
      // Should only be used for single base variants where we expect 5 'alleles'
      // for A, C, G, T and Del
      assert_eq!(cts.len(), 5);
      let indel_flags = [false, false, false, false, true];
      let ref_ix = (self.ref_seq[x].base()) as usize;

      self.est_allele_freqs(x, ref_ix, cts, qcts, &indel_flags)
   }

   pub fn get_mallele_freqs(&self, x: usize, cts: &[[usize; 2]], qcts: &[[usize; N_QUAL]], indel_flags: &[bool]) -> VcfRes {
      self.est_allele_freqs(x,0, cts, qcts, indel_flags)
   }

   pub fn est_allele_freqs(&self, x: usize, ref_ix: usize, cts: &[[usize; 2]], qcts: &[[usize; N_QUAL]], indel_flags: &[bool]) -> VcfRes {

      let n_alls = cts.len();
      assert_eq!(n_alls, qcts.len());

      let qual_model = self.cfg.qual_table();

      // Fold across strands
      let jcts: Vec<usize> = cts.iter().map(|x| x[0] + x[1]).collect();

      // Sort possible alleles by reverse numeric order on counts except that the reference base is always first
      let alleles: Vec<usize> = {
         let mut ix: Vec<usize> = (0..n_alls).filter(|x| *x != ref_ix).collect();
         ix.sort_unstable_by(|a, b| jcts[*b].cmp(&jcts[*a]));

         // Remove alleles where alleles not seen on both strands (apart from reference)
         let mut ix1 = Vec::with_capacity(n_alls);
         ix1.push(ref_ix); // Reference allele;
         for &k in ix.iter() {
            if cts[k][0] > 0 && cts[k][1] > 0 && cts[k][0] + cts[k][1] > 2 {
               ix1.push(k)
            }
         }
         ix1
      };

      let mut mr = freq_mle(&alleles, qcts, qual_model);
      if log_enabled!(Trace) {
         trace!("mle freq. estimates");
         for &k in alleles.iter() {
            trace!("{}\t{}\t{}", k, mr.alleles[k].freq, indel_flags[k]);
         }
      }
      // Remove non-reference alleles where the frequencies were estimated below the thresholds
      let snv_lim = self.cfg.snv_threshold(ThresholdType::Hard);
      let indel_lim = self.cfg.indel_threshold(ThresholdType::Hard);

      for ar in mr.alleles.iter_mut() { ar.flag = false }

      let alleles: Vec<_> = alleles.iter().enumerate()
         .filter(|(i, k)| *i == 0 || mr.alleles[**k].freq >= match indel_flags[**k] {
            true => indel_lim,
            false => snv_lim,
         })
         .map(|(_, &k)| k).collect();

      for &i in alleles.iter() {
         mr.alleles[i].flag =true;
      }

      // Adjust the frequencies to account for any alleles that have been filtered
      let tot = alleles.iter().fold(0.0, |s, &x| s + mr.alleles[x].freq);

      // Rescale if necessary
      if tot < 1.0 {
         assert!(tot > 0.0);
         for ar in mr.alleles.iter_mut() {
            if ar.flag {
               ar.freq /= tot
            } else {
               ar.freq = 0.0;
            }
         }
      }

      let ModelRes{alleles: all, phred} = mr;
      let (all, phred) = (all, phred);

      let alleles: Vec<_> = alleles.iter().filter(|&&k| all[k].flag)
         .map(|&k| Allele{ix: k, res: all[k]}).collect();

      VcfRes{alleles, adesc: None, x, phred}
   }

   // Generate VCF output line for large deletions
   pub fn del_output(&self, del: &LargeDeletion) -> String {
      let mut f = String::new();

      let cts = del.counts;

      let fq = del.fq();
      let sd = (fq * (1.0 - fq) / (del.n() as f64)).sqrt();
      let z = fq / sd;
      let phred = if z > 10.0 { MAX_PHRED }
      else {
         (pnormc(z).log10()*-10.0).round().min(MAX_PHRED as f64) as u8
      };

      let flt = if phred >= 30 { 0 } else { FLT_Q30 };

      // ALT, QUAL, FILTER
      let _ = write!(f, "<DEL>\t{}\t{}", phred, Filter(flt));
      // INFO
      let _ = write!(f, "\tSVTYPE=DEL;END={};SVLEN={};CIPOS={},{};CILEN={},{}", del.end() + 1,
                     del.length, del.pos_ci(0), del.pos_ci(1), del.len_ci(0), del.len_ci(1));
      // FORMAT
      let _ = write!(f, "\tGT:ADF:ADR:HPL\t0/1:{},{}:{},{}:{:.5}",
                     cts[0][0], cts[1][0], cts[0][1], cts[1][1], fq);
      f
   }

   // Generate Optional String with VCF output line
   pub fn output(&self, vr: &mut VcfRes, cts: &[[usize; 2]], qcts: &[[usize; N_QUAL]]) -> Option<String> {

      let x = vr.x;

      let raw_depth = cts.iter().fold(0, |t, x| t + x[0] + x[1]);

      let thresh = 0.05 / (self.seq_len as f64);

      // Sort alleles by frequency (keeping reference alleles at position 0)

      vr.alleles[1..].sort_unstable_by(|a1, a2| a2.res.freq.partial_cmp(&a1.res.freq).unwrap());

      // Find index of major allele
      let (major_idx, mj_idx) = vr.alleles.iter().enumerate().max_by(|(_, ar1), (_, ar2)| ar1.res.freq.partial_cmp(&ar2.res.freq).unwrap())
         .map(|(i, ar)| (ar.ix, i)).unwrap();

      // Reference allele
      let ref_ix = vr.alleles[0].ix;

      // Filter cutoffs

      let snv_soft_lim =  self.cfg.snv_threshold(ThresholdType::Soft);
      let indel_soft_lim =  self.cfg.indel_threshold(ThresholdType::Soft);

      let desc = vr.adesc.as_ref().unwrap_or(&self.all_desc);

      // Extra per allele results
      let mut all_res: Vec<_> = vr.alleles.iter().map(|ar| {
         // Average quality
         let (n, s) = qcts[ar.ix].iter().enumerate().fold((0, 0), |(n, s), (q, &ct)| {
            (n + ct, s + ct * q)
         });
         let avg_qual = if n > 0 { s as f64 / n as f64 } else { 0.0 };

         let mut flt = 0;

         // Fisher strand test
         let fisher_strand = if ar.ix != major_idx {
            let k = ar.ix;
            let k1 = major_idx;
            let all_cts =  [cts[k][0], cts[k1][0], cts[k][1], cts[k1][1]];
            let fs = self.ftest.fisher(&all_cts);
            if ar.res.freq < 0.75 && fs <= thresh {
               flt |= FLT_FS
            }
            fs
         } else { 1.0 };

         // Wilcoxon-Mann-Whitney test for quality bias between the major (most frequent)
         // allele and all minor alleles
         let wilcox = if ar.ix != major_idx {
            mann_whitney(qcts, major_idx, ar.ix).unwrap_or(1.0)
         } else { 1.0 };

         // Set allele freq. flags
         let (lim, pos_adjust) = if desc[ar.ix].len() != desc[ref_ix].len() {
            // This is an indel (size difference from reference)
            let hp_size = (self.ref_seq[x].hpoly() & 0xf).max(self.ref_seq[x + 1].hpoly() & 0xf) + 1;
            if hp_size >= self.homopolymer_limit {
               flt |= FLT_HOMO_POLY
            }
            (indel_soft_lim, 0)
         } else {
            // In a complex variant, a SNV could start a few bases after the location of the variant
            let x = if ar.ix == ref_ix { 0 } else {
               // Find position of first base that differs between this allele and the reference
               desc[ref_ix].iter().zip(desc[ar.ix].iter()).enumerate()
                  .find(|(_, (c1, c2))| *c1 != *c2).map(|(ix, _)| ix as u32).unwrap()
            };
            (snv_soft_lim, x)
         };
         if ar.res.freq < lim {
            flt |= FLT_LOW_FREQ
         }

         // LR flags
         if ar.res.lr_test < 30 {
            flt |= FLT_Q30
         }

         // Blacklist
         if self.cfg.blacklist(x + 1 + pos_adjust as usize) {
            flt |= FLT_BLACKLIST
         }

         ExtraRes{ flt, avg_qual, fisher_strand, wilcox}
      }).collect();

      // Set wilcox allele flags
      let mjr_qual = all_res[mj_idx].avg_qual;

      for res in all_res.iter_mut() {
         if res.wilcox <= thresh && mjr_qual - res.avg_qual > 2.0 {
            res.flt |= FLT_QUAL_BIAS
         }
      }

      // Genotype call
      let f0 = vr.alleles[0].res.freq >= 1.0e-5;
      let gt = match vr.alleles.len() {
         1 => String::from("0"),
         n => if f0 {
            let mut s = String::from("0");
            for i in 1..n {
               s = format!("{}/{}", s, i);
            }
            s
         } else {
            let mut s = String::from("1");
            for i in 2..n {
               s = format!("{}/{}", s, i);
            }
            s
         },
      };

      // Collect global flags
      let mut flt = if vr.phred < 30 { FLT_Q30 } else { 0 };

      // If no non-reference allele has no flags set, then set global flags
      // to union of all allele flags
      if all_res[1..].iter().all(|ar| ar.flt != 0) {
         for ar in all_res[1..].iter() {
            flt |= ar.flt
         }
      }

      if vr.alleles.len() > 1 {
         let mut f = String::new();
         write!(f, "{}\t{}", desc[ref_ix], desc[vr.alleles[1].ix]).ok()?;

         for s in vr.alleles[2..].iter().map(|a| &desc[a.ix]) {
            write!(f, ",{}", s).ok()?;
         }
         write!(f, "\t{}\t{}", vr.phred, Filter(flt)).ok()?;

         // INFO field
         write!(f, "\tDP={}", raw_depth).ok()?;

         for ar in vr.alleles[1..].iter() {
            if desc[ar.ix].len() != desc[ref_ix].len() {
               write!(f, ";INDEL").ok()?;
               break
            }
         }

         // FORMAT field
         write!(f, "\tGT:ADF:ADR:HPL:FQSE:AQ:AFLT:QAVG:FSB:QBS").ok()?;

         // GT field
         write!(f, "\t{}:{}", gt, cts[vr.alleles[0].ix][0]).ok()?;
         // ADF
         for all in vr.alleles[1..].iter() {
            write!(f, ",{}", cts[all.ix][0]).ok()?;
         }
         // ADR
         write!(f, ":{}", cts[vr.alleles[0].ix][1]).ok()?;
         for all in vr.alleles[1..].iter() {
            write!(f, ",{}", cts[all.ix][1]).ok()?;
         }
         // HPL
         write!(f, ":{:.5}", vr.alleles[1].res.freq).ok()?;
         for all in vr.alleles[2..].iter() {
            write!(f, ",{:.5}", all.res.freq).ok()?;
         }
         // FQSE
         write!(f, ":{:.5}", vr.alleles[0].res.se).ok()?;
         for all in vr.alleles[1..].iter() {
            write!(f, ",{:.5}", all.res.se).ok()?;
         }
         // AQ
         write!(f, ":{}", vr.alleles[0].res.lr_test).ok()?;
         for all in vr.alleles[1..].iter() {
            write!(f, ",{}", all.res.lr_test).ok()?;
         }
         // AFLT
         write!(f, ":{}", Filter(all_res[0].flt)).ok()?;
         for ar in &all_res[1..] {
            write!(f, ",{}", Filter(ar.flt)).ok()?;
         }
         // QAVG
         write!(f, ":{:.2}", all_res[0].avg_qual).ok()?;
         for ar in &all_res[1..] {
            write!(f, ",{:.2}", ar.avg_qual).ok()?;
         }
         // FS
         write!(f, ":{:.2e}", all_res[0].fisher_strand).ok()?;
         for ar in &all_res[1..] {
            write!(f, ",{:.2e}", ar.fisher_strand).ok()?;
         }
         // QSB
         write!(f, ":{:.2e}", all_res[0].wilcox).ok()?;
         for ar in &all_res[1..] {
            write!(f, ",{:.2e}", ar.wilcox).ok()?;
         }
         Some(f)
      } else {
         None
      }
   }
}

#[derive(Default, Copy, Clone)]
struct Filter(u32);

impl fmt::Display for Filter {
   fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      if self.0 == 0 {
         write!(f, "PASS")
      } else {
         let mut first = true;
         let mut x = self.0;
         for (s, _) in FLT_STR.iter() {
            if (x & 1) == 1 {
               if !first {
                  write!(f, ";{}", s)?
               } else {
                  first = false;
                  write!(f, "{}", s)?
               }
            }
            x >>= 1;
            if x == 0 {
               break;
            }
         }
         Ok(())
      }
   }
}

pub fn write_vcf_header(sam_hdr: &SamHeader, cfg: &Config) -> io::Result<BufWriter<Writer>> {

   let vcf_output = format!("{}.vcf", cfg.output_prefix());
   let reg = cfg.region();

   let mut vcf_wrt = CompressIo::new()
      .path(vcf_output)
      .ctype(CompressType::Bgzip)
      .bufwriter()?;

   let sample = cfg.sample().unwrap_or("SAMPLE");

   // Write VCF Headers

   writeln!(vcf_wrt, "##fileformat=VCFv4.2")?;
   writeln!(vcf_wrt, "##contig=<ID={},length={}>", sam_hdr.tid2name(reg.tid()), reg.ctg_size())?;
   writeln!(vcf_wrt, "##FILTER=<ID=PASS,Description=\"Site contains at least one allele that passes filters\">")?;
   for (s1, s2) in FLT_STR.iter() { writeln!(vcf_wrt, "##FILTER=<ID={},Description=\"{}\">", s1, s2)?; }
   writeln!(vcf_wrt, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand (high-quality bases)\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand (high-quality bases)\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=HPL,Number=A,Type=Float,Description=\"Estimate of heteroplasmy frequency for alternate alleles\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=FQSE,Number=R,Type=Float,Description=\"Standard errors of allele frequency estimates per allele\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=AQ,Number=R,Type=Float,Description=\"Phred scaled likelihood ratio  for each allele (H0: allele freq is zero)\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=AFLT,Number=R,Type=String,Description=\"Filters per allele\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=FSB,Number=R,Type=Float,Description=\"Fisher test of allele strand bias per allele\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=QAVG,Number=R,Type=Float,Description=\"Average allele base quality scores\">")?;
   writeln!(vcf_wrt, "##FORMAT=<ID=QBS,Number=R,Type=Float,Description=\"Mann-Whitney-Wilcoxon test of minor allele base quality scores\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"95% confidence interval around POS for structural variants\">")?;
   writeln!(vcf_wrt, "##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"95% confidence interval around SVLEN for structural variants\">")?;
   writeln!(vcf_wrt, "##ALT=<ID=DEL, Description=\"Deletion\">")?;
   writeln!(vcf_wrt, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample)?;
   Ok(vcf_wrt)
}
