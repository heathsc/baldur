use nalgebra::{DVector, DMatrix};
use log::Level::Trace;
use crate::stat_funcs::chisq1;

pub const N_QUAL: usize = 64;
const MAX_PHRED: u8 = 40;
const MAX_ITER: usize = 5000;
const ZERO_LIM: f64 = 1.0e-8;

// Convergence criteria for maximization iterations
const CONVERGENCE_CRITERION: f64 = 1.0e-12;

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum State {
   In, Test, Out
}

impl Default for State {
   fn default() -> Self { Self::Out }
}

#[derive(Default, Copy, Clone)]
pub struct AlleleRes {
   pub freq: f64, // Allele frequency estimate
   pub se: f64, // Standard error of frequency estimate
   pub lr_test: u8, // Likelihood ratio test for allele frequency being > 0.0
   pub flag: bool,
}

pub struct ModelRes {
   pub alleles: Vec<AlleleRes>,
   pub phred: u8,
}

struct ObsCount {
   e: f64, // Error probability for observations
   total: f64, // Total counts at this error level
   cts: Vec<(f64, f64)>, // allele counts at this error level
}

/// Collect vector of observations for quality values where at least one of the active alleles has observations
fn collect_obs(qcts: &[[usize; N_QUAL]], qual_model: &[f64; N_QUAL], counts: &mut [f64], index: &[usize]) -> Vec<ObsCount> {
   let mut flag = vec!(false; qcts.len());
   for &ix in index.iter() { flag[ix] = true }
   let obs_set: Vec<_> = (0..N_QUAL)
      .filter(|&i| index.iter().any(|&ix| qcts[ix][i] > 0))
      .map(|i| {
         let mut total = 0.0;
         let mut cts: Vec<_> = qcts.iter().zip(flag.iter()).enumerate().map(|(ix, (ct, &fg))| {
            if fg {
               let z = ct[i] as f64;
               counts[ix] += z;
               total += z;
               (z, 0.0)
            } else {
               (0.0, 0.0)
            }
          }).collect();
         cts.push((0.0, 0.0));
         let e = qual_model[i];
         ObsCount{e, total, cts}
      }).collect();
   obs_set
}

/// Check list of active parameters.  Regenerate index and rescale frequencies if the
/// number of active parameters has changed
fn check_active_set(flag: &[State], index: &mut Vec<usize>, freq: &mut[f64]) -> usize {
   let n_active = flag.iter().fold(0, |s, &st| if st != State::Out { s + 1 } else { s });

   if n_active != index.len() {
      // Active parameters have changed.  Adjust frequency estimates
      assert!(n_active > 0);
      // Update index vec.  Collect sum of active frequency estimates
      let mut tot = 0.0;
      index.clear();
      for (ix, (&fg, fq)) in flag.iter().zip(freq.iter_mut()).enumerate() {
         if fg != State::Out {
            tot += *fq;
            index.push(ix);
         }
      }
      // Rescale frequencies
      for fq in freq.iter_mut() {
         *fq /= tot
      }
   }

   // If first allele in index is in State::Test, try and swap with another allele.
   // This is because first allele in index is not directly in model (as it is estimated as 1 - the sum of the
   // other frequencies) so it is harder to check the constraints
   if flag[index[0]] == State::Test {
      if let Some((i, &j)) = index[1..].iter().enumerate().find(|(_, &j)| flag[j] == State::In) {
         index[i + 1] = index[0];
         index[0] = j
      }
   }

   n_active
}

/// Updates prob. values for each observation and returns log likelihood
fn calc_like(obs_set: &mut [ObsCount], freq: &[f64], index: &[usize]) -> f64 {
   let mut like = 0.0;
   for obs in obs_set.iter_mut() {
      let e = obs.e;
      let ix0 = index[0];
      let p0 = freq[ix0] * (1.0 - e) + (1.0 - freq[ix0]) * e;
      let (cts, ps) = &mut obs.cts[ix0];
      *ps = p0;
      like += *cts * p0.ln();
      let mut flt = obs.total - *cts;
      for &ix in index[1..].iter() {
         let p = freq[ix] * (1.0 - e) + (1.0 - freq[ix]) * e;
         let (cts, ps) = &mut obs.cts[ix];
         *ps = p;
         like += *cts * p.ln();
         flt -= *cts;
      }
      like += flt * e.ln();
   }
   like
}

/// Updates prob. values for each observation and
/// calculate score vector (vector of first partial derivatives of log likelihood)
/// Returns current value of log likelihood
fn calc_score_vec(obs_set: &mut [ObsCount], freq: &[f64], index: &[usize], score: &mut[f64]) -> f64 {
   let mut like = 0.0;
   for sc in score[0..index.len() - 1].iter_mut() { *sc = 0.0 }
   for obs in obs_set.iter_mut() {
      let e = obs.e;
      let ix0 = index[0];
      let z = 1.0 - 2.0 * e;
      let p0 = freq[ix0] * (1.0 - e) + (1.0 - freq[ix0]) * e;
      let (cts, ps) = &mut obs.cts[ix0];
      *ps = p0;
      like += *cts * p0.ln();
      let mut flt = obs.total - *cts;
      let konst0 = *cts / p0;
      for (&ix, sc) in index[1..].iter().zip(score.iter_mut()) {
         let p = freq[ix] * (1.0 - e) + (1.0 - freq[ix]) * e;
         let (cts, ps) = &mut obs.cts[ix];
         *ps = p;
         like += *cts * p.ln();
         *sc += z * (*cts / p - konst0);
         flt -= *cts;
      }
      like += flt * e.ln();
   }
   like
}

/// Calculate Fisher Information matrix (Expected matrix of second derivatives)
/// All off diagonal elements are the same, so we just calculate the diagonal
/// and store in info_diag, and return the off diagonal element
fn calc_fisher_inf(obs_set: &[ObsCount], index: &[usize], info_diag: &mut [f64]) -> f64 {
   for elem in info_diag.iter_mut() { *elem = 0.0 }
   let mut off = 0.0;
   for obs in obs_set.iter() {
      let e = obs.e;
      let total = index.iter().fold(0.0, |s, &ix| s + obs.cts[ix].0 );
      let ix0 = index[0];
      let z = (1.0 - 2.0 * e) * (1.0 - 2.0 * e);
      let q0 = total / obs.cts[ix0].1;
      off += z * q0;
      for (&ix, elem) in index[1..].iter().zip(info_diag.iter_mut()) {
         *elem += z * (q0 +  total / obs.cts[ix].1);
      }
   }
   off
}

/// Calculate parameter change vector.
/// `info_diag` has on input the diagonal elements of the fisher information matrix J
/// and on output has the change vector.
/// `info_off` has the off diagonal elements (they are all identical).
/// `score` has on input the score vector S
fn calc_delta<'a>(info_diag: &'a mut [f64], info_off: f64, score: &[f64]) -> &'a [f64] {
   // Compute S / J
   match score.len() {
      0 => panic!("calc_delta called with empty score vector"),
      1 => info_diag[0] = score[0] / info_diag[0], // 1D case
      2 => {
         // 2 x 2 case.
         // Calculate determinant
         let det = info_diag[0] * info_diag[1] - info_off * info_off;
         // Off diagonal element of inverse of J
         let off = -info_off / det;
         // And put elements of S / J into info_diag
         let d0 = score[0] * info_diag[1] / det + score[1] * off;
         let d1 = score[0] * off + score[1] * info_diag[0] / det;
         info_diag[0] = d0;
         info_diag[1] = d1;
      },
      n => {
         // n x n case
         let mut svec = DVector::from_column_slice(&score);
         let jac = DMatrix::from_fn(n, n, |i, j| {
            if i == j { info_diag[i] }
            else { info_off }
         });
         let chol = jac.cholesky().expect("J not positive definite!");
         chol.solve_mut(&mut svec);
         for i in 0..n {
            info_diag[i] = svec[i]
         }
      },
   }
   &info_diag[0..score.len()]
}

/// Calculate asymptotic standard errors for parameter estimates by inverting
/// Fisher information matrix supplied as the diagonal elements and one off-diagonal
/// element (as they are all identical for this problem
fn calc_se(info_diag: &[f64], info_off: f64, se: &mut [f64], index: &[usize]) {
   match info_diag.len() {
      0 => panic!("calc_se() called with too few parameters"),
      1 => se[index[1]] = (1.0 / info_diag[0]).sqrt(),
      2 => {
         let det = info_diag[0] * info_diag[1] - info_off * info_off;
         se[index[1]] = (info_diag[1] / det).sqrt();
         se[index[2]] = (info_diag[0] / det).sqrt();
      },
      n => {
         let jac = DMatrix::from_fn(n, n, |i, j| {
            if i == j { info_diag[i] }
            else { info_off }
         });
         let chol = jac.cholesky().expect("J not positive definite!");
         let vcov = chol.inverse();
         for (&v, &ix) in vcov.diagonal().iter().zip(index[1..].iter()) {
            se[ix] = v.sqrt()
         }
      },
   }
}

/// Find ML estimates of allele frequencies using Fisher's scoring method
///
/// `obs_set` are the observations per allele at each observed quality level
///
/// `freq` at input has the initial values for the frequencies, and on exit has
/// the ML frequency estimates.
///
/// `index` has at input the indices of alleles to consider, and at output has
/// the indices of alleles retained (as nonzero) in the model
///
/// `se` at output has the asymptotic standard errors for the frequencies
///
fn ml_estimation(obs_set: &mut [ObsCount], freq: &mut[f64], mut se: Option<&mut[f64]>, index: &mut Vec<usize>) -> f64 {

   // (Maximum) number of parameters is one less than the number of active alleles (because the
   // frequencies of all active alleles have to sum to 1)

   assert!(!index.is_empty());
   let np = index.len() + 1;
   let n_all = freq.len();

   // Storage for score vector
   let mut score = vec!(0.0; np);

   // Storage for Diagonal of Fisher Information matrix
   let mut fisher_inf = vec!(0.0; np);

   // Storage for temporary frequency estimates
   let mut fq = vec!(0.0; n_all);

   // Allele flags
   let mut flag = vec!(State::default(); n_all);
   for &ix in index.iter() {
      flag[ix] = State::In
   }

   // Main iteration loop
   let mut n_active = check_active_set(&flag, index, freq); // Number of allele currently active in model
   let mut like: f64 = 0.0; // Current log likeihood
   for it in 0..MAX_ITER {
      loop {
         // Update observation probs.
         like = calc_score_vec(obs_set, freq, index, &mut score[..n_active - 1]);
         if n_active <= 1 { break }

         // If any of the alleles are in State::Test (which means that their current value os zero)
         // Check to see if the relevant element of the score vector is negative.  If so then we
         // know the maximum will be < 0, so we can remove the allele permanently from the model
         let mut changed = false;
         for (&ix, &d) in index[1..].iter().zip(score.iter()) {
            if flag[ix] == State::Test && d < 0.0 {
               flag[ix] = State::Out;
               changed = true;
            }
         }
         if !changed { break }
         // Update active set
         n_active = check_active_set(&flag, index, freq);
      }
      trace!("like: {}, n_active = {}, np = {}", like, n_active, np);
      if n_active <= 1 { break }

      // Calculate Fisher information matrix
      let off = calc_fisher_inf(obs_set, index, &mut fisher_inf[..n_active - 1]);

      // Calculate update vector for frequencies
      let delta = calc_delta(&mut fisher_inf[..n_active - 1], off, &score[..n_active - 1]);

      // Find largest values between 0 and 1 for alpha where freq + alpha * delta is >=0
      // for all frequencies
      let mut alpha = check_freq(delta, freq, index);

      // The scoring update can result in the likelihood decreasing.  We start by adding alpha * delta
      // to the current frequencies and testing the likelihood.  If the likelihood has decreased
      // we divide alpha by 2 and retry while alpha > 0
      let like1 = like;
      while alpha > 0.0 {
         update_freq(delta, alpha, freq, &mut fq, index);
         like = calc_like(obs_set, &fq, index);
         if like < like1 {
            alpha *= 0.5
         } else {
            for (&f, fp) in fq.iter().zip(freq.iter_mut()) { *fp = f }
            break
         }
      }
      if log_enabled!(Trace) {
         eprint!("It:{}\t{}\t{}", it, like, alpha);
         for &f in freq.iter() { eprint!("\t{:.6}", f) }
         eprintln!();
      }
      let mut changed = false;

      for (f, fg) in freq.iter_mut().zip(flag.iter_mut()) {
         match *fg {
            State::In => if *f <= ZERO_LIM {
               *fg = State::Test;
               *f = 0.0;
               changed = true
            },
            State::Test =>  if *f <= ZERO_LIM {
               *fg = State::Out;
               *f = 0.0;
               changed = true
            } else { *fg = State::In },
            _ => (),
         }
      }
      if changed {
         // Update active set
         n_active = check_active_set(&flag, index, freq);
      } else if like - like1 < CONVERGENCE_CRITERION { break }
   }

   // To get standard errors for all non-zero frequency estimates we will add back in a parameter
   // with value 0 to the first position of index (and then remove it afterwards).
   if let Some(se) = se.take() {
      trace!("Obtaining SE estimates");
      // Find an unused parameter
      let k = flag.iter().enumerate().find(|(_, &fg)| fg == State::Out).map(|(ix, _)| ix).expect("Couldn't find place for dummy parameter");
      // Swap with first parameter (index[0]); push the previous first parameter onto the back of index
      let ix0 = index[0];
      index[0] = k;
      freq[k] = 0.0;
      index.push(ix0);
      // Update observation probabilities
      calc_like(obs_set, freq, index);
      // Calculate Fisher Information
      let fi = &mut fisher_inf[..n_active];
      let off = calc_fisher_inf(obs_set, index, fi);
      calc_se(fi, off, se, index);
      index[0] = index.pop().unwrap();
      if log_enabled!(Trace) {
         for (&f, &s) in freq.iter().zip(se.iter()) { trace!("{:.6}\t{:.6}", f, s) }
      }
   }
   like
}

/// Find the maximum value of alpha between 0 and 1 where F + alpha * D
/// (F = current frequency estimates and D = change vector) results in all
/// elements of F being >=0.  Return alpha.
fn check_freq(delta: &[f64], freq: &[f64], index: &[usize]) -> f64 {
   // delta for first parameter (which is just 1 - the others)
   let d0 = -delta.iter().sum::<f64>();
   // Current frequency for first parameter
   let f0 = freq[index[0]];
   let mut alpha = if f0 + d0 < 0.0 { -f0 / d0 } else { 1.0 };
   for (&ix, &d) in index[1..].iter().zip(delta.iter()) {
      let f = freq[ix];
      if f + d < 0.0 {
         alpha = alpha.min(-f / d)
      }
   }
   alpha
}

/// Update frequency estimates using `new_freq` = `freq` + `alpha` * `delta`
fn update_freq(delta: &[f64], alpha: f64, freq: &[f64], new_freq: &mut[f64], index: &[usize]) {
   // delta for first parameter (which is just 1 - the others)
   let d0 = -delta.iter().sum::<f64>();
   // update first parameter
   let ix0 = index[0];
   new_freq[ix0] = (freq[ix0] + alpha * d0).max(0.0);
   for (&ix, &d) in index[1..].iter().zip(delta.iter()) {
      new_freq[ix] = (freq[ix] + alpha * d).max(0.0);
   }
}

/// Find ML estimates of allele frequencies.  Counts by quality for each possible allele are
/// given in `qcts`.  `alls` has the indexes in `qcts` of the subset of alleles that should be
/// considered (the frequencies of any other alleles is assumed to be zero)
///
pub fn freq_mle(alls: &[usize], qcts: &[[usize; N_QUAL]], qual_model: &[f64; N_QUAL]) -> ModelRes {
   let n_active = alls.len();
   let mut n_alls = qcts.len();
   let ref_ix = alls[0];

   // Sanity checks
   assert!(n_active > 0 && n_active <= n_alls && alls.iter().all(|&i| i < n_alls));

   // Make local copy of alls
   let mut alls = alls.to_vec();

   // If all alleles are in alls, add a dummy allele for the estimation of s.e.
   if n_active == n_alls { n_alls += 1 }

   let mut all_res = vec!(AlleleRes::default(); n_alls);
   for &ix in alls.iter() { all_res[ix].flag = true }

   let phred = if n_active == 1 {
      // One allele case
      all_res[alls[0]] = AlleleRes { freq: 1.0, lr_test: MAX_PHRED, se: 0.0, flag: true };
      MAX_PHRED
   } else {
      // Multiple alleles

      // Make list of all required values for the iterations
      let mut counts = vec!(0.0; n_alls);
      let mut obs_set = collect_obs(qcts, qual_model, &mut counts, &alls);
      let total_obs :f64 = counts.iter().sum();

      // Initial frequency estimates
      let mut freq: Vec<_> = counts.iter().map(|&ct| ct / total_obs).collect();
      // Storage for standard error estimates
      let mut se = vec!(0.0; n_alls);

      trace!("Initial ML maximization");
      let log_like = ml_estimation(&mut obs_set, &mut freq, Some(&mut se), &mut alls);

      let n_active = alls.len();

      // If multiple alleles have been retained, calculate LR test for each such non-reference allele in turn
      if n_active > 1 {
         trace!("Obtaining LR ratios");
         let mut max_ph = 0;
         for &i in alls.iter() {
            // Check if we need to do calculate LR.  If the approximate Z test indicates a p value
            // that corresponds with the phred score > N_QUAL then we simply set the phred to MAX_PHRED
            let s = se[i];
            assert!(s > 0.0);
            let z = freq[i] / s;
            let ph = if z < 5.0 {
               // Zero the frequency for each retained allele in turn, rescaling the others
               let z = 1.0 / (1.0 - freq[i]);
               let mut fq1: Vec<_> = freq.iter().enumerate()
                  .map(|(k, &fq)| if k == i { 0.0 } else { fq * z }).collect();

               let mut alls1: Vec<usize> = alls.iter()
                  .filter(|&ix| *ix != i).copied().collect();
               let log_like1 = ml_estimation(&mut obs_set, &mut fq1, None, &mut alls1);
               let lr = (log_like - log_like1).max(0.0);
               if lr < 13.0 { // A LR of 13 gives a phred score of >64, which is the maximum quality value we allow
                  ((chisq1(2.0 * lr).log10() * -10.0).round() as u8).min(MAX_PHRED)
               } else { MAX_PHRED }
            } else { MAX_PHRED };
            trace!("Allele {}\tFreq: {}\tLR: {}\tSE: {}", i, freq[i], ph, se[i]);
            all_res[i].lr_test = ph;
            all_res[i].freq = freq[i];
            all_res[i].se = se[i];
            if ref_ix != i && ph > max_ph {
               max_ph = ph;
            }
         }
         max_ph
      } else {
         let ar = &mut all_res[alls[0]];
         ar.freq = 1.0;
         ar.se = se[alls[0]];
         ar.lr_test = MAX_PHRED;
         ar.flag = true;
         MAX_PHRED
      }
   };

   ModelRes {
      alleles: all_res,
      phred,
   }
}

pub fn setup_qual_model() -> [f64; 64] {
   let mut a = [0.0; 64];
   let l10: f64 = f64::ln(10.0);
   for (q, ap) in a.iter_mut().enumerate() {
      *ap = (-0.1 * (q as f64) * l10).exp();
   }
   a
}

