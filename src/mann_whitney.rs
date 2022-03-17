use crate::model::N_QUAL;
use crate::stat_funcs::erfc;

pub fn mann_whitney(qcts: &[[usize; N_QUAL]], ix1: usize, ix2: usize) -> Option<f64> {
   // Calculate n1, n2, mean rank for group ix2, and sum of (t_k^3 - t_k) for tied ranks
   let mut next_rank = 1;
   let mut rank_sum = 0.0;
   let (mut n1, mut n2) = (0, 0);
   let mut z_tie = 0.0;
   let ct1 = qcts[ix1];
   let ct2 = qcts[ix2];
   for (&a, &b) in ct1.iter().zip(ct2.iter()) {
      if a + b > 0 {
         n1 += a;
         n2 += b;
         let rank = if a + b > 1 {
            // Tied ranks
            let t = (a + b) as f64;
            z_tie += t * (t * t - 1.0);
            0.5 * ((next_rank + a + b - 1) as f64)
         } else {
            next_rank as f64
         };
         rank_sum += (b as f64) * rank;
         next_rank += a + b;
      }
   }
   let u = rank_sum - ((n2 * (n2 + 1)) >> 1) as f64;
   let n = (n1 + n2) as f64;
   let mu = 0.5 * ((n1 * n2) as f64);
   let v = mu * (n + 1.0 - z_tie / (n * (n - 1.0))) / 6.0;
   if v > 0.0 {
      let z = (u - mu).abs() / v.sqrt();
      let p = erfc(z / std::f64::consts::SQRT_2);
      Some(p)
   } else {
      None
   }
}
