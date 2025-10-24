use super::ReadDels;

const EMF_FQ_VALID: u8 = 1;
const EMF_LK_VALID: u8 = 2;

pub struct EmParam<'a> {
    dels: &'a ReadDels,
    fq: Box<[f64]>,
    lk: f64,
    valid_flag: u8,
}

impl<'a> EmParam<'a> {
    pub fn new(dels: &'a ReadDels) -> Self {
        // Add extra entry for wildtype
        let n = dels.n_dels() + 1;
        assert!(n > 1, "Invalid size for EmParam");
        let fq = vec![0.0; n].into_boxed_slice();
        Self {
            dels,
            fq,
            lk: 0.0,
            valid_flag: 0,
        }
    }
    
    pub fn dels(&self) -> &'a ReadDels {
        self.dels
    }
    
    pub fn fq(&self) -> Option<&[f64]> {
        if (self.valid_flag & EMF_FQ_VALID) == 0 {
            None
        } else {
            Some(&self.fq)
        }
    }

    pub fn wt_fq(&self) -> Option<f64> {
        self.fq().and_then(|q| q.last().copied())
    }

    pub fn lk(&self) -> Option<f64> {
        if (self.valid_flag & EMF_LK_VALID) == 0 {
            None
        } else {
            Some(self.lk)
        }
    }

    pub fn set_fq(&mut self, fq: &[f64]) {
        self.cp_fq(fq);
        self.valid_flag = (self.valid_flag & !EMF_LK_VALID) | EMF_FQ_VALID;
    }

    pub fn copy_into(&mut self, other: &Self) {
        self.cp_fq(&other.fq);
        self.lk = other.lk;
        self.dels = other.dels;
        self.valid_flag = other.valid_flag;
    }
    
    fn cp_fq(&mut self, fq: &[f64]) {
        assert_eq!(fq.len(), self.fq.len(), "Frequency vec lengths are unequal");
        let s = fq.iter().sum::<f64>();
        assert!((s - 1.0).abs() < 1.0e-12, "Frequency vector not normalized");
        self.fq.copy_from_slice(fq);
    }
    
    pub fn rescale(&mut self, q: f64) {
        assert!((0.0..1.0).contains(&q));
        let n = self.fq.len();
        let scale = (1.0 - q) / (1.0 - self.fq[n - 1]);
        for p in self.fq.iter_mut() {
            *p *= scale
        }
        self.fq[n - 1] = q;
        self.valid_flag &= !EMF_LK_VALID;
    }

    pub fn set_init_fq(&mut self, obs_cts: &[(u64, u64)]) {
        let n = self.fq.len();
        assert_eq!(n, obs_cts.len(), "Incorrect size of obs vector");
        let mut t = 0.0;

        for ((a, b), q) in obs_cts[..n - 1].iter().zip(self.fq[..n - 1].iter_mut()) {
            let p = *a as f64 / *b as f64;
            *q = p;
            t += p;
        }
        self.fq[n - 1] = t;
        
        // Rescale
        for p in self.fq.iter_mut() {
            *p /= 2.0 * t
        }
        
        self.valid_flag = (self.valid_flag & !EMF_LK_VALID) | EMF_FQ_VALID;
    }
    
     pub fn est_freq(&mut self, obs_cts: &[(u64, u64)]) -> anyhow::Result<f64> {
         self.set_init_fq(obs_cts);
         self.lk = self.dels.est_freq(&mut self.fq)?;
         self.valid_flag |= EMF_LK_VALID;
         Ok(self.lk)
     }
     
     pub(super) fn em_max(
         &mut self,
         cts: &mut [f64],
         tmp_dels: &mut Vec<usize>,
         fix_wt: bool,
     ) -> anyhow::Result<f64> {
         
         let lk = self.dels.em_max(cts, tmp_dels, &mut self.fq, fix_wt)?;
         self.lk = lk;
         self.valid_flag |= EMF_LK_VALID;
         Ok(lk)
     }
     
     pub fn calc_profile_like(&mut self, fq_mle: &[f64], np: usize) -> anyhow::Result<Vec<(f64, f64)>> {
         assert!(np > 1);
         self.set_fq(fq_mle);
 
         let mut prof_like = Vec::with_capacity(np + 1);
         let nd = self.fq.len();
         let mut tmp_dels: Vec<usize> = Vec::with_capacity(nd - 1);
         let mut cts = vec![0.0f64; nd];
         let z = 1.0 / (np + 1) as f64;
 
         for i in 0..np {
             let wt_fq = (i + 1) as f64 * z;
             self.rescale(wt_fq);
             self.lk = self.em_max(&mut cts, &mut tmp_dels, true)?;
             prof_like.push((wt_fq, self.lk))
         }
         self.valid_flag |= EMF_LK_VALID;
         
         Ok(prof_like)
     }
     
     pub fn take(self) -> (Box<[f64]>, Option<f64>) {
         let lk = self.lk();
         (self.fq, lk)
     }
}
