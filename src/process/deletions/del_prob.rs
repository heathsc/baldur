use rstar::{AABB, RTree};

use super::{Deletion, ReadDels};

pub fn calc_del_prob(x: usize, rt: &RTree<Deletion>, fq: &[f64], target_size: usize) -> f64 {
    let mut z = 0.0;

    let bb = AABB::from_point([x as isize, 0]);
    for d in rt.locate_in_envelope_intersecting(&bb) {
        z += fq[d.ix as usize]
    }

    let bb = AABB::from_point([(x + target_size) as isize, 0]);
    for d in rt.locate_in_envelope_intersecting(&bb) {
        z += fq[d.ix as usize]
    }

    z
}

pub struct LikeContrib {
    obs_counts: Vec<(u64, u64)>,
    read_dels: ReadDels,
}

impl LikeContrib {
    pub fn new(obs_counts: Vec<(u64, u64)>, read_dels: ReadDels) -> Self {
        Self {
            obs_counts,
            read_dels,
        }
    }

    pub fn est_freq(&self) -> anyhow::Result<Vec<f64>> {
        info!("Estimating deletion frequencies");
        let nd = self.obs_counts.len();
        assert!(nd > 1);
        let mut fq = self.calc_initial_freq();
        let del_mask = self.read_dels.del_mask();
        let mask = del_mask.dset(0);

        let mut tmp_dels: Vec<usize> = Vec::with_capacity(nd - 1);

        let mut cts = vec![0.0f64; nd];
        // gradient vector
        let mut g = vec![0.0f64; nd - 1];

        let mut converged = false;
        // EM steps
        for it in 0..1000000 {
            cts.fill(0.0);
            g.fill(0.0);
            let mut lk = 0.0;
            for (rd, ct) in self.read_dels.iter() {
                let ct = *ct as f64;

                // For all reads we now add contricutions from the excluded deletions. For these we add 'dummy reads'
                // to compensate for the reads that have been excluded
                tmp_dels.clear();
                let mut excl_freq = 0.0;
                handle_mask(rd.dset(1).iter().copied(), |ix| {
                    tmp_dels.push(ix);
                    excl_freq += fq[ix];
                });

                let excl_freq1 = 1.0 - excl_freq;
                
                if excl_freq > 0.0 {
                    let ct1 = ct / (1.0 - excl_freq);
                    for i in tmp_dels.drain(..) {
                        cts[i] += ct1 * fq[i];
                    }
                }

                if let Some(ix) = rd.obs_index() {
                    // Get read contributions for observed deletions
                    cts[ix] += ct;
                    assert!(fq[ix] < excl_freq1);
                    let z = fq[ix] / excl_freq1;

                    lk += z.ln() * ct;
                    g[ix] += ct / fq[ix];
                    
                    // Since we are dividing by 1 - exlc_freq, we also have ontributions to the gradients of
                    // the excluded deletions
                    let zg = ct / excl_freq1;
                    handle_mask(rd.dset(1).iter().copied(), |i| {
                        g[i] += zg;
                    });
                } else {
                    // Get read contributions where no deletion was observed
                    //
                    // Go through all deletions that could have been observed if the read was full length.
                    // This means all deletions that are not covered by the read and do not overlap the
                    // physical start of the read.
                    // To get this we do the bitwise or of the covered and excluded deletions, then
                    // the bitwise xor of this with the mask of all deletions.
                    tmp_dels.clear();
                    let mut z = fq[nd - 1];
                    // eprintln!("Adding wt {} {}", nd - 1, z);
                    tmp_dels.push(nd - 1);
                    handle_mask(
                        rd.dset(0)
                            .iter()
                            .zip(rd.dset(1).iter())
                            .zip(mask.iter())
                            .map(|((m1, m2), m3)| (m1 | m2) ^ m3),
                        |ix| {
                            tmp_dels.push(ix);
                            z += fq[ix];
                        },
                    );

                    let y = ct / z;
                    for i in tmp_dels.iter().copied() {
                        cts[i] += fq[i] * y;
                    }
                    assert!(z - excl_freq1 < 1.0e-12, "{z} {excl_freq}");
                    let zz = z / excl_freq;

                    lk += zz.ln() * ct;

                    let mut zc = 0.0;
                    handle_mask(rd.dset(0).iter().copied(), |ix| {
                        g[ix] -= y;
                        zc += fq[ix];
                    });
                
                    
                    let zg = -ct * zc / (z * excl_freq1);
                    handle_mask(rd.dset(1).iter().copied(), |ix| {
                        g[ix] += zg;
                    });
                }
            }

            // Get new estimates of frequencies
            let mut ss = 0.0;
            let z = cts.iter().sum::<f64>();

            for (c, f) in cts.iter().zip(fq.iter_mut()) {
                let t = *f;

                *f = *c / z;
                ss += (t - *f).powi(2);
            }
            let delta = (ss / nd as f64).sqrt();
            let gsq = g.iter().map(|x| x * x).sum::<f64>();
            let grms = (gsq / (nd - 1) as f64).sqrt();

            eprintln!(
                "{it}\t{lk}\t{delta}\t{grms}\t{z}\t{}\t{}",
                fq[0],
                fq[nd - 1]
            );
            if delta < 1.0e-12 {
                debug!(
                    "Deletion frequency estimation: convergence criteria achieved after {} iterations",
                    it + 1
                );
                converged = true;
                break;
            }
        }
        if !converged {
            Err(anyhow!("Deletion frequency estimates did not converge"))
        } else {
            for (i, (f, (c, n))) in fq.iter().zip(self.obs_counts.iter()).enumerate() {
                let z = *c as f64 / *n as f64;
                eprintln!("{i}\t{f}\t{z}\t{c}\t{n}\t{}", g[i]);
            }
            Ok(fq)
        }
    }

    fn calc_initial_freq(&self) -> Vec<f64> {
        debug!("Get initial frequency estimates");
        let n = self.obs_counts.len();
        let mut fq = Vec::with_capacity(n);
        let mut t = 0.0;
        for (a, b) in self.obs_counts[..n - 1].iter() {
            let p = *a as f64 / *b as f64;
            fq.push(p);
            t += p;
        }
        if t < 1.0 {
            fq.push(1.0 - t)
        } else {
            // Rescale
            let p = 1.0 / n as f64;
            fq.push(p);
            t += p;
            for p in fq.iter_mut() {
                *p /= t
            }
        }
        fq
    }
}

fn handle_mask<I, F>(mask_itr: I, mut f: F)
where
    I: Iterator<Item = u64>,
    F: FnMut(usize),
{
    for (j, mut x) in mask_itr.enumerate() {
        let mut i = 0;
        while x != 0 {
            if (x & 1) == 1 {
                let ix = (j << 6) | i;
                f(ix)
            }
            i += 1;
            x >>= 1
        }
    }
}
