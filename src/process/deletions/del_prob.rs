use rstar::{AABB, RTree};

use super::{Deletion, ReadDelSet, ReadDels};

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

impl ReadDels {
    pub fn est_freq(&self, obs_cts: &[(u64, u64)]) -> anyhow::Result<Vec<f64>> {
        info!("Estimating deletion frequencies");
        let nd = obs_cts.len();
        assert!(nd > 1);
        let mut fq = calc_initial_freq(obs_cts);

        let mut tmp_dels: Vec<usize> = Vec::with_capacity(nd - 1);

        let mut cts = vec![0.0f64; nd];
        // gradient vector
        // let mut g = vec![0.0f64; nd - 1];

        let mut converged = false;

        let mut old_lk: Option<f64> = None;
        // EM steps
        for _ in 0..1000000 {
            let lk = self.em_update(&mut cts, &mut tmp_dels, &mut fq);

            let diff = old_lk.take().map(|x| lk - x);

            // eprintln!("{it}\t{lk}\t{diff:?}\t{}\t{}", fq[0], fq[nd - 1]);

            // Check convergence
            if diff.map(|d| d < 1.0e-8).unwrap_or(false) {
                converged = true;
                break;
            }

            old_lk = Some(lk);
        }

        if !converged {
            Err(anyhow!("Deletion frequency estimates did not converge"))
        } else {
            Ok(fq)
        }
    }

    fn em_update(&self, cts: &mut [f64], tmp_dels: &mut Vec<usize>, fq: &mut [f64]) -> f64 {
        // Initialize tmp storage and likelihood
        let mut lk = 0.0;
        cts.fill(0.0);

        for (rd, ct) in self.iter() {
            let ct = *ct as f64;

            // Handle contributions from excluded deletions. Get total frequency
            // of excluded alleles
            let excl_fq = em_excl_contrib(rd, ct, cts, tmp_dels, fq);

            if let Some(ix) = rd.obs_index() {
                lk += em_handle_observed_del(ct, ix, cts, fq, excl_fq)
            } else {
                lk += em_handle_unobserved_del(rd, ct, cts, tmp_dels, fq, excl_fq)
            }
        }

        // Update frequency estimates
        let z = cts.iter().sum::<f64>();

        for (c, f) in cts.iter().zip(fq.iter_mut()) {
            *f = *c / z;
        }

        lk
    }
}

fn calc_initial_freq(obs_cts: &[(u64, u64)]) -> Vec<f64> {
    debug!("Get initial frequency estimates");
    let n = obs_cts.len();
    let mut fq = Vec::with_capacity(n);
    let mut t = 0.0;
    for (a, b) in obs_cts[..n - 1].iter() {
        let p = *a as f64 / *b as f64;
        fq.push(p);
        t += p;
    }
    fq.push(t);

    // Rescale
    for p in fq.iter_mut() {
        *p /= 2.0 * t
    }
    fq
}
/// Get read contributions for observed deletions
fn em_handle_observed_del(ct: f64, ix: usize, cts: &mut [f64], fq: &[f64], excl_fq: f64) -> f64 {
    cts[ix] += ct;
    (fq[ix] / (1.0 - excl_fq)).ln() * ct
}

/// Get read contributions where no deletion was observed
///
/// Go through all deletions that could have been observed if the read was full length.
/// This means all deletions that are not covered by the read and do not overlap the
/// physical start of the read.
/// To get this we do the bitwise or of the covered and excluded deletions, then
/// the bitwise xor of this with the mask of all deletions.
fn em_handle_unobserved_del(
    rd: &ReadDelSet,
    ct: f64,
    cts: &mut [f64],
    tmp_dels: &mut Vec<usize>,
    fq: &[f64],
    excl_fq: f64,
) -> f64 {
    tmp_dels.clear();
    let nd = cts.len();
    let mut z = fq[nd - 1];

    tmp_dels.push(nd - 1);
    handle_mask(rd.dset(2).iter().copied(), |ix| {
        tmp_dels.push(ix);
        z += fq[ix];
    });

    let y = ct / z;
    for i in tmp_dels.iter().copied() {
        cts[i] += fq[i] * y;
    }

    (z / (1.0 - excl_fq)).ln() * ct
}

/// Get contributions from deletions are excluded (due to them overlapping the physical start of the read)
fn em_excl_contrib(
    rd: &ReadDelSet,
    ct: f64,
    cts: &mut [f64],
    tmp_dels: &mut Vec<usize>,
    fq: &[f64],
) -> f64 {
    tmp_dels.clear();

    let mut excl_fq = 0.0;
    handle_mask(rd.dset(1).iter().copied(), |ix| {
        tmp_dels.push(ix);
        excl_fq += fq[ix];
    });
    
    for i in tmp_dels.drain(..) {
        cts[i] += ct * fq[i];
    }

    excl_fq
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
