use std::collections::HashSet;

use log::Level::Trace;

use crate::{
    alleles::{get_large_deletions, AllDesc, LargeDeletion, Trunc},
    depth::*,
    model::N_QUAL,
    process::ProcWork,
    vcf::{VcfCalc, VcfRes},
};

/// Estimate single base allele frequencies
pub(crate) fn estimate_single_base_freq(
    pw: &mut ProcWork,
    res: &mut Vec<VcfRes>,
    vc: &VcfCalc,
    ix: Option<usize>,
    iy: Option<usize>,
) {
    debug!(
        "Estimating single base allele frequencies ({:?} - {:?})",
        ix, iy
    );
    let reg = vc.cfg.region();
    let ref_seq = &pw.ref_seq;
    let counts = &mut pw.depth.counts;
    let ix = ix.unwrap_or(0);
    let iy = iy.unwrap_or(counts.len());
    for (i, dep) in counts[ix..iy].iter_mut().enumerate() {
        let x = i + reg.start() + ix;
        trace!("At: {}", x + 1);

        // Check if we have one or more insertion alleles at this locus to consider
        let opt_ins = pw.depth.ins_hash.get(&x).filter(|&ins| {
            ins.hash
                .iter()
                .map(|(_, &i)| &ins.alleles[i])
                .any(|ct| ct.cts[0] > 1 && ct.cts[1] > 1)
        });

        // Frequency estimation
        let vr = if let Some(ins_all) = opt_ins {
            // Handle insertion alleles

            // Collect regular and insertion alleles
            let mut adesc: Vec<AllDesc> = Vec::with_capacity(8);
            let ref_ix = ref_seq[x].base();

            // The counts for the regular alleles
            let mut new_cts = Vec::with_capacity(6);
            let mut new_qcts = Vec::with_capacity(6);

            let mut add_allele = |x: u8| {
                adesc.push(AllDesc::make(vec![x], None, Trunc::No));
                new_cts.push(dep.cts[x as usize]);
                new_qcts.push(dep.qcts[x as usize]);
            };

            // Add regular bases + del
            add_allele(ref_ix);
            for b in (0..5).filter(|&x| x != ref_ix) {
                add_allele(b)
            }

            // And the indel flags
            let mut indel_flag = vec![false, false, false, false, true];

            // Add possible insertion alleles and their counts
            for (all, ix) in ins_all.hash.iter().map(|(all, &i)| (all, i)) {
                let ct = &ins_all.alleles[ix];
                if ct.cts[0] > 1 && ct.cts[1] > 1 {
                    adesc.push(AllDesc::make(all.to_vec(), Some(ix), Trunc::No));
                    new_cts.push(ct.cts);
                    new_qcts.push(ct.qcts);
                    indel_flag.push(true);
                }
            }

            if log_enabled!(Trace) {
                trace!("Insertion counts: {}", x + 1);
                for (a, c) in adesc.iter().zip(new_cts.iter()) {
                    trace!("{}\t{}\t{}", a, c[0], c[1]);
                }
            }

            // Get mle estimates of alleles frequencies
            let mut vr = vc.get_mallele_freqs(x, &new_cts, &new_qcts, &indel_flag);

            if log_enabled!(Trace) {
                trace!("Insertion freq. estimates: {}", x + 1);
                for ar in vr.alleles.iter() {
                    trace!("{}\t{}\t{}", &adesc[ar.ix], ar.res.freq, ar.res.lr_test);
                }
            }
            dep.cts = new_cts;
            dep.qcts = new_qcts;
            vr.adesc = Some(adesc);
            vr
        } else {
            // Handle standard alleles
            if log_enabled!(Trace) {
                trace!("Standard counts: {}", x + 1);
                for (a, c) in vc.all_desc[..5].iter().zip(dep.cts.iter()) {
                    trace!("{}\t{}\t{}", a, c[0], c[1]);
                }
            }
            let mut vr = vc.get_allele_freqs(x, &dep.cts[..5], &dep.qcts[..5]);
            if log_enabled!(Trace) {
                trace!("Standard freq. estimates: {}", x + 1);
                for ar in vr.alleles.iter() {
                    trace!(
                        "{}\t{}\t{}",
                        &vc.all_desc[ar.ix],
                        ar.res.freq,
                        ar.res.lr_test
                    );
                }
            }
            vr.adesc = None;
            vr
        };

        res.push(vr);
    }
}

fn recalc_linear_block(pw: &mut ProcWork, start: usize, stop: usize, read_hash: &HashSet<usize>) {
    let depth = &mut pw.depth;
    let dalign = depth.dalign.as_mut().unwrap();
    let ins_hash = &mut depth.ins_hash;
    let dep_cts = &mut depth.counts;

    // Zero allele counts
    for (x, ct) in (start..=stop).zip(dep_cts[start..=stop].iter_mut()) {
        for c in ct.cts.iter_mut() {
            *c = [0; 2]
        }
        for qc in ct.qcts.iter_mut() {
            *qc = [0; N_QUAL]
        }
        let ihash = ins_hash.get_mut(&x);
        if let Some(ins) = ihash {
            for all in ins.alleles.iter_mut() {
                all.cts = [0; 2];
                all.qcts = [0; N_QUAL]
            }
        }
    }
    // Collect individual read data
    for (_, da) in dalign
        .iter()
        .enumerate()
        .filter(|(k, _)| !read_hash.contains(k))
    {
        da.process_allele(start, stop, |x, c, qual| {
            if (c as usize) % HIDX == BASES_INS {
                let st = if c as usize > HIDX { 1 } else { 0 };
                let ins_alls = ins_hash.get_mut(&x).expect("Missing insertion!");
                let (ix, q1) = da
                    .ins_hash
                    .get(&x)
                    .expect("Missing insertion index for read");
                let all = &mut ins_alls.alleles[*ix];
                all.cts[st] += 1;
                all.qcts[*q1 as usize] += 1;
            } else {
                let (ix, st) = base_to_ix_strand(c as usize);
                let ct = &mut dep_cts[x];
                if ix < 4 {
                    // Don't put back residual deletions
                    ct.cts[ix][st] += 1;
                    ct.qcts[ix][qual as usize] += 1;
                }
            }
        });
    }
}

fn recalc_block(pw: &mut ProcWork, x: usize, y: usize, active_del: &[&LargeDeletion]) {
    if log_enabled!(Trace) {
        trace!("Recalculating counts for block {} - {}", x, y);
        for d in active_del {
            trace!("Active del: {}", d);
        }
    }
    // Prepare HashSet with reads in active deletion set
    let nreads: usize = active_del.iter().map(|d| d.n()).sum();
    let mut read_hash = HashSet::with_capacity(nreads);
    for d in active_del.iter() {
        for &ix in d.reads.iter() {
            read_hash.insert(ix);
        }
    }

    let ref_len = pw.ref_seq.len();
    let x = x % ref_len;
    let y = y % ref_len;
    if y < x {
        recalc_linear_block(pw, x, ref_len - 1, &read_hash);
        recalc_linear_block(pw, 0, y, &read_hash);
    } else {
        recalc_linear_block(pw, x, y, &read_hash);
    }
}

// Recalculated allele counts after accounting for reads showing large deletions
fn recalc_allele_counts(pw: &mut ProcWork, large_del: &[LargeDeletion]) {
    trace!("Recalculating allele counts under large deletion(s)");

    // Active set of deletions
    let mut active = vec![&large_del[0]];

    let mut x = active[0].start;
    let mut y = active.iter().map(|d| d.end()).min().unwrap();
    for del in large_del[1..].iter() {
        while del.start > x {
            let y1 = if y >= del.start { del.start - 1 } else { y };
            recalc_block(pw, x, y1, &active);
            x = y1 + 1;
            active.retain(|d| d.end() >= x);
            if active.is_empty() {
                break;
            }
            y = active.iter().map(|d| d.end()).min().unwrap();
        }
        active.push(del);
        y = y.min(del.end());
    }
    // Handle remaining deletions
    while !active.is_empty() {
        y = active.iter().map(|d| d.end()).min().unwrap();
        recalc_block(pw, x, y, &active);
        x = y + 1;
        active.retain(|d| d.end() >= x);
    }
}

pub(crate) fn process_large_deletions(
    pw: &mut ProcWork,
    res: &mut Vec<VcfRes>,
    vc: &VcfCalc,
) -> Vec<LargeDeletion> {
    debug!("Checking for large deletions");

    let ref_len = vc.cfg.region().len();
    let is_del = |vr: &VcfRes| vr.alleles.iter().any(|ar| ar.ix == 4);
    let mut del_vec = Vec::new();

    let mut block_start = None;
    let mut spanning_block = None;
    let mut blocks = Vec::new();
    for (ix, vr) in res.iter().enumerate() {
        if is_del(vr) {
            if block_start.is_none() {
                block_start = Some(ix)
            }
        } else if let Some(i) = block_start.take() {
            if res[i].x == 0 {
                assert_eq!(i, 0);
                spanning_block = Some(ix - 1)
            } else {
                blocks.push((i - 1, ix - 1));
            }
        }
    }
    let l = res.len();

    match (block_start, spanning_block) {
        (Some(i), Some(k)) => blocks.push((i - 1, k)),
        (Some(i), None) => blocks.push((i - 1, l - 1)),
        (None, Some(k)) => blocks.push((0, k)),
        _ => (),
    }

    for (j, k) in blocks.drain(..) {
        assert!(j != k);
        let sz = if k >= j {
            res[k].x + 1 - res[j].x
        } else {
            let lr = res.last().unwrap();
            lr.x + 2 - res[j].x + res[k].x
        };
        let mut recalc = None;
        if sz > vc.cfg.small_deletion_limit() {
            trace!(
                "Large deletion block: {} - {}, len = {}",
                res[j].x + 1,
                res[k].x + 1,
                sz
            );
            let mut large_del = get_large_deletions(
                res,
                j,
                k,
                ref_len,
                pw.depth.dalign.as_mut().unwrap(),
                vc.cfg,
            );
            // Recalculate allele counts under any large deletions
            if !large_del.is_empty() {
                recalc_allele_counts(pw, &large_del);
                let a = large_del[0].start % l;
                let b = large_del.iter().map(|d| d.end()).max().unwrap() % l;
                recalc = Some((a, b));
                for d in large_del.drain(..) {
                    del_vec.push(d)
                }
            }
        }
        if let Some((a, b)) = recalc {
            let mut nres = Vec::with_capacity(sz);
            if b >= a {
                estimate_single_base_freq(pw, &mut nres, vc, Some(a), Some(b + 1));
                for (v1, v2) in res[a..=b].iter_mut().zip(nres.drain(..)) {
                    *v1 = v2
                }
            } else {
                estimate_single_base_freq(pw, &mut nres, vc, Some(a), None);
                for (v1, v2) in res[a..].iter_mut().zip(nres.drain(..)) {
                    *v1 = v2
                }
                estimate_single_base_freq(pw, &mut nres, vc, None, Some(b + 1));
                for (v1, v2) in res[..=b].iter_mut().zip(nres.drain(..)) {
                    *v1 = v2
                }
            }
        }
    }
    del_vec
}
