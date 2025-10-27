use std::{
    fs::File,
    io::{BufWriter, Write},
};

use rstar::{AABB, RTree, RTreeObject};

use super::{DelSet, DelType, Deletion, Deletions, ReadDels, ReadExtent, del_set::calc_n_u64};

use crate::cli::Guides;

/// We assume that a DEL can not be observed for technical reasons if it occurs within
/// DEL_LIMIT of the end of a read
const DEL_LIMIT: isize = 20;

impl RTreeObject for Deletion {
    type Envelope = AABB<[isize; 2]>;
    fn envelope(&self) -> Self::Envelope {
        AABB::from_corners([self.coords.start, 0], [self.coords.end, 1])
    }
}

impl<'a> Deletions<'a> {
    fn build_rtree(&mut self) -> Option<(RTree<Deletion>, usize)> {
        debug!("Building RTree from deletions");
        self.del_hash.take().map(|dh| {
            let del_vec: Vec<_> = dh.into_values().collect();
            let n_dels = del_vec.len();
            assert!(n_dels > 0);
            (RTree::bulk_load(del_vec), n_dels)
        })
    }

    pub(super) fn get_like_data(&mut self) -> Option<(Vec<(u64, u64)>, ReadDels, RTree<Deletion>)> {
        info!("Collecting information to estimate deletion frequencies");
        if let Some((rt, n_dels)) = self.build_rtree()
            && let Some(re) = self.read_extents.take()
        {
            let (obs, rd) = self.collect_contrib(n_dels, &rt, &re);

            Some((obs, rd, rt))
        } else {
            None
        }
    }

    /// Collect contributions to likelihood for the iterative steps
    ///
    /// Observed deletions contributes counts for each deletion while reads with no observed deletions
    /// contribute counts to all deletions that are not covered by the reads
    fn collect_contrib(
        &self,
        n_dels: usize,
        rt: &RTree<Deletion>,
        re: &[ReadExtent],
    ) -> (Vec<(u64, u64)>, ReadDels) {
        debug!("Collecting contributions from reads");
        // Vector of (observations, reads_covering) for each observed deletion.
        // obs_counts[n_dels] is for the no deletion (i.e. wildtype) case
        let mut obs_counts = vec![(0u64, 0u64); n_dels + 1];

        // Index for wildtype allele
        let wildtype = n_dels;

        // Add a dummy observation for the wildtype allele so that the mle whould be >0
        obs_counts[wildtype] = (0, 1);
        let mut read_dels = ReadDels::new(n_dels);
        let mask = read_dels.del_mask();

        for r in re.iter() {
            if let Some(i) = r.obs_deletion() {
                obs_counts[i as usize].0 += 1;
            }

            // Build up bitmaps for observed, excluded and non-observed deletions
            let mut ds = DelSet::new(n_dels, 3);

            let k = calc_n_u64(n_dels);
            if !get_contained_dels(rt, r, &mut obs_counts, &mut ds) {
                continue;
            }

            get_excluded_dels(rt, r, &mut ds);

            // Remove excluded dels from contained set

            let (d0, d1) = ds.all_dset().split_at_mut(k);
            for (m0, m1) in d0.iter_mut().zip(d1.iter()) {
                // Turn off bits (dels) in the observed set that are also in the excluded set
                *m0 ^= *m0 & *m1
            }

            // Set non-observed set as the complement of the union of the contained and excluded sets
            let (d0, d2) = ds.all_dset().split_at_mut(2 * k);
            let msk = mask.dset(0);
            for (i, m2) in d2.iter_mut().enumerate() {
                *m2 = msk[i] ^ (d0[i] | d0[i + k]);
            }

            /* 
            let msg = ["Contained", "Excluded", "Not observed"];

            if let Some(k) = r.guide.map(|g| g.index())
                && k == 8
            {
                for (ix, s) in msg.iter().enumerate() {
                    eprintln!("OOOK! {r:?}\n\t{s}:");
                    let d0 = ds.dset(ix);
                    for (j, m) in d0.iter().enumerate() {
                        let mut x = *m;
                        let mut i = 0;
                        while x != 0 {
                            if (x & 1) == 1 {
                                let (s, e) = self.coord_vec[(j << 6) | i].coords();
                                eprintln!("\t\t{s}\t{e}");
                            }
                            x >>= 1;
                            i += 1;
                        }
                    }
                }
            } */
            for (((a, b), c), d) in ds
                .dset(0)
                .iter()
                .zip(ds.dset(1).iter())
                .zip(ds.dset(2).iter())
                .zip(msk.iter())
            {
                if (a & b) != 0 || (a & c) != 0 || (b & c) != 0 {
                    panic!("Overlap between deletion sets")
                }
                if (a | b | c) != *d {
                    panic!("Union of sets does not include all deletions")
                }
            }
            read_dels.add(ds, r.obs_deletion());
        }
        for (i, o) in obs_counts[..n_dels - 1].iter().enumerate() {
            assert!(o.0 > 0 && o.0 <= o.1, "Odd obs for deletion {i} {o:?}");
        }

        (obs_counts, read_dels)
    }
}

/// Get deletions totally contained within read
fn get_contained_dels(
    rt: &RTree<Deletion>,
    r: &ReadExtent,
    obs_counts: &mut [(u64, u64)],
    ds: &mut DelSet,
) -> bool {
    if let Some((bb, bb2)) = get_bb(r) {
        // Collect deletions that are completely covered by this read
        for d in rt.locate_in_envelope_intersecting(&bb) {
            let i = d.ix;
            obs_counts[i as usize].1 += 1;
            ds.add_index(i as usize, 0);
        }

        // If necessary, check in alternate spae of dual reference
        if let Some(bb) = bb2.as_ref() {
            for d in rt.locate_in_envelope_intersecting(bb) {
                let i = d.ix;
                obs_counts[i as usize].1 += 1;
                ds.add_index(i as usize, 0);
            }
        }
        true
    } else {
        false
    }
}

fn get_excluded_dels(rt: &RTree<Deletion>, r: &ReadExtent, ds: &mut DelSet) {
    let (bb1, bb2) = if let Some(g) = r.guide {
        g.get_bb(r.reversed(), r.target_size)
    } else {
        (
            get_bb_from_point(r.physical_start(), r.reversed()),
            get_bb_from_point(r.physical_start() + r.target_size, r.reversed()),
        )
    };

    // Find deletions that overlap
    for d in rt.locate_in_envelope_intersecting(&bb1) {
        let i = d.ix;
        ds.add_index(i as usize, 1);
    }
    // And check in 2nd copy of dual referenced
    for d in rt.locate_in_envelope_intersecting(&bb2) {
        let i = d.ix;
        ds.add_index(i as usize, 1);
    }
}

type BBox = AABB<[isize; 2]>;

fn get_bb(rd: &ReadExtent) -> Option<(BBox, Option<BBox>)> {
    let (s, e) = rd.coords();
    match (get_bb_from_coords(s, e), rd.coords2()) {
        (Some(bb), Some((s, e))) => Some((bb, get_bb_from_coords(s, e))),
        (Some(bb), None) => Some((bb, None)),
        _ => None,
    }
}

fn get_bb_from_coords(s: isize, e: isize) -> Option<AABB<[isize; 2]>> {
    let x = s + DEL_LIMIT;
    let y = e - DEL_LIMIT;
    if x < y {
        Some(AABB::from_corners([x, 0], [y, 1]))
    } else {
        None
    }
}

fn get_bb_from_point(s: isize, rev: bool) -> AABB<[isize; 2]> {
    if DEL_LIMIT < 2 {
        AABB::from_point([s, 0])
    } else {
        let (x, y) = if rev { (s, s + 50) } else { (s - 50, s) };
        AABB::from_corners([x, 0], [y, 1])
    }
}

pub fn write_deletions(
    rt: &RTree<Deletion>,
    prefix: &str,
    fq: &[f64],
    wt_ci: &[f64; 2],
    guides: Option<&Guides>,
) -> anyhow::Result<()> {
    let file_name = format!("{}_del.txt", prefix);
    let mut wrt = BufWriter::new(File::create(&file_name)?);

    let mut v: Vec<_> = rt.iter().collect();
    v.sort_unstable_by_key(|d| d.coords);

    writeln!(
        wrt,
        "Start\tEnd\tSize\tType\tIndex\tFwd\tRev\tFreq\tGuides(fwd)\tGuides(neg)"
    )?;
    for d in v.into_iter() {
        write!(wrt, "{d}\t{:.6}", fq[d.ix as usize])?;
        if let Some(g) = guides {
            d.write_guides(g, &mut wrt)?;
        }
        writeln!(wrt)?;
    }
    writeln!(
        wrt,
        "-\t-\t-\t{}\t-\t-\t-\t{:.6} ({:.6}-{:.6})",
        DelType::Wildtype,
        fq.last().unwrap(),
        wt_ci[0],
        wt_ci[1]
    )?;
    Ok(())
}
