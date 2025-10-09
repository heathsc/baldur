use std::{
    fs::File,
    io::{BufWriter, Write},
};

use rstar::{AABB, RTree, RTreeObject};

use super::{DelSet, DelType, Deletion, Deletions, LikeContrib, ReadDels, ReadExtent};

/// We assume that a DEL can not be observed for technical reasons if it occurs within
/// DEL_LIMIT of the end of a read
const DEL_LIMIT: isize = 20;

impl RTreeObject for Deletion {
    type Envelope = AABB<[isize; 2]>;
    fn envelope(&self) -> Self::Envelope {
        AABB::from_corners([self.coords.start, 0], [self.coords.end, 1])
    }
}

impl Deletions {
    fn build_rtree(&mut self) -> Option<(RTree<Deletion>, usize)> {
        debug!("Building RTree from deletions");
        self.del_hash.take().map(|dh| {
            let del_vec: Vec<_> = dh.into_values().collect();
            let n_dels = del_vec.len();
            assert!(n_dels > 0);
            (RTree::bulk_load(del_vec), n_dels)
        })
    }

    pub(super) fn get_like_data(&mut self) -> Option<(LikeContrib, RTree<Deletion>)> {
        info!("Collecting information to estimate deletion frequencies");
        if let Some((rt, n_dels)) = self.build_rtree()
            && let Some(re) = self.read_extents.take()
        {
            Some((collect_contrib(n_dels, &rt, &re), rt))
        } else {
            None
        }
    }
}

/// Collect contributions to likelihood for the iterative steps
///
/// Observed deletions contributes counts for each deletion while reads with no observed deletions
/// controbute counts to all deletions that are not covered by the reads
fn collect_contrib(n_dels: usize, rt: &RTree<Deletion>, re: &[ReadExtent]) -> LikeContrib {
    debug!("Collecting contributions from reads");
    // Vector of (observations, reads_covering) for each observed deletion.
    // obs_counts[n_dels] is for the no deletion (i.e. wildtype) case
    let mut obs_counts = vec![(0u64, 0u64); n_dels + 1];

    // Index for wildtype allele
    let wildtype = n_dels;

    // Add a dummy observation for the wildtype allele so that tje mle whould be >0
    obs_counts[wildtype] = (0, 1);
    let mut read_dels = ReadDels::new(n_dels);

    for (r, (bb, bb2)) in re.iter().filter_map(|r| get_bb(r).map(|bb| (r, bb))) {
        if let Some(i) = r.obs_deletion() {
            obs_counts[i as usize].0 += 1;
        }

        // Build up bitmaps for covered and excluded deletions
        let mut ds = DelSet::new(n_dels, 2);

        // Collect deletions that are completely covered by this read
        for d in rt.locate_in_envelope(&bb) {
            let i = d.ix;
            obs_counts[i as usize].1 += 1;
            ds.add_index(i as usize, 0);
        }

        // If necessary, check in alternate spae of dual reference
        if let Some(bb) = bb2.as_ref() {
            for d in rt.locate_in_envelope(bb) {
                let i = d.ix;
                obs_counts[i as usize].1 += 1;
                ds.add_index(i as usize, 0);
            }
        }

        // Find deletions that overlap physical start point of read
        let bb = get_bb_from_point(r.physical_start());
        for d in rt.locate_in_envelope_intersecting(&bb) {
            let i = d.ix;
            ds.add_index(i as usize, 1);
        }
        // And check in 2nd copy of dual referenced
        let bb = get_bb_from_point(r.physical_start() + r.target_size);
        for d in rt.locate_in_envelope_intersecting(&bb) {
            let i = d.ix;
            ds.add_index(i as usize, 1);
        }

        /* 
        let mut overlap = false;
        for (j, (a, b)) in ds.dset(0).iter().zip(ds.dset(1).iter()).enumerate() {
            let mut x = a & b;
            if x != 0 {
                if !overlap {
                    eprintln!("OOOK! {r:?}");
                    overlap = true;
                }

                let mut i = 0;
                while x != 0 {
                    if (x & 1) != 0 {
                        eprintln!("  Del {}", (j << 6) | i)
                    }
                    x >>= 1;
                    i += 1;
                }
            }
        } */

        assert!(
            ds.dset(0)
                .iter()
                .zip(ds.dset(1).iter())
                .all(|(a, b)| (a & b) == 0),
            "Covered and excluded sets overlap"
        ); 

        read_dels.add(ds, r.obs_deletion());
    }
    for (i, o) in obs_counts[..n_dels - 1].iter().enumerate() {
        assert!(o.0 > 0 && o.0 <= o.1, "Odd obs for deletion {i} {o:?}");
    }

    LikeContrib::new(obs_counts, read_dels)
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

fn get_bb_from_point(s: isize) -> AABB<[isize; 2]> {
    if DEL_LIMIT < 2 {
        AABB::from_point([s, 0])
    } else {
        let x = s - DEL_LIMIT + 2;
        let y = s + DEL_LIMIT - 2;
        AABB::from_corners([x, 0], [y, 1])
    }
}

pub fn write_deletions(rt: &RTree<Deletion>, prefix: &str, fq: &[f64]) -> anyhow::Result<()> {
    let file_name = format!("{}_del.txt", prefix);
    let mut wrt = BufWriter::new(File::create(&file_name)?);

    let mut v: Vec<_> = rt.iter().collect();
    v.sort_unstable_by_key(|d| d.coords);

    writeln!(wrt, "Start\tEnd\tSize\tType\tIndex\tFwd\tRev\tFreq")?;
    for d in v.into_iter() {
        writeln!(wrt, "{d}\t{}", fq[d.ix as usize])?;
    }
    writeln!(
        wrt,
        "-\t-\t-\t{}\t-\t-\t-\t{}",
        DelType::Wildtype,
        fq.last().unwrap()
    )?;
    Ok(())
}
