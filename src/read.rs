use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
};

use r_htslib::*;

use crate::{
    align_store::AlignStore, alleles::*, cli::Config, context::*, deletions::DelType,
    process::ProcWork,
};

const QUAL_SKIP: u8 = 20;
const MAX_SPLIT: usize = 8;

/// Process the BAM record(s) coming from a single read
fn handle_read<W: Write>(
    blst: &mut [BamRec],
    cfg: &Config,
    proc_work: &mut ProcWork<'_>,
    wrt: Option<&mut W>,
    wrt_rej: Option<&mut W>,
) -> anyhow::Result<()> {
    let reg = cfg.region();
    let qt = cfg.qual_threshold();
    let max_qual = cfg.max_qual();
    let max_indel_qual = cfg.max_indel_qual();
    let pe = cfg.paired_end();

    let fg = blst[0].flag() & BAM_FREVERSE;
    let calib_reqd = cfg.have_qual_calib() || cfg.output_qual_calib();

    // Sort BAM records by starting position within read
    if !pe {
        blst.sort_unstable_by_key(|b| {
            let cigar = b.cigar().expect("No Cigar!");
            let cig = if fg != 0 {
                cigar.iter().last()
            } else {
                cigar.iter().next()
            }
            .expect("Empty cigar");
            match cig.op() {
                CigarOp::SoftClip | CigarOp::HardClip => cig.op_len(),
                _ => 0,
            }
        });
    }

    // Only process read if all mappings are to the same strand of the same contig
    if blst
        .iter()
        .all(|rec| rec.tid() == Some(reg.tid()) && (pe || (rec.flag() & BAM_FREVERSE) == fg))
    {
        let mut prev = None;

        let mut align_store = AlignStore::new(cfg);
        let mut ins_allele = Vec::new();
        let mut ins_hash: HashMap<usize, (usize, u8)> = HashMap::new();

        let mut tot_mm = [0; 2];
        for b in blst.iter() {
            let reverse = (b.flag() & BAM_FREVERSE) != 0;
            if b.qual() >= cfg.mapq_threshold() {
                if let Some(ref_start) = b.pos() {
                    let cigar = b.cigar().expect("No Cigar!");

                    // Check that cigar does not have an Ins that does not follow a match.  We do not
                    // handle this (odd!) combination, so if we find it we will just skip the read
                    if cigar
                        .windows(2)
                        .any(|v| v[1].op() == CigarOp::Ins && v[0].op() != CigarOp::Match)
                    {
                        warn!("Bad CIGAR: Ins not preceded by a Match,  Skipping read");
                        continue;
                    }

                    if !(pe || reverse) {
                        if prev.is_some() {
                            if let Some(d) = proc_work.dels.as_mut() {
                                d.add_del(align_store.pos(), ref_start, false, DelType::Split)
                            }
                            align_store.fill_to(b'S', QUAL_SKIP, ref_start);
                        } else {
                            prev = Some(0)
                        }
                    }
                    let seq_qual = b.get_seq_qual()?;

                    // Evaluate context at each position in read
                    // (we could do this on the fly but it gets complicated with skips etc
                    // so I prefer to do it this way)

                    let read_len = seq_qual.len();
                    let mut context = Vec::with_capacity(read_len);
                    if calib_reqd {
                        let mut ctxt = Context5::new(reverse);
                        for (i, c) in seq_qual.iter().enumerate() {
                            ctxt.add_base(*c);
                            if i >= 2 {
                                context.push(ctxt)
                            }
                        }
                        ctxt = Context5::new(reverse);
                        context.push(ctxt);
                        context.push(ctxt);
                    } else {
                        let ctxt = Context5::new(reverse);
                        for _ in 0..read_len {
                            context.push(ctxt)
                        }
                    }
                    let mut it = seq_qual.iter().zip(context.iter()).peekable();

                    align_store.set_pos(ref_start);
                    let tr = if reverse {
                        [b'a', b'c', b'g', b't', b'b', b'd', b'h', b'u']
                    } else {
                        [b'A', b'C', b'G', b'T', b'B', b'D', b'H', b'U']
                    };
                    let (del, del_low, n) = if reverse {
                        (b'_', b'&', b'n')
                    } else {
                        (b'-', b'^', b'N')
                    };

                    let transform = |c: u8| match c >> 2 {
                        0 => n,
                        x if x < qt => tr[((c & 3) | 4) as usize],
                        _ => tr[(c & 3) as usize],
                    };

                    let chk_max = |c: &u8, mq: u8| {
                        let q = *c >> 2;
                        if q > mq && q != 61 {
                            (mq << 2) | (*c & 3)
                        } else {
                            *c
                        }
                    };

                    let qual_cal = |c: &u8, mq: u8, ix: usize, ctx: &Context5| match (
                        cfg.have_qual_calib(),
                        ctx.context3(),
                    ) {
                        (false, _) | (_, None) => chk_max(c, mq),
                        (true, Some(ct)) => cfg
                            .qual_calib(ct, *c >> 2)
                            .map(|x| (*c & 3) | (x[ix] << 2))
                            .unwrap_or_else(|| chk_max(c, mq)),
                    };

                    let mut prev_base: Option<(u8, Context5)> = None;
                    let mut cigar_it = cigar.iter().peekable();
                    while let Some(cig) = cigar_it.next() {
                        let l = cig.op_len() as usize;
                        assert!(l > 0, "Zero length cigar element");
                        let mut pbase = prev_base.take();
                        match cig.op() {
                            CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                                let mut v: Vec<(u8, u8, Context5)> = Vec::with_capacity(l);
                                let ins = if let Some(cig1) = cigar_it.peek() {
                                    cig1.op() == CigarOp::Ins
                                } else {
                                    false
                                };
                                for k in 0..l {
                                    let (c, ctxt) =
                                        it.next().expect("Mismatch between Cigar and sequence");
                                    let c1 = qual_cal(c, max_qual, 0, ctxt);
                                    if ins && k == l - 1 {
                                        prev_base = Some((c1, *ctxt));
                                    } else {
                                        let base = transform(c1);
                                        v.push((base, c1, *ctxt));
                                    }
                                }
                                let ct = align_store.add(&v, proc_work.ref_seq, proc_work);
                                tot_mm[0] += ct[0];
                                tot_mm[1] += ct[1];
                            }
                            CigarOp::SoftClip => {
                                it.nth(l - 1);
                            }
                            CigarOp::Ins => {
                                let (pbase, ctxt) =
                                    pbase.take().expect("Insertion not following a match!");
                                if let Some(pos) = align_store.in_region() {
                                    // Collect inserted allele
                                    ins_allele.clear();
                                    let mut min_q = pbase >> 2;
                                    ins_allele.push(pbase & 3);
                                    for _ in 0..l {
                                        let (c, ctxt) =
                                            it.next().expect("Mismatch between Cigar and sequence");
                                        let c1 = qual_cal(c, max_indel_qual, 0, ctxt);
                                        min_q = min_q.min(c1 >> 2);
                                        ins_allele.push(c1 & 3);
                                    }
                                    let mut base = if reverse { b'i' } else { b'I' };
                                    if min_q >= qt {
                                        let ins_all = proc_work
                                            .depth
                                            .ins_hash
                                            .entry(pos)
                                            .or_insert_with(InsAlleles::default);
                                        let ins_ix = ins_all.get_or_insert(&ins_allele);
                                        ins_all.add_obs(ins_ix, min_q, reverse);
                                        ins_hash.insert(pos, (ins_ix, min_q));
                                    } else {
                                        base += 1
                                    }
                                    // Adds previous base as an I (or i, J, j)
                                    align_store.add(
                                        &[(base, pbase, ctxt)],
                                        proc_work.ref_seq,
                                        proc_work,
                                    );
                                } else {
                                    it.nth(l - 1);
                                }
                            }
                            CigarOp::Del => {
                                let (c, ctxt) =
                                    it.peek().expect("Mismatch between Cigar and sequence");
                                let c1 = qual_cal(*c, max_indel_qual, 1, *ctxt);
                                let b = if (c1 >> 2) < qt { del_low } else { del };
                                let v = vec![(b, c1, **ctxt); l];
                                if let Some(d) = proc_work.dels.as_mut() {
                                    let x = align_store.pos();
                                    let y = x + l - 1;
                                    if reverse {
                                        d.add_del(y, x, true, DelType::Del)
                                    } else {
                                        d.add_del(x, y, false, DelType::Del)
                                    }
                                }
                                align_store.add(&v, proc_work.ref_seq, proc_work);
                            }
                            CigarOp::RefSkip => align_store.advance(l),
                            _ => (),
                        }
                    }
                    if !pe && reverse {
                        if let Some(x) = prev {
                            if let Some(d) = proc_work.dels.as_mut() {
                                d.add_del(x, align_store.pos(), true, DelType::Split)
                            }
                            align_store.fill_to(b's', QUAL_SKIP, x)
                        }
                        prev = Some(ref_start);
                    }
                }
            }
        }
        if align_store.changed() {
            let depth = &mut proc_work.depth;
            depth.add_obs_vec(
                align_store.seq(),
                align_store.qual(),
                ins_hash,
                blst[0].qname()?,
            );
            if let Some(w) = wrt {
                w.write_all(align_store.seq())?;
                if let Some(b) = blst.get(0) {
                    write!(w, "\t{}", b.qname()?)?;
                }
                writeln!(w)?;
            }
            if let Some(b) = blst.get(0) {
                let z = (tot_mm[1] as f64) / ((tot_mm[0] + tot_mm[1]) as f64);
                //            println!("{}\t{:.6}\t{}\t{}", b.qname()?, z, tot_mm[0], tot_mm[1]);
                if z >= 0.1 {
                    warn!(
                        "High mismatch rate for read {}\t{:.6}\t{}\t{}",
                        b.qname()?,
                        z,
                        tot_mm[0],
                        tot_mm[1]
                    );
                }
            }
        }
    } else if let Some(w) = wrt_rej {
        writeln!(w, "{}", blst[0].qname().expect("Missing query name"))?;
    }
    Ok(())
}

pub(crate) fn read_file(hts_file: &mut Hts, cfg: &Config, pw: &mut ProcWork) -> anyhow::Result<()> {
    let reg = cfg.region();

    let mut wrt = if cfg.view() {
        let view_output = format!("{}_view.txt", cfg.output_prefix());
        Some(BufWriter::new(File::create(&view_output)?))
    } else {
        None
    };

    let mut wrt_rej = if cfg.rejected() {
        let view_output = format!("{}_rejected.txt", cfg.output_prefix());
        Some(BufWriter::new(File::create(&view_output)?))
    } else {
        None
    };
    
    let mut blst = Vec::new();
    for _ in 0..MAX_SPLIT {
        blst.push(BamRec::new()?)
    }
    let mut b = BamRec::new()?;
    let mut idx = 0;

    while b.read(hts_file)? {
        if b.tid() == Some(reg.tid()) && b.qual() >= cfg.mapq_threshold() {
            if idx > 0 {
                if !b.qnames_eq(&blst[0])? {
                    if idx <= MAX_SPLIT {
                        handle_read(&mut blst[0..idx], cfg, pw, wrt.as_mut(), wrt_rej.as_mut())?
                    }
                    b.swap(&mut blst[0]);
                    idx = 1;
                } else {
                    if idx < MAX_SPLIT {
                        b.swap(&mut blst[idx]);
                    }
                    idx += 1;
                }
            } else {
                b.swap(&mut blst[0]);
                idx = 1;
            }
        }
    }

    // Process remaining reads (if any)
    if idx > 0 && idx <= MAX_SPLIT {
        handle_read(&mut blst[0..idx], cfg, pw, wrt.as_mut(), wrt_rej.as_mut())?
    }
    Ok(())
}
