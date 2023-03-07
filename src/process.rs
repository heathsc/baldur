use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufWriter, Write},
};

use log::Level::Trace;

use r_htslib::*;

use crate::{
    alleles::*,
    cli::Config,
    context::*,
    deletions::*,
    depth::*,
    fisher::FisherTest,
    freq::{estimate_single_base_freq, process_large_deletions},
    read::read_file,
    reference::RefPos,
    vcf::{write_vcf_header, VcfCalc, VcfRes},
};

type Qhist = [[usize; 4]; 64];

pub(crate) struct ProcWork<'a> {
    pub(crate) ref_seq: &'a [RefPos],
    pub(crate) depth: Depth,
    pub(crate) qual_hist: Qhist,
    pub(crate) ctxt_hist: [Qhist; N_CTXT],
    pub(crate) dels: Option<Deletions>,
}

fn output_calibration_data(cfg: &Config, pw: &ProcWork) -> io::Result<()> {
    if cfg.output_qual_calib() {
        let calc_q = |q: &[usize]| {
            let p = (q[1] + 1) as f64 / ((q[0] + q[1] + 2) as f64);
            -10.0 * p.log10()
        };
        let qual_cal_output = format!("{}_qcal.txt", cfg.output_prefix());
        let mut wrt = BufWriter::new(File::create(&qual_cal_output)?);
        write!(
            wrt,
            "Qual\tEmp_Qual\tEmp_Qual_Del\tMatch\tMismatch\tMatch_Del\tMismatch_Del"
        )?;
        for i in 0..N_CTXT {
            write!(wrt, "\t{:#}\t\t\t\t\t", Ctxt5((i as u16) << 2))?;
        }
        writeln!(wrt)?;
        for (q, qc) in pw.qual_hist.iter().enumerate() {
            if qc[0] + qc[1] + qc[2] + qc[3] > 0 {
                let empirical_q = calc_q(&qc[..2]);
                let empirical_qd = calc_q(&qc[2..]);
                write!(
                    wrt,
                    "{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}",
                    q, empirical_q, empirical_qd, qc[0], qc[1], qc[2], qc[3]
                )?;
                for cc in pw.ctxt_hist.iter() {
                    let qc1 = &cc[q];
                    let empirical_q = calc_q(&qc1[..2]);
                    let empirical_qd = calc_q(&qc1[2..]);
                    write!(
                        wrt,
                        "\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}",
                        empirical_q, empirical_qd, qc1[0], qc1[1], qc1[2], qc1[3]
                    )?;
                }
                writeln!(wrt)?;
            }
        }
    }
    Ok(())
}

pub fn process_data(mut hts_file: Hts, cfg: Config) -> io::Result<()> {
    let reg = cfg.region();
    let ref_seq = cfg.reference().contig(reg.tid()).unwrap().seq();
    let dels = if cfg.output_deletions() {
        Some(Deletions::new(
            cfg.region().ctg_size(),
            cfg.small_deletion_limit(),
        ))
    } else {
        None
    };

    let mut pw = ProcWork {
        depth: Depth::new(reg.len(), !cfg.no_call()),
        qual_hist: [[0; 4]; 64],
        ctxt_hist: [[[0; 4]; 64]; N_CTXT],
        ref_seq,
        dels,
    };

    let threadpool = HtsThreadPool::new(1);
    if let Some(th) = threadpool.as_ref() {
        hts_file.hts_file_mut().set_thread_pool(th)?
    }

    // Read in input file, collect information, write out view file
    read_file(&mut hts_file, &cfg, &mut pw)?;

    let hdr = hts_file.header().unwrap();

    // Output calibration data if requested
    output_calibration_data(&cfg, &pw)?;

    let sam_hdr = if let HtsHdr::Sam(h) = hdr {
        h
    } else {
        panic!("Incorrect header")
    };

    let depth_output = format!("{}_depth.txt", cfg.output_prefix());
    let mut wrt = BufWriter::new(File::create(&depth_output)?);
    let mut vcf_wrt = if !cfg.no_call() {
        Some(write_vcf_header(sam_hdr, &cfg)?)
    } else {
        None
    };

    let ref_base = |x: usize| BASES.as_bytes()[(ref_seq[x].base()) as usize] as char;

    // Print out depth records
    for (i, dep) in pw.depth.counts.iter().enumerate() {
        let x = i + reg.start();

        // Write out depth record
        writeln!(
            wrt,
            "{}\t{}\t{}\t{}",
            sam_hdr.tid2name(reg.tid()),
            x + 1,
            ref_base(x),
            dep
        )?;
    }

    if let Some(dels) = pw.dels.as_ref() {
        let view_output = format!("{}_del.txt", cfg.output_prefix());
        let mut wrt = BufWriter::new(File::create(&view_output)?);
        for (d, x) in dels.iter() {
            writeln!(wrt, "{}\t{}", d, x)?;
        }
    }

    let fisher_test = FisherTest::new();
    let mut vr_cache: Option<(&mut VcfRes, &DepthCounts)> = None;

    // If variant calling required, perform frequency estimation
    let mut vcf_data = if vcf_wrt.is_some() {
        let all_desc: Vec<AllDesc> = (0..5)
            .map(|i| AllDesc::make(vec![i], None, Trunc::No))
            .collect();
        let vc = VcfCalc {
            ftest: &fisher_test,
            ref_seq,
            homopolymer_limit: cfg.homopolymer_limit(),
            seq_len: ref_seq.len(),
            all_desc, // Allele descriptions for the standard 5 bases (A C G T Del)
            cfg: &cfg,
        };

        let mut res = Vec::with_capacity(pw.depth.counts.len());

        // Calculate Individual base frequency estimates
        estimate_single_base_freq(&mut pw, &mut res, &vc, None, None);

        let large_dels = process_large_deletions(&mut pw, &mut res, &vc);

        Some((vc, res, large_dels))
    } else {
        None
    };

    // Output VCF data
    debug!("Starting VCF output");
    if let Some(vwrt) = vcf_wrt.as_mut() {
        let (vc, mut res, large_dels) = vcf_data.take().unwrap();
        // Find first locus that does not have a retained deletion
        if let Some((ix, _)) = res
            .iter()
            .enumerate()
            .find(|(_, vr)| vr.alleles.iter().all(|a| a.ix != 4))
        {
            trace!("Starting output at {}", res[ix].x + 1);
            let counts = &mut pw.depth.counts;
            let mut del_iter = large_dels.iter().peekable();
            for vr in res[ix..].iter_mut() {
                let x = vr.x;
                if vr.alleles.iter().all(|ar| ar.ix != 4) {
                    // None of the retained alleles is a deletion, so we can output the previous locus
                    if let Some((vr1, dep1)) = vr_cache.take() {
                        if let Some(del) = del_iter.next_if(|d| d.start <= vr1.x) {
                            let s = vc.del_output(del);
                            writeln!(
                                vwrt,
                                "{}\t{}\t.\t{}\t{}",
                                sam_hdr.tid2name(reg.tid()),
                                del.start + 1,
                                ref_base(del.start),
                                s
                            )?;
                        }
                        if let Some(s) = vc.output(vr1, &dep1.cts, &dep1.qcts) {
                            let rs = cfg.rs(vr1.x + 1).unwrap_or(".");
                            writeln!(
                                vwrt,
                                "{}\t{}\t{}\t{}",
                                sam_hdr.tid2name(reg.tid()),
                                vr1.x + 1,
                                rs,
                                s
                            )?;
                        }
                    }
                    // Store current locus in cache
                    vr_cache = Some((vr, &counts[x]))
                } else if let Some((mut vr1, dep1)) = vr_cache.take() {
                    trace!("Merging {} to {}", x + 1, vr1.x + 1);
                    let desc = vr1.adesc.as_ref().unwrap_or(&vc.all_desc);
                    let mut adesc: Vec<AllDesc> = Vec::new();

                    // Make list of possible new alleles
                    let mut hset: HashSet<Vec<u8>> = HashSet::new();
                    let ad = vr.adesc.as_ref().unwrap_or(&vc.all_desc);
                    let mut orig_all_list = Vec::with_capacity(ad.len() * desc.len());
                    let mut trim_possible = true;
                    for all1 in vr.alleles.iter() {
                        let ds1 = &ad[all1.ix];
                        for all in vr1.alleles.iter() {
                            let ds = &desc[all.ix];
                            let mut new_ds = ds.to_vec();
                            new_ds.extend_from_slice(ds1);
                            if !hset.contains(&new_ds) {
                                hset.insert(new_ds.to_vec());
                                let new_ds = AllDesc::make(new_ds, None, Trunc::No);
                                adesc.push(new_ds);
                                orig_all_list.push((all.ix, all1.ix));
                            } else {
                                trim_possible = false;
                            }
                        }
                    }

                    let (new_cts, new_qcts, indel_flag) = get_allele_counts(
                        &adesc,
                        vr1.x,
                        x,
                        ref_seq.len(),
                        pw.depth.dalign.as_ref().unwrap(),
                        &pw.depth.ins_hash,
                    );
                    if log_enabled!(Trace) {
                        trace!("Deletion counts : {}", x + 1);
                        for (a, c) in adesc.iter().zip(new_cts.iter()) {
                            trace!("{}\t{}\t{}", a, c[0], c[1]);
                        }
                    }

                    // Get allele frequency estimates
                    let vr2 = vc.get_mallele_freqs(vr1.x, &new_cts, &new_qcts, &indel_flag);

                    // Check whether we can split up the new locus
                    let first_all = &adesc[vr2.alleles[0].ix];
                    let mut indel = false;
                    let mut n_orig_alleles = 0;
                    for ar in vr2.alleles.iter() {
                        if ar.ix < vr1.alleles.len() {
                            n_orig_alleles += 1
                        }
                        if adesc[ar.ix].len() != first_all.len() {
                            indel = true
                        }
                    }

                    if log_enabled!(Trace) {
                        trace!("Deletion freq. estimates : {}.  No alleles: {}, No orig alleles: {}, indel: {}", x + 1, vr2.alleles.len(), n_orig_alleles, indel);
                        for ar in vr2.alleles.iter() {
                            trace!("{}\t{}\t{}", &adesc[ar.ix], ar.res.freq, ar.res.lr_test);
                        }
                    }

                    let no_change =
                        vr2.alleles.len() == n_orig_alleles && vr1.alleles.len() == n_orig_alleles;

                    // See if deletion allele from allele 2 still present in retained joint allele
                    if vr2.alleles.len() == 1
                        || no_change
                        || (trim_possible
                            && vr2.alleles.iter().all(|ar| orig_all_list[ar.ix].1 != 4))
                    {
                        trace!("Deletion alleles not retained.  Print out previous cluster and start a new one");
                        trace!("Removing eliminated alleles from current cluster");
                        let mut allele_flag = vec![false; ad.len()];
                        for ar in vr2.alleles.iter() {
                            let ix = orig_all_list[ar.ix].1;
                            allele_flag[ix] = true;
                        }
                        let mut alleles: Vec<_> = vr
                            .alleles
                            .drain(..)
                            .filter(|ar| allele_flag[ar.ix])
                            .collect();
                        let z = alleles.iter().fold(0.0, |s, ar| s + ar.res.freq);
                        if z < 1.0 {
                            for ar in alleles.iter_mut() {
                                ar.res.freq /= z;
                            }
                        }
                        vr.alleles = alleles;

                        if vr2.alleles.len() == 1 {
                            vr1.alleles.truncate(1);
                            vr1.alleles[0].res.freq = 1.0;
                        }

                        // Write out previous locus
                        if let Some(del) = del_iter.next_if(|d| d.start <= vr1.x) {
                            let s = vc.del_output(del);
                            writeln!(
                                vwrt,
                                "{}\t{}\t.\t{}\t{}",
                                sam_hdr.tid2name(reg.tid()),
                                del.start + 1,
                                ref_base(del.start),
                                s
                            )?;
                        }
                        if let Some(s) = vc.output(&mut vr1, &dep1.cts, &dep1.qcts) {
                            let rs = cfg.rs(vr1.x + 1).unwrap_or(".");
                            writeln!(
                                vwrt,
                                "{}\t{}\t{}\t{}",
                                sam_hdr.tid2name(reg.tid()),
                                vr1.x + 1,
                                rs,
                                s
                            )?;
                        }
                        // Add current locus to cache
                        vr_cache = Some((vr, &counts[x]))
                    } else {
                        // Can not simplify, so add merged locus to cache
                        counts[vr2.x].cts = new_cts;
                        counts[vr2.x].qcts = new_qcts;
                        vr.alleles = vr2.alleles;
                        vr.x = vr2.x;
                        vr.phred = vr2.phred;
                        vr.adesc = Some(adesc);
                        vr_cache = Some((vr, &counts[vr2.x]))
                    }
                } else {
                    panic!("Should not get here");
                }
            }
            if let Some((mut vr1, dep1)) = vr_cache.take() {
                if let Some(del) = del_iter.next_if(|d| d.start <= vr1.x) {
                    let s = vc.del_output(del);
                    writeln!(
                        vwrt,
                        "{}\t{}\t.\t{}\t{}",
                        sam_hdr.tid2name(reg.tid()),
                        del.start + 1,
                        ref_base(del.start),
                        s
                    )?;
                }
                if let Some(s) = vc.output(&mut vr1, &dep1.cts, &dep1.qcts) {
                    let rs = cfg.rs(vr1.x + 1).unwrap_or(".");
                    writeln!(
                        vwrt,
                        "{}\t{}\t{}\t{}",
                        sam_hdr.tid2name(reg.tid()),
                        vr1.x + 1,
                        rs,
                        s
                    )?;
                }
            }
        }
    }
    Ok(())
}
