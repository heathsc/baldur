use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufWriter, Write},
};

use log::Level::Trace;

use r_htslib::*;

use crate::{
    cli::Config,
    fisher::FisherTest,
    reference::RefPos,
    vcf::{write_vcf_header, VcfCalc, VcfRes},
    depth::*,
    context::*,
    alleles::*,
    read::read_file,
};

type Qhist = [[usize; 4]; 64];

pub(crate) struct ProcWork<'a> {
    pub(crate) ref_seq: &'a [RefPos],
    pub(crate) depth: Depth,
    pub(crate) qual_hist: Qhist,
    pub(crate) ctxt_hist: [Qhist; N_CTXT],
}

fn output_calibration_data(cfg: &Config, pw: &ProcWork) -> io::Result<()> {
    if cfg.output_qual_calib() {
        let calc_q = |q: &[usize]| {
            let p = (q[1] + 1) as f64 / ((q[0] + q[1] + 2) as f64);
            -10.0 * p.log10()
        };
        let qual_cal_output = format!("{}_qcal.txt", cfg.output_prefix());
        let mut wrt = BufWriter::new(File::create(&qual_cal_output)?);
        write!(wrt, "Qual\tEmp_Qual\tEmp_Qual_Del\tMatch\tMismatch\tMatch_Del\tMismatch_Del")?;
        for i in 0..N_CTXT {
            write!(wrt, "\t{:#}\t\t\t\t\t", Ctxt5((i as u16) << 2))?;
        }
        writeln!(wrt)?;
        for (q, qc) in pw.qual_hist.iter().enumerate() {
            if qc[0] + qc[1] + qc[2] + qc[3] > 0 {
                let empirical_q = calc_q(&qc[..2]);
                let empirical_qd = calc_q(&qc[2..]);
                write!(wrt, "{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}", q, empirical_q, empirical_qd,
                       qc[0], qc[1], qc[2], qc[3])?;
                for cc in pw.ctxt_hist.iter() {
                    let qc1 = &cc[q];
                    let empirical_q = calc_q(&qc1[..2]);
                    let empirical_qd = calc_q(&qc1[2..]);
                    write!(wrt, "\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}", empirical_q, empirical_qd,
                           qc1[0], qc1[1], qc1[2], qc1[3])?;
                }
                writeln!(wrt)?;
            }
        }
    }
    Ok(())
}

pub fn process_data(mut sam_file: HtsFile, mut sam_hdr: SamHeader, cfg: Config) -> io::Result<()> {

    let reg = cfg.region();
    let ref_seq = cfg.reference().contig(reg.tid()).unwrap().seq();
    let mut pw = ProcWork{
        depth: Depth::new(reg.len(), !cfg.no_call()),
        qual_hist: [[0; 4]; 64],
        ctxt_hist: [[[0; 4]; 64]; N_CTXT],
        ref_seq,
    };

    // Read in input file, collect information, write out view file
    read_file(&mut sam_file, &mut sam_hdr, &cfg, &mut pw)?;

    // Output calibration data if requested
    output_calibration_data(&cfg, &pw)?;

    let depth_output = format!("{}_depth.txt", cfg.output_prefix());
    let mut wrt = BufWriter::new(File::create(&depth_output)?);
    let mut vcf_wrt = if !cfg.no_call() {
        Some(write_vcf_header(&sam_hdr, &cfg)?)
    } else {
        None
    };

    let fisher_test = FisherTest::new();
    let mut vr_cache: Option<(usize, VcfRes, &mut DepthCounts, Option<Vec<AllDesc>>)> = None;
    let all_desc: Vec<AllDesc> = (0..5).map(|i| if i < 4 { AllDesc(vec!(i)) } else { AllDesc(vec!()) }).collect();

    let vcf_calc = VcfCalc {
        ftest: &fisher_test,
        ref_seq,
        homopolymer_limit: cfg.homopolymer_limit(),
        seq_len: ref_seq.len(),
        cfg: &cfg,
    };

    for (i, dep) in pw.depth.counts.iter_mut().enumerate() {
        let x = i + reg.start();
        let ref_base = BASES.as_bytes()[(ref_seq[x].base()) as usize] as char;

        // Write out depth record
        writeln!(wrt, "{}\t{}\t{}\t{}", sam_hdr.tid2name(reg.tid()), x + 1, ref_base, dep)?;

        // Handle VCF output
        if let Some(vwrt) = vcf_wrt.as_mut() {
            trace!("At: {}", x + 1);
            let (mut vr, all_desc1) = if let Some(ins_all) = pw.depth.ins_hash.get(&x) {

                // Check if there is at least 1 insertion with the required number of observatations
                if ins_all.hash.iter().map(|(_, &i)| &ins_all.alleles[i]).any(|ct| ct.cts[0] > 1 && ct.cts[1] > 1) {

                    // If so, then collect regular and insertion alleles
                    let mut adesc: Vec<AllDesc> = Vec::with_capacity(8);
                    let ref_ix = ref_seq[x].base();
                    // Add regular bases + del
                    adesc.push(AllDesc(vec!(ref_ix)));
                    for b in (0..5).filter(|&x| x != ref_ix) {
                        adesc.push(AllDesc(vec!(b)));
                    }
                    // Add possible insertion alleles
                    for (all, ct) in ins_all.hash.iter().map(|(all, &i)| (all, &ins_all.alleles[i])) {
                        if ct.cts[0] > 1 && ct.cts[1] > 1 {
                            adesc.push(AllDesc(all.to_vec()))
                        }
                    }
                    let (new_cts, new_qcts, indel_flag) =
                       get_allele_counts(&adesc, x, x, ref_seq.len(), pw.depth.dalign.as_ref().unwrap(), &pw.depth.ins_hash);
                    for a in adesc.iter_mut() {
                        if a.len() == 1 && a[0] == 4 { a.0.clear() }
                    }
                    if log_enabled!(Trace) {
                        trace!("Insertion counts: {}", x + 1);
                        for (a, c) in adesc.iter().zip(new_cts.iter()) {
                            trace!("{}\t{}\t{}", a, c[0], c[1]);
                        }
                    }
                    let vr = vcf_calc.get_mallele_freqs(&new_cts, &new_qcts, &indel_flag);
                    if log_enabled!(Trace) {
                        trace!("Insertion freq. estimates: {}", x + 1);
                        for ar in vr.alleles.iter() {
                            trace!("{}\t{}\t{}", &adesc[ar.ix], ar.res.freq, ar.res.lr_test);
                        }
                    }
                    dep.cts = new_cts;
                    dep.qcts = new_qcts;
                    (vr, Some(adesc))
                } else {

                    // No insertions, so just handle the regular alleles
                    if log_enabled!(Trace) {
                        trace!("Standard counts: {}", x + 1);
                        for (a, c) in all_desc[..5].iter().zip(dep.cts.iter()) {
                            trace!("{}\t{}\t{}", a, c[0], c[1]);
                        }
                    }
                    let (vr1, ad) = (vcf_calc.get_allele_freqs(x, &dep.cts[..5], &dep.qcts[..5]), None);
                    if log_enabled!(Trace) {
                        trace!("Standard freq. estimates: {}", x + 1);
                        for ar in vr1.alleles.iter() {
                            trace!("{}\t{}\t{}", &all_desc[ar.ix], ar.res.freq, ar.res.lr_test);
                        }
                    }
                    (vr1, ad)
                }
            } else {
                if log_enabled!(Trace) {
                    trace!("Standard counts: {}", x + 1);
                    for (a, c) in all_desc[..5].iter().zip(dep.cts.iter()) {
                        trace!("{}\t{}\t{}", a, c[0], c[1]);
                    }
                }
                let (vr1, ad) = (vcf_calc.get_allele_freqs(x, &dep.cts[..5], &dep.qcts[..5]), None);
                if log_enabled!(Trace) {
                    trace!("Standard freq. estimates: {}", x + 1);
                    for ar in vr1.alleles.iter() {
                        trace!("{}\t{}\t{}", &all_desc[ar.ix], ar.res.freq, ar.res.lr_test);
                    }
                }
                (vr1, ad)
            };

            // For the current locus, the deletion allele will always have the ix field set to 4, so
            // we can check if any of the retained alleles is a deletion by checking if any have ix == 4
            if !vr.alleles.iter().any(|ar| ar.ix == 4) {
                // None of the retained alleles is a deletion, so we can output the previous locus
                if let Some((x1, mut vr1, dep1, odesc)) = vr_cache.take() {
                    let desc = odesc.as_ref().unwrap_or(&all_desc);
                    if let Some(s) = vcf_calc.output(x1, &mut vr1, &dep1.cts, &dep1.qcts, desc) {
                        let rs = cfg.rs(x1 + 1).unwrap_or(".");
                        writeln!(vwrt, "{}\t{}\t{}\t{}", sam_hdr.tid2name(reg.tid()), x1 + 1, rs, s)?;
                    }
                }
                // Store current locus in cache
                vr_cache = Some((x, vr, dep, all_desc1))
            } else if let Some((x1, mut vr1, mut dep1, odesc)) = vr_cache.take() {
                // Current locus contains deletion allele, so we need to combine with the previous locus
                trace!("Merging {} to {}", x1 + 1, x + 1);
                let desc = odesc.as_ref().unwrap_or(&all_desc);
                let mut adesc: Vec<AllDesc> = Vec::new();

                // Make list of possible new alleles
                let mut hset: HashSet<Vec<u8>> = HashSet::new();
                let ad = all_desc1.as_ref().unwrap_or(&all_desc);
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
                            let new_ds = AllDesc(new_ds);
                            adesc.push(new_ds);
                            orig_all_list.push((all.ix, all1.ix));
                        } else {
                            trim_possible = false;
                        }
                    }
                }
                let (new_cts, new_qcts, indel_flag) = get_allele_counts(&adesc, x1, x, ref_seq.len(), pw.depth.dalign.as_ref().unwrap(), &pw.depth.ins_hash);
                if log_enabled!(Trace) {
                    trace!("Deletion counts : {}", x + 1);
                    for (a, c) in adesc.iter().zip(new_cts.iter()) {
                        trace!("{}\t{}\t{}", a, c[0], c[1]);
                    }
                }

                // Get allele frequency estimates
                let vr2 = vcf_calc.get_mallele_freqs(&new_cts, &new_qcts, &indel_flag);

                // Check whether we can split up the new locus
                let first_all = &adesc[vr2.alleles[0].ix];
                let mut indel = false;
                let mut n_orig_alleles = 0;
                for ar in vr2.alleles.iter() {
                    if ar.ix < vr1.alleles.len() { n_orig_alleles += 1 }
                    if adesc[ar.ix].len() != first_all.len() { indel = true }
                }

                if log_enabled!(Trace) {
                    trace!("Deletion freq. estimates : {}.  No alleles: {}, No orig alleles: {}, indel: {}", x + 1, vr2.alleles.len(), n_orig_alleles, indel);
                    for ar in vr2.alleles.iter() {
                        trace!("{}\t{}\t{}", &adesc[ar.ix], ar.res.freq, ar.res.lr_test);
                    }
                }

                let no_change = vr2.alleles.len() == n_orig_alleles && vr1.alleles.len() == n_orig_alleles;

                // See if deletion allele from allele 2 still present in retained joint allele
                if vr2.alleles.len() == 1 || no_change || (trim_possible && vr2.alleles.iter().all(|ar| orig_all_list[ar.ix].1 != 4)) {
                    trace!("Deletion alleles not retained.  Print out previous cluster and start a new one");
                    trace!("Removing eliminated alleles from current cluster");
                    let mut allele_flag = vec!(false; ad.len());
                    for ar in vr2.alleles.iter() {
                        let ix = orig_all_list[ar.ix].1;
                        allele_flag[ix] = true;
                    }
                    let mut alleles: Vec<_> = vr.alleles.drain(..).filter(|ar| allele_flag[ar.ix]).collect();
                    let z = alleles.iter().fold(0.0, |s, ar| s + ar.res.freq );
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
                    if let Some(s) = vcf_calc.output(x1, &mut vr1, &dep1.cts, &dep1.qcts, desc) {
                        let rs = cfg.rs(x1 + 1).unwrap_or(".");
                        writeln!(vwrt, "{}\t{}\t{}\t{}", sam_hdr.tid2name(reg.tid()), x1 + 1, rs, s)?;
                    }
                    // Add current locus to cache
                    vr_cache = Some((x, vr, dep, all_desc1))
                } else {
                    // Can not simplify, so add merged locus to cache
                    dep1.cts = new_cts;
                    dep1.qcts = new_qcts;
                    vr_cache = Some((x1, vr2, dep1, Some(adesc)))
                }
            } else {
                // This occurs if we have a deletion in the current locus but there is no subsequent locus
                // which can only happen if there is a deletion at the start of the region
                panic!("Deletions across origin not yet handled")
            }
        }
    }
    // Print out remaining (cached) locus if required

    if let Some(vwrt) = vcf_wrt.as_mut() {

        if let Some((x1, mut vr1, dep1, odesc)) = vr_cache.take() {
            let desc = odesc.as_ref().unwrap_or(&all_desc);
            if let Some(s) = vcf_calc.output(x1, &mut vr1, &dep1.cts, &dep1.qcts, desc) {
                let rs = cfg.rs(x1 + 1).unwrap_or(".");
                writeln!(vwrt, "{}\t{}\t{}\t{}", sam_hdr.tid2name(reg.tid()), x1 + 1, rs, s)?;
            }
        }
    }
    Ok(())
}
