use std::{
    collections::{HashMap, HashSet},
    io::BufRead,
    path::{Path, PathBuf},
};

use clap::ArgMatches;
use compress_io::compress::CompressIo;
use r_htslib::{Hts, HtsHdr};

use crate::{
    model::{N_QUAL, setup_qual_model},
    reference::Reference,
};

use super::{super::Region, Config};

impl Config {
    pub fn from_matches(m: &ArgMatches) -> anyhow::Result<(Hts, Self)> {
        let mapq_threshold = *m.try_get_one::<u8>("mapq_threshold")?.unwrap();
        let max_qual = *m.try_get_one::<u8>("max_qual")?.unwrap();
        let max_indel_qual = *m.try_get_one::<u8>("max_indel_qual")?.unwrap();
        let qual_threshold = *m
            .try_get_one::<u8>("qual_threshold")?
            .map(|x: &u8| x.min(&max_qual))
            .unwrap();
        let qual_table = setup_qual_model();

        let homopolymer_limit = *m.try_get_one::<u8>("homopolymer_limit")?.unwrap();
        let adjust = *m.try_get_one::<usize>("adjust")?.unwrap();
        let small_deletion_limit = *m.try_get_one::<usize>("small_deletion_limit")?.unwrap();
        let large_deletion_limit = *m.try_get_one::<usize>("large_deletion_limit")?.unwrap();
        let paired_end = m.get_flag("paired_end");
        let no_call = m.get_flag("no_call");
        let rejected = m.get_flag("rejected");
        let view = m.get_flag("view");
        let output_deletions = m.get_flag("output_deletions");
        let output_qual_calib = m.get_flag("output_qual_calib");

        let sample = m.get_one::<String>("sample").map(|s| Box::from(s.as_str()));

        let output_prefix = m
            .get_one::<String>("output_prefix")
            .map(|s| Box::from(s.as_str()))
            .unwrap();

        let input = m
            .get_one::<String>("input")
            .map(|s| s.as_str())
            .unwrap_or("-");

        let hts = Hts::open(Some(input), "r")?;

        let hdr = if let Some(HtsHdr::Sam(hdr)) = hts.header() {
            hdr
        } else {
            return Err(anyhow!(
                "Wrong input file type (expected SAM/BAM/CRAM)".to_string(),
            ));
        };

        debug!("Opened file {} for input", &input);

        let ref_file = m
            .get_one::<PathBuf>("reference")
            .expect("Missing reference"); // Should be enforced by clap

        let rdr = CompressIo::new().path(ref_file).bufreader()?;
        debug!("Opened {} for input", ref_file.display());
        let reference = Reference::from_reader(rdr, hdr)?;
        if reference.n_contigs() == 0 {
            return Err(anyhow!(
                "No contigs read in from reference file {}",
                ref_file.display()
            ));
        }
        debug!(
            "Reference read in successfully with {} contigs",
            reference.n_contigs()
        );

        let region = {
            if let Some(reg_str) = m.get_one::<String>("region") {
                let region = Region::from_str(reg_str, &reference)?;
                if reference.contig(region.tid()).is_none() {
                    return Err(anyhow!(
                        "Region {} not found in reference file {}",
                        reg_str,
                        ref_file.display()
                    ));
                }
                region
            } else {
                if reference.n_contigs() > 1 {
                    return Err(anyhow!(
                        "Region must be specified when reference file contains multiple contigs"
                            .to_string(),
                    ));
                }
                let tid = reference.contigs().keys().next().expect("Empty reference");
                Region::from_str(hdr.tid2name(*tid), &reference)?
            }
        };

        let blacklist = m.get_one::<PathBuf>("blacklist").map(|file| {
            read_blacklist(file, hdr.tid2name(region.tid()))
                .expect("Could not read in blacklist from file")
        });

        let qual_calib = m.get_one::<PathBuf>("qual_calib").map(|file| {
            read_qual_calib(file, max_qual).expect("Could not read in bias list from file")
        });

        let rs = m.get_one::<PathBuf>("rs_list").map(|file| {
            read_rslist(file, hdr.tid2name(region.tid()))
                .expect("Could not read in rs list from file")
        });

        let snv_thresholds = m
            .get_many::<f64>("snv_thresholds")
            .expect("Missing snv_thresholds")
            .copied()
            .collect::<Vec<_>>()
            .try_into()
            .expect("Unable to convert to array");

        let indel_thresholds = m
            .get_many::<f64>("indel_thresholds")
            .expect("Missing indel_thresholds")
            .copied()
            .collect::<Vec<_>>()
            .try_into()
            .expect("Unable to convert to array");

        Ok((
            hts,
            Config {
                output_prefix,
                sample,
                region,
                reference,
                adjust,
                small_deletion_limit,
                large_deletion_limit,
                blacklist,
                qual_calib,
                rs,
                mapq_threshold,
                qual_threshold,
                snv_thresholds,
                indel_thresholds,
                homopolymer_limit,
                qual_table,
                max_qual,
                max_indel_qual,
                paired_end,
                no_call,
                rejected,
                output_qual_calib,
                view,
                output_deletions,
            },
        ))
    }
}

fn read_blacklist<S: AsRef<Path>>(file: S, ctg: &str) -> anyhow::Result<HashSet<usize>> {
    let file = file.as_ref();
    let mut rdr = CompressIo::new().path(file).bufreader()?;
    let mut buf = String::new();
    let mut hs = HashSet::new();
    debug!("Reading in blacklist from {}", file.display());
    let mut cts = (0, 0);
    loop {
        if rdr.read_line(&mut buf)? == 0 {
            break;
        }
        let fields: Vec<_> = buf.trim().split('\t').collect();
        if fields[0] == ctg {
            cts.0 += 1;
            let i = fields[1].parse::<usize>();
            let j = fields[2].parse::<usize>();
            if let (Ok(i1), Ok(j1)) = (i, j)
                && j1 > i1
            {
                cts.1 += j1 - i1;
                for ix in i1..j1 {
                    hs.insert(ix + 1);
                }
            }
        }
        buf.clear();
    }
    debug!("Entries read: {}, sites: {}", cts.0, cts.1);
    Ok(hs)
}

fn read_rslist<S: AsRef<Path>>(file: S, ctg: &str) -> anyhow::Result<HashMap<usize, Box<str>>> {
    let file = file.as_ref();
    let mut rdr = CompressIo::new().path(file).bufreader()?;
    let mut buf = String::new();
    let mut hs = HashMap::new();
    debug!("Reading in rs list from {}", file.display());
    let mut cts = 0;
    loop {
        if rdr.read_line(&mut buf)? == 0 {
            break;
        }
        let fields: Vec<_> = buf.split('\t').collect();
        if fields[0] == ctg {
            let i = fields[1].parse::<usize>();
            let j = fields[2].parse::<usize>();
            if let (Ok(i1), Ok(j1)) = (i, j)
                && j1 == i1 + 1
            {
                hs.insert(j1, Box::from(fields[3]));
                cts += 1;
            }
        }
        buf.clear();
    }
    debug!("Entries read: {}", cts);
    Ok(hs)
}

fn read_qual_calib<S: AsRef<Path>>(
    file: S,
    max_qual: u8,
) -> anyhow::Result<[[[u8; 2]; N_QUAL]; N_QUAL]> {
    let file = file.as_ref();
    let mut rdr = CompressIo::new().path(file).bufreader()?;
    let mut buf = String::new();
    let mut cal = [[[0, 0]; N_QUAL]; N_QUAL];

    debug!("Reading in quality calibration from {}", file.display());
    let mut line = 0;
    loop {
        if rdr.read_line(&mut buf)? == 0 {
            break;
        }
        if line > 0 {
            let fields: Vec<_> = buf.trim_end().split('\t').collect();
            if fields.len() != 6 * 65 + 1 {
                return Err(anyhow!(
                    "Wrong number of fields at line {} of file {} (seen {}, expected {})",
                    line + 1,
                    file.display(),
                    fields.len(),
                    6 * 65 + 1
                ));
            }
            let q = fields[0].parse::<usize>().map_err(|e| {
                anyhow!(
                    "Could not parse quality score from column 1 at line {} of file {}: {}",
                    line + 1,
                    file.display(),
                    e
                )
            })?;

            for (ctxt, calp) in cal.iter_mut().enumerate() {
                let ix = ctxt * 6 + 7;
                let mut q1 = [0, 0];
                for k in 0..2 {
                    let z = fields[ix + k].parse::<f64>().map_err(|e| {
                        anyhow!(
                            "Could not parse float from column {} at line {} of file {}: {}",
                            ix + 1,
                            line + 1,
                            file.display(),
                            e
                        )
                    })?;
                    q1[k] = z.round().min(max_qual as f64) as u8;
                }
                calp[q] = q1;
            }
        }
        line += 1;
        buf.clear();
    }
    Ok(cal)
}
