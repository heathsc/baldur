use std::{
    collections::{HashMap, HashSet},
    io::BufRead,
    path::{Path, PathBuf},
};

use clap::{command, value_parser, Arg, ArgAction};
use compress_io::compress::CompressIo;
use r_htslib::{Hts, HtsHdr};
use regex::{Match, Regex};

use crate::{
    log_utils::LogLevel,
    model::{setup_qual_model, N_QUAL},
    reference::Reference,
};

pub struct Region {
    tid: usize,
    start: usize,
    stop: usize,
    ctg_size: usize,
}

fn parse_usize_with_commas(s: &str) -> Option<usize> {
    s.replace(',', "").parse::<usize>().ok()
}

impl Region {
    pub fn from_str(reg_str: &str, reference: &Reference) -> anyhow::Result<Self> {
        let err = |s| Err(anyhow!("Could not parse region string '{}'", s,));

        let parse_x = |s: Match| parse_usize_with_commas(s.as_str());

        let reg = Regex::new(r#"^([^:]+):?([0-9,]+)?-?([0-9,]+)?"#).unwrap();
        if let Some(cap) = reg.captures(reg_str) {
            match (cap.get(1), cap.get(2), cap.get(3)) {
                (Some(c), None, None) => Self::new(c.as_str(), None, None, reference),
                (Some(c), Some(p), None) => Self::new(c.as_str(), parse_x(p), None, reference),
                (Some(c), Some(p), Some(q)) => {
                    Self::new(c.as_str(), parse_x(p), parse_x(q), reference)
                }
                _ => err(reg_str),
            }
        } else {
            err(reg_str)
        }
    }

    pub fn new(
        chrom: &str,
        start: Option<usize>,
        stop: Option<usize>,
        reference: &Reference,
    ) -> anyhow::Result<Self> {
        let ctg = reference
            .name2contig(chrom)
            .ok_or_else(|| anyhow!("Contig {} not present in input file", chrom))?;
        let start = start.unwrap_or(1).max(1) - 1;
        let ctg_size = ctg.size();
        assert!(ctg_size > 0, "Zero sized contig!");
        let stop = stop.unwrap_or(ctg_size).max(1).min(ctg_size) - 1;
        debug!("Chromosome region: {}:{}-{}", ctg.name(), start, stop);
        if stop >= start {
            Ok(Region {
                tid: ctg.tid(),
                start,
                stop,
                ctg_size,
            })
        } else {
            Err(anyhow!("Invalid range - stop < start".to_string(),))
        }
    }

    pub fn tid(&self) -> usize {
        self.tid
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn len(&self) -> usize {
        self.stop + 1 - self.start
    }

    pub fn ctg_size(&self) -> usize {
        self.ctg_size
    }
}

#[derive(Debug, Copy, Clone)]
pub enum ThresholdType {
    Hard,
    Soft,
}

impl ThresholdType {
    fn idx(&self) -> usize {
        match self {
            Self::Hard => 0,
            Self::Soft => 1,
        }
    }
}

pub struct Config {
    region: Region,
    reference: Reference,
    output_prefix: Box<str>,
    sample: Option<Box<str>>,
    adjust: usize,
    small_deletion_limit: usize,
    blacklist: Option<HashSet<usize>>,
    qual_calib: Option<[[[u8; 2]; N_QUAL]; N_QUAL]>,
    rs: Option<HashMap<usize, Box<str>>>,
    snv_thresholds: [f64; 2],
    indel_thresholds: [f64; 2],
    qual_table: [f64; N_QUAL],
    mapq_threshold: u8,
    qual_threshold: u8,
    max_qual: u8,
    max_indel_qual: u8,
    homopolymer_limit: u8,
    paired_end: bool,
    no_call: bool,
    rejected: bool,
    output_qual_calib: bool,
    view: bool,
    output_deletions: bool,
}

impl Config {
    pub fn region(&self) -> &Region {
        &self.region
    }

    pub fn mapq_threshold(&self) -> u8 {
        self.mapq_threshold
    }

    pub fn qual_threshold(&self) -> u8 {
        self.qual_threshold
    }

    pub fn snv_threshold(&self, t: ThresholdType) -> f64 {
        self.snv_thresholds[t.idx()]
    }

    pub fn indel_threshold(&self, t: ThresholdType) -> f64 {
        self.indel_thresholds[t.idx()]
    }

    pub fn max_qual(&self) -> u8 {
        self.max_qual
    }

    pub fn max_indel_qual(&self) -> u8 {
        self.max_indel_qual
    }

    pub fn homopolymer_limit(&self) -> u8 {
        self.homopolymer_limit
    }

    pub fn paired_end(&self) -> bool {
        self.paired_end
    }

    pub fn no_call(&self) -> bool {
        self.no_call
    }

    pub fn rejected(&self) -> bool {
        self.rejected
    }
    
    pub fn view(&self) -> bool {
        self.view
    }

    pub fn output_deletions(&self) -> bool {
        self.output_deletions
    }

    pub fn have_qual_calib(&self) -> bool {
        self.qual_calib.is_some()
    }

    pub fn output_qual_calib(&self) -> bool {
        self.output_qual_calib
    }

    pub fn output_prefix(&self) -> &str {
        &self.output_prefix
    }

    pub fn sample(&self) -> Option<&str> {
        self.sample.as_ref().map(|x| x as &str)
    }

    pub fn adjust(&self) -> usize {
        self.adjust
    }

    pub fn small_deletion_limit(&self) -> usize {
        self.small_deletion_limit
    }

    pub fn reference(&self) -> &Reference {
        &self.reference
    }

    pub fn blacklist(&self, x: usize) -> bool {
        match &self.blacklist {
            Some(m) => m.contains(&x),
            None => false,
        }
    }

    pub fn qual_calib(&self, ctxt: usize, q: u8) -> Option<&[u8; 2]> {
        self.qual_calib.as_ref().map(|v| {
            &v[ctxt as usize][if q != 61 {
                q.min(self.max_qual) as usize
            } else {
                61
            }]
        })
    }

    pub fn qual_table(&self) -> &[f64; N_QUAL] {
        &self.qual_table
    }

    pub fn rs(&self, x: usize) -> Option<&str> {
        self.rs.as_ref().and_then(|h| h.get(&x).map(|s| s as &str))
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
            if let (Ok(i1), Ok(j1)) = (i, j) {
                if j1 > i1 {
                    cts.1 += j1 - i1;
                    for ix in i1..j1 {
                        hs.insert(ix + 1);
                    }
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
            if let (Ok(i1), Ok(j1)) = (i, j) {
                if j1 == i1 + 1 {
                    hs.insert(j1, Box::from(fields[3]));
                    cts += 1;
                }
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

pub fn handle_cli() -> anyhow::Result<(Hts, Config)> {
    let m = command!()
        .arg(
            Arg::new("mapq_threshold")
                .short('q')
                .long("mapq-threshold")
                .default_value("0")
                .value_parser(value_parser!(u8))
                .value_name("MAPQ")
                .help("MAPQ quality threshold"),
        )
        .arg(
            Arg::new("max_qual")
                .short('M')
                .long("max-qual")
                .default_value("30")
                .value_parser(value_parser!(u8).range(..N_QUAL as i64))
                .value_name("QUAL")
                .help("Quality values are capped at this value"),
        )
        .arg(
            Arg::new("max_indel_qual")
                .short('I')
                .long("max-indel-qual")
                .default_value("20")
                .value_parser(value_parser!(u8).range(..N_QUAL as i64))
                .value_name("QUAL")
                .help("Quality values for indels are capped at this value"),
        )
        .arg(
            Arg::new("snv_thresholds")
                .long("snv-thresholds")
                .num_args(2)
                .value_parser(value_parser!(f64))
                .value_delimiter(',')
                .default_values(&["0.0005", "0.0025"])
                .value_name("FREQ1, FREQ2")
                .help("SNV's with frequency below FREQ1 (hard limit) are not reported, and below FREQ2 (soft limit) have a low_freq warning"),
        )
        .arg(
            Arg::new("indel_thresholds")
                .long("indel-thresholds")
                .num_args(2)
                .value_delimiter(',')
                .value_parser(value_parser!(f64))
                .default_values(&["0.025", "0.1"])
                .value_name("FREQ1, FREQ2")
                .help("Indels with freq below FREQ1 (hard limit) are not reported, and below FREQ2 (soft limit) have a low_freq warning"),
        )
        .arg(
            Arg::new("qual_threshold")
                .short('Q')
                .long("qual-threshold")
                .default_value("0")
                .value_parser(value_parser!(u8).range(..N_QUAL as i64))
                .value_name("QUAL")
                .help("Base quality threshold"),
        )
        .arg(
            Arg::new("homopolymer_limit")
                .short('P')
                .long("homopolymer-limit")
                .default_value("4")
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .help("Minimum size of homopolymer runs flagged as problematic"),
        )
        .arg(
            Arg::new("qual_calib")
                .hide(true)
                .long("qual-calib")
                .value_name("QCAL file")
                .value_parser(value_parser!(PathBuf))
                .help("File with quality calibration estimates"),
        )
        .arg(
            Arg::new("small_deletion_limit")
                .long("small-deletion-limit")
                .value_parser(value_parser!(usize))
                .value_name("SIZE")
                .default_value("64")
                .help("Maximum size for a small (explicit) deletion"),
        )
        .next_help_heading("Operation")
        .arg(
            Arg::new("adjust")
                .short('a')
                .long("adjust")
                .value_parser(value_parser!(usize))
                .value_name("ADJUST")
                .default_value("0")
                .help("Adjustment to genomic position"),
        )
        .arg(
            Arg::new("region")
                .short('r')
                .long("region")
                .value_parser(value_parser!(String))
                .value_name("REGION")
                .help("Genomic region to consider"),
        )
        .next_help_heading("Input/Output")
        .arg(
            Arg::new("reference")
                .short('T')
                .long("reference")
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .value_name("FASTA File")
                .help("Reference FASTA file"),
        )
        .arg(
            Arg::new("rs_list")
                .long("rs-list")
                .value_parser(value_parser!(PathBuf))
                .value_name("BED file")
                .help("BED file with dbSNP rs identifiers"),
        )
        .arg(
            Arg::new("blacklist")
                .long("blacklist")
                .value_parser(value_parser!(PathBuf))
                .value_name("BED file")
                .help("BED file with list of blacklisted sites"),
        )
        .arg(
            Arg::new("loglevel")
                .short('l')
                .long("loglevel")
                .value_name("LOGLEVEL")
                .value_parser(value_parser!(LogLevel))
                .ignore_case(true)
                .default_value("info")
                .help("Set log level"),
        )
        .arg(
            Arg::new("output_prefix")
                .short('o')
                .long("output-prefix")
                .default_value("baldur")
                .value_parser(value_parser!(String))
                .value_name("PREFIX")
                .help("Output prefix"),
        )
        .arg(
            Arg::new("no_call")
                .long("no-call")
                .action(ArgAction::SetTrue)
                .help("Do not call variants"),
        )
        .arg(
            Arg::new("view")
                .short('V')
                .long("view")
                .action(ArgAction::SetTrue)
                .help("Generate pileup view"),
        )
        .arg(
            Arg::new("rejected")
                .long("rejected")
                .action(ArgAction::SetTrue)
                .help("Write out list of rejected read IDs"),
        )
        .arg(
            Arg::new("output_deletions")
                .short('D')
                .action(ArgAction::SetTrue)
                .long("output-deletions")
                .help("Output large deletions"),
        )
        .arg(
            Arg::new("output_qual_calib")
                .hide(true)
                .long("output-qual-calib")
                .action(ArgAction::SetTrue)
                .help("Output quality calibration values"),
        )
        .arg(
            Arg::new("sample")
                .short('n')
                .long("sample")
                .value_parser(value_parser!(String))
                .value_name("NAME")
                .help("Sample name (for VCF file)"),
        )
        .arg(
            Arg::new("input")
                .value_parser(value_parser!(String))
                .value_name("INPUT")
                .help("Input SAM/BAM/CRAM file(s)"),
        )
        .arg(
            Arg::new("paired_end")
                .long("paired-end")
                .hide(true)
                .action(ArgAction::SetTrue)
                .help("Illumina paired end reads"),
        )
        .get_matches();

    super::log_utils::init_log(&m);

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
        read_blacklist(file, hdr.tid2name(region.tid))
            .expect("Could not read in blacklist from file")
    });

    let qual_calib = m.get_one::<PathBuf>("qual_calib").map(|file| {
        read_qual_calib(file, max_qual).expect("Could not read in bias list from file")
    });

    let rs = m.get_one::<PathBuf>("rs_list").map(|file| {
        read_rslist(file, hdr.tid2name(region.tid)).expect("Could not read in rs list from file")
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
