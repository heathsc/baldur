use std::{
   collections::{HashMap, HashSet},
   io::{self, BufRead, Error, ErrorKind},
};

use lazy_static::lazy_static;
use regex::{Match, Regex};

use clap::{Command, Arg, crate_version};
use compress_io::compress::CompressIo;
use r_htslib::{HtsFile, SamHeader};

use super::io_err;
use crate::reference::Reference;
use crate::model::{setup_qual_model, N_QUAL};

lazy_static! {
    static ref RE_REGION: Regex = Regex::new(r#"^([^:]+):?([0-9,]+)?-?([0-9,]+)?"#).unwrap();
}

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
   pub fn from_str(reg_str: &str, reference: &Reference) -> io::Result<Self> {
      let err = |s| {
         Err(Error::new(
            ErrorKind::Other,
            format!("Could not parse region string '{}'", s),
         ))
      };

      let parse_x = |s: Match| parse_usize_with_commas(s.as_str());

      if let Some(cap) = RE_REGION.captures(reg_str) {
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
   ) -> io::Result<Self> {
      let ctg = reference.name2contig(chrom).ok_or_else(|| {
         Error::new(
            ErrorKind::Other,
            format!("Contig {} not present in input file", chrom),
         )
      })?;
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
         Err(Error::new(
            ErrorKind::Other,
            "Invalid range - stop < start".to_string(),
         ))
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
pub enum ThresholdType {Hard, Soft}

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
   blacklist: Option<HashSet<usize>>,
   qual_calib: Option<[[[u8; 2]; N_QUAL] ;N_QUAL]>,
   rs: Option<HashMap<usize, Box<str>>>,
   snv_thresholds: [f64; 2],
   indel_thresholds: [f64; 2],
   qual_table: [f64; N_QUAL],
   mapq_threshold: u8,
   qual_threshold: u8,
   max_qual: u8,
   homopolymer_limit: u8,
   paired_end: bool,
   no_call: bool,
   output_qual_calib: bool,
   view: bool,
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
      self.snv_thresholds[ t.idx() ]
   }

   pub fn indel_threshold(&self, t: ThresholdType) -> f64 {
      self.indel_thresholds[ t.idx() ]
   }

   pub fn max_qual(&self) -> u8 {
      self.max_qual
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

   pub fn view(&self) -> bool {
      self.view
   }

   pub fn have_qual_calib(&self) -> bool { self.qual_calib.is_some() }

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
         &v[ctxt as usize][if q != 61 { q.min(self.max_qual) as usize } else { 61 } ]
      })
   }

   pub fn qual_table(&self) -> &[f64; N_QUAL] { &self.qual_table }

   pub fn rs(&self, x: usize) -> Option<&str> {
      self.rs
         .as_ref()
         .and_then(|h| h.get(&x).map(|s| s as &str))
   }
}

fn read_blacklist(file: &str, ctg: &str) -> io::Result<HashSet<usize>> {
   let mut rdr = CompressIo::new().path(file).bufreader()?;
   let mut buf = String::new();
   let mut hs = HashSet::new();
   debug!("Reading in blacklist from {}", file);
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

fn read_rslist(file: &str, ctg: &str) -> io::Result<HashMap<usize, Box<str>>> {
   let mut rdr = CompressIo::new().path(file).bufreader()?;
   let mut buf = String::new();
   let mut hs = HashMap::new();
   debug!("Reading in rs list from {}", file);
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

fn read_qual_calib(file: &str, max_qual: u8) -> io::Result<[[[u8; 2]; N_QUAL]; N_QUAL]> {
   let mut rdr = CompressIo::new().path(file).bufreader()?;
   let mut buf = String::new();
   let mut cal =  [[[0, 0];N_QUAL];N_QUAL];

   debug!("Reading in quality calibration from {}", file);
   let mut line = 0;
   loop {
      if rdr.read_line(&mut buf)? == 0 {
         break;
      }
      if line > 0
      {
         let fields: Vec<_> = buf.trim_end().split('\t').collect();
         if fields.len() != 6 * 65 + 1 {
            return  Err(Error::new(
               ErrorKind::Other,
               format!("Wrong number of fields at line {} of file {} (seen {}, expected {})",
                       line + 1, file, fields.len(), 6 * 65 + 1),
            ))
         }
         let q = fields[0].parse::<usize>().map_err(|e|
            Error::new(
               ErrorKind::Other,
               format!("Could not parse quality score from column 1 at line {} of file {}: {}", line + 1, file, e
            )))?;
         for (ctxt, calp) in cal.iter_mut().enumerate() {
            let ix = ctxt * 6 + 7;
            let mut q1 = [0, 0];
            for k in 0..2 {
               let z = fields[ix + k].parse::<f64>().map_err(|e|
                  Error::new(
                     ErrorKind::Other,
                     format!("Could not parse float from column {} at line {} of file {}: {}", ix + 1, line + 1, file, e
                     )))?;
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

pub fn handle_cli() -> io::Result<(HtsFile, SamHeader, Config)> {
   let m = Command::new("ont_align_view")
      .version(crate_version!())
      .author("Simon Heath")
      .about("Visualize ONT alignments in a given genomic region")
      .arg(
         Arg::new("loglevel")
            .short('l')
            .long("loglevel")
            .takes_value(true)
            .value_name("LOGLEVEL")
            .help("Set log level")
            .possible_values(&["none", "error", "warn", "info", "debug", "trace"])
            .ignore_case(true),
      )
      .arg(
         Arg::new("mapq_threshold")
            .short('q')
            .long("mapq-threshold")
            .default_value("0")
            .takes_value(true)
            .value_name("MAPQ")
            .help("MAPQ quality threshold"),
      )
      .arg(
         Arg::new("max_qual")
            .short('M')
            .long("max_qual")
            .default_value("40")
            .takes_value(true)
            .value_name("QUAL")
            .help("Quality values are capped at this value"),
      )
      .arg(
         Arg::new("snv_thresholds")
            .long("snv-thresholds")
            .number_of_values(2)
            .require_value_delimiter(true)
            .default_values(&["0.0005", "0.0025"])
            .value_name("FREQ")
            .help("SNV's with frequency below first value (hard limit) are not reported, and below second value (soft limit) have a low_freq warning"),
      )
      .arg(
         Arg::new("indel_thresholds")
            .long("indel-thresholds")
            .number_of_values(2)
            .require_value_delimiter(true)
            .default_values(&["0.1", "0.2"])
            .value_name("FREQ")
            .help("Indels with freq below first value (hard limit) are not reported, and below second value (soft limit) have a low_freq warning"),
      )
      .arg(
         Arg::new("qual_threshold")
            .short('Q')
            .long("qual-threshold")
            .default_value("0")
            .takes_value(true)
            .value_name("QUAL")
            .help("Base quality threshold"),
      )
      .arg(
         Arg::new("homopolymer_limit")
            .short('P')
            .long("homopolymer-limit")
            .default_value("4")
            .takes_value(true)
            .value_name("INT")
            .help("Minimum size of homopolymer runs flagged as problematic"),
      )
      .arg(
         Arg::new("region")
            .short('r')
            .long("region")
            .takes_value(true)
            .value_name("REGION")
            .help("Genomic region to consider"),
      )
      .arg(
         Arg::new("reference")
            .short('T')
            .long("reference")
            .takes_value(true)
            .required(true)
            .value_name("FASTA File")
            .help("Reference FASTA file"),
      )
      .arg(
         Arg::new("sample")
            .short('n')
            .long("sample")
            .takes_value(true)
            .value_name("NAME")
            .help("Sample name (for VCF file)"),
      )
      .arg(
         Arg::new("blacklist")
            .long("blacklist")
            .takes_value(true)
            .value_name("BED file")
            .help("BED file with list of blacklisted sites"),
      )
      .arg(
         Arg::new("qual_calib")
            .hide(true)
            .long("qual-calib")
            .takes_value(true)
            .value_name("QCAL file")
            .help("File with quality calibration estimates"),
      )
      .arg(
         Arg::new("rs_list")
            .long("rs-list")
            .takes_value(true)
            .value_name("BED file")
            .help("BED file with dbSNP rs identifiers"),
      )
      .arg(
         Arg::new("adjust")
            .short('a')
            .long("adjust")
            .takes_value(true)
            .value_name("ADJUST")
            .default_value("0")
            .help("Adjustment to genomic position"),
      )
      .arg(
         Arg::new("paired_end")
            .short('p')
            .long("paired-end")
            .help("Illumina paired end reads"),
      )
      .arg(
         Arg::new("output_prefix")
            .short('o')
            .long("output-prefix")
            .default_value("ont-align-view")
            .takes_value(true)
            .value_name("PREFIX")
            .help("Output prefix"),
      )
      .arg(
         Arg::new("no_call")
            .long("no-call")
            .help("Do not call variants"),
      )
      .arg(
         Arg::new("view")
            .short('V')
            .long("view")
            .help("Generate pileup view"),
      )
      .arg(
         Arg::new("output_qual_calib")
            .hide(true)
            .long("output-qual-calib")
            .help("Output quality calibration values"),
      )
      .arg(
         Arg::new("input")
            .takes_value(true)
            .value_name("INPUT")
            .help("Input SAM/BAM/CRAM file(s)"),
      )
      .get_matches();

   super::log_utils::init_log(&m);

   let mapq_threshold: u8 = m.value_of_t("mapq_threshold").unwrap();
   let max_qual: u8 = m.value_of_t("max_qual").map(|x: u8| x.min((N_QUAL - 1) as u8)).unwrap();
   let qual_threshold: u8 = m.value_of_t("qual_threshold").map(|x: u8| x.min(max_qual)).unwrap();
   let qual_table = setup_qual_model();

   let homopolymer_limit: u8 = m.value_of_t("homopolymer_limit").unwrap();
   let adjust: usize = m.value_of_t("adjust").unwrap();
   let paired_end = m.is_present("paired_end");
   let no_call = m.is_present("no_call");
   let view = m.is_present("view");

   let output_qual_calib = m.is_present("output_qual_calib");

   let sample = m.value_of("sample").map(Box::from);

   let output_prefix = m.value_of("output_prefix").map(Box::from).unwrap();

   let input = m.value_of("input").unwrap_or("-");

   let mut sam_file = HtsFile::new(input, "r")?;
   debug!("Opened file {} for input", &input);

   let sam_hdr = SamHeader::read(&mut sam_file)?;
   debug!("Read SAM header");

   let ref_file = m.value_of("reference").expect("Missing reference"); // Should be enforced by clap

   let rdr = CompressIo::new().path(ref_file).bufreader()?;
   debug!("Opened {} for input", ref_file);
   let reference = Reference::from_reader(rdr, &sam_hdr)?;
   if reference.n_contigs() == 0 {
      return Err(io_err(format!(
         "No contigs read in from reference file {}",
         ref_file
      )));
   }
   debug!(
        "Reference read in successfully with {} contigs",
        reference.n_contigs()
    );

   let region = {
      if let Some(reg_str) = m.value_of("region") {
         let region = Region::from_str(reg_str, &reference)?;
         if reference.contig(region.tid()).is_none() {
            return Err(io_err(format!(
               "Region {} not found in reference file {}",
               reg_str, ref_file
            )));
         }
         region
      } else {
         if reference.n_contigs() > 1 {
            return Err(io_err(
               "Region must be specified when reference file contains multiple contigs"
                  .to_string(),
            ));
         }
         let tid = reference.contigs().keys().next().expect("Empty reference");
         Region::from_str(sam_hdr.tid2name(*tid), &reference)?
      }
   };

   let blacklist = m.value_of("blacklist").map(|file| {
      read_blacklist(file, sam_hdr.tid2name(region.tid))
         .expect("Could not read in blacklist from file")
   });

   let qual_calib = m.value_of("qual_calib").map(|file| {
      read_qual_calib(file, max_qual)
         .expect("Could not read in bias list from file")
   });

   let rs = m.value_of("rs_list").map(|file| {
      read_rslist(file, sam_hdr.tid2name(region.tid))
         .expect("Could not read in rs list from file")
   });

   let snv_thresholds: [f64; 2] = m.values_of_t::<f64>("snv_thresholds")
      .map(|v| v.try_into().expect("Couldn't convert to array"))
      .expect("Missing snv_thresholds");

   let indel_thresholds: [f64; 2] = m.values_of_t::<f64>("indel_thresholds")
      .map(|v| v.try_into().expect("Couldn't convert to array"))
      .expect("Missing indel_thresholds");

   Ok((
      sam_file,
      sam_hdr,
      Config {
         output_prefix,
         sample,
         region,
         reference,
         adjust,
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
         paired_end,
         no_call,
         output_qual_calib,
         view,
      },
   ))
}
