use std::collections::{HashMap, HashSet};

use crate::{model::N_QUAL, reference::Reference};
use super::{Guides, Region};

mod getters;
mod mk_config;

pub struct Config {
    region: Region,
    reference: Reference,
    guides: Option<Guides>,
    output_prefix: Box<str>,
    sample: Option<Box<str>>,
    adjust: usize,
    small_deletion_limit: usize,
    large_deletion_limit: usize,
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
