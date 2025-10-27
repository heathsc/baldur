use crate::reference::Reference;

use super::{Config, Guides, Region, super::ThresholdType};

use crate::model::N_QUAL;

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
    
    pub fn profile_like(&self) -> bool {
        self.profile_like
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
    
    pub fn large_deletion_limit(&self) -> usize {
        self.large_deletion_limit
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
            &v[ctxt][if q != 61 {
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
    
    pub fn guides(&self) -> Option<&Guides> {
        self.guides.as_ref()
    }
}
