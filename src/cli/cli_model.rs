use std::path::PathBuf;

use clap::{Arg, ArgAction, Command, command, value_parser};

use crate::{
    log_utils::LogLevel,
    model::N_QUAL,
};

pub(super) fn cli_model() -> Command {
    command!()
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
            .default_values(["0.0005", "0.0025"])
            .value_name("FREQ1, FREQ2")
            .help("SNV's with frequency below FREQ1 (hard limit) are not reported, and below FREQ2 (soft limit) have a low_freq warning"),
    )
    .arg(
        Arg::new("indel_thresholds")
            .long("indel-thresholds")
            .num_args(2)
            .value_delimiter(',')
            .value_parser(value_parser!(f64))
            .default_values(["0.025", "0.1"])
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
    .arg(
        Arg::new("large_deletion_limit")
            .long("large-deletion-limit")
            .value_parser(value_parser!(usize))
            .value_name("SIZE")
            .default_value("64")
            .help("Minimum size for a large deletion"),
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
        Arg::new("guides")
            .short('g')
            .long("guides")
            .value_parser(value_parser!(PathBuf))
            .value_name("Guide File")
            .help("Guide definition file"),
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
}
