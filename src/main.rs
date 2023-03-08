#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

mod align_store;
mod alleles;
mod cli;
mod context;
mod deletions;
mod depth;
mod detailed_align;
mod fisher;
mod freq;
mod log_utils;
mod mann_whitney;
mod model;
mod process;
mod read;
mod reference;
mod stat_funcs;
mod vcf;

fn main() -> anyhow::Result<()> {
    let (infile, cfg) = cli::handle_cli()?;
    process::process_data(infile, cfg)
}
