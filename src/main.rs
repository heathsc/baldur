#[macro_use]
extern crate log;

use std::io;

mod cli;
mod fisher;
mod log_utils;
mod process;
mod reference;
mod model;
mod mann_whitney;
mod vcf;
mod depth;
mod context;
mod detailed_align;
mod alleles;
mod align_store;
mod read;

pub fn io_err(s: String) -> io::Error {
    io::Error::new(io::ErrorKind::Other, s)
}

fn main() -> io::Result<()> {
    let (infile, sam_hdr, cfg) = cli::handle_cli()?;
    process::process_data(infile, sam_hdr, cfg)
}
