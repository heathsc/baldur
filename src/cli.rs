use r_htslib::Hts;

mod cli_model;
mod config;
mod guide;
mod region;
mod threshold_type;

pub use config::Config;
pub use region::Region;
pub use guide::{Guide, Guides};

pub use threshold_type::ThresholdType;

pub fn handle_cli() -> anyhow::Result<(Hts, Config)> {
    let m = cli_model::cli_model().get_matches();

    super::log_utils::init_log(&m);
    
    Config::from_matches(&m)
}