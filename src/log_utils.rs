use clap::ArgMatches;

/// Set up stderrlog using options from clap::ArgMatches
///
/// Relevant options:
///
///   loglevel: can be error, warn, info, trace or none.  Default is info
///   quiet: turns off all loggin irrespective of loglevel
///   timestamp: sets timestamp options for logger
///
///   loglevel none is a synonym for the quiet option
///
pub fn init_log(m: &ArgMatches) {
    let verbose = m
        .value_of("loglevel")
        .map(|s| match s.to_lowercase().as_str() {
            "error" => Some(0),
            "warn" => Some(1),
            "info" => Some(2),
            "debug" => Some(3),
            "trace" => Some(4),
            "none" => None,
            _ => Some(0),
        })
        .unwrap_or(Some(2));

    stderrlog::new()
        .verbosity(verbose.unwrap_or(0))
        .init()
        .unwrap();
}
