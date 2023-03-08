use std::fmt;

use clap::{builder::PossibleValue, ArgMatches, ValueEnum};

/// LogLevel
///
/// Represents minimum level of messages that will be logged
///
#[derive(Debug, Clone, Copy)]
pub enum LogLevel {
    Error = 0,
    Warn,
    Info,
    Debug,
    Trace,
    None,
}

impl ValueEnum for LogLevel {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Self::Error,
            Self::Warn,
            Self::Info,
            Self::Debug,
            Self::Trace,
            Self::None,
        ]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::Error => Some(PossibleValue::new("error")),
            Self::Warn => Some(PossibleValue::new("warn")),
            Self::Info => Some(PossibleValue::new("info")),
            Self::Debug => Some(PossibleValue::new("debug")),
            Self::Trace => Some(PossibleValue::new("trace")),
            Self::None => Some(PossibleValue::new("none")),
        }
    }
}

impl LogLevel {
    fn level(&self) -> usize {
        *self as usize
    }
    pub fn is_none(&self) -> bool {
        self.level() > 4
    }
    pub fn get_level(&self) -> usize {
        if self.level() > 4 {
            0
        } else {
            self.level()
        }
    }
}

impl fmt::Display for LogLevel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let level_str = ["error", "warn", "info", "debug", "trace", "none"];
        if self.level() < 6 {
            write!(f, "{}", level_str[self.level()])
        } else {
            write!(f, "unknown")
        }
    }
}

/// Initialize logging from command line arguments
pub fn init_log(m: &ArgMatches) {
    let verbose = m
        .get_one::<LogLevel>("loglevel")
        .copied()
        .expect("Missing default log level");
    let quiet = verbose.is_none();

    stderrlog::new()
        .quiet(quiet)
        .verbosity(verbose.get_level())
        .init()
        .unwrap();
}
