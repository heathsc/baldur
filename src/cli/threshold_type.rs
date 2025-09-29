#[derive(Debug, Copy, Clone)]
pub enum ThresholdType {
    Hard,
    Soft,
}

impl ThresholdType {
    pub fn idx(&self) -> usize {
        match self {
            Self::Hard => 0,
            Self::Soft => 1,
        }
    }
}
