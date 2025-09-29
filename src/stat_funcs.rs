#[link(name = "m")]
unsafe extern "C" {
   #[link_name = "erfc"]
   fn libm_erfc(x: f64) -> f64;
   #[link_name = "erf"]
   fn libm_erf(x: f64) -> f64;
}

use std::f64::consts::SQRT_2;

pub fn erf(x: f64) -> f64 {
   unsafe { libm_erf(x) }
}

pub fn erfc(x: f64) -> f64 {
   unsafe { libm_erfc(x) }
}

/// Lower tail of standard normal
pub fn pnorm(z: f64) -> f64 {
   0.5 * (1.0 + erf(z / SQRT_2))
}

/// Upper tail of standard normal
pub fn pnormc(z: f64) -> f64 {
    0.5 * erfc(z / SQRT_2)
}

/// Upper tail of chi-squared prob. of x with 1df
pub fn chisq1(x: f64) -> f64 {
   2.0 * pnorm(-x.sqrt())
}
