use std::cmp;

use libc::c_double;

#[link(name = "m")]
unsafe extern "C" {
    unsafe fn lgamma(x: c_double) -> c_double;
}

const LFACT_STORE_SIZE: usize = 256;

fn addlog(x1: f64, x2: f64) -> f64 {
    if x1 > x2 {
        let diff = x2 - x1;
        if diff < -745.0 {
            x1
        } else {
            x1 + diff.exp().ln_1p()
        }
    } else {
        let diff = x1 - x2;
        if diff < -745.0 {
            x2
        } else {
            x2 + diff.exp().ln_1p()
        }
    }
}

pub struct FisherTest {
    lfact_store: [f64; LFACT_STORE_SIZE],
}

impl FisherTest {
    pub fn new() -> Self {
        let mut lfact_store = [0.0; LFACT_STORE_SIZE];
        for i in 2..LFACT_STORE_SIZE {
            lfact_store[i] = lfact_store[i - 1] + (i as f64).ln();
        }
        Self { lfact_store }
    }

    pub fn lfact(&self, x: usize) -> f64 {
        if x < LFACT_STORE_SIZE {
            self.lfact_store[x]
        } else {
            unsafe { lgamma((x + 1) as f64) }
        }
    }

    pub fn fisher(&self, ftab: &[usize; 4]) -> f64 {
        let row = [(ftab[0] + ftab[1]) as f64, (ftab[2] + ftab[3]) as f64];
        let col = [(ftab[0] + ftab[2]) as f64, (ftab[1] + ftab[3]) as f64];
        let n = row[0] + row[1];
        let mut c = ftab.to_vec();
        if n < 1.0 {
            1.0
        } else {
            let delta = (ftab[0] as f64) - row[0] * col[0] / n;
            let konst = self.lfact(c[0] + c[2])
                + self.lfact(c[1] + c[3])
                + self.lfact(c[0] + c[1])
                + self.lfact(c[2] + c[3])
                - self.lfact(c[0] + c[1] + c[2] + c[3]);
            let mut llike =
                konst - self.lfact(c[0]) - self.lfact(c[1]) - self.lfact(c[2]) - self.lfact(c[3]);
            let mut lprob = llike;
            if delta > 0.0 {
                // Decrease counter diagonal elements until zero (this will increase delta)
                let min = cmp::min(c[1], c[2]);
                for i in 0..min {
                    llike += (((c[1] - i) * (c[2] - i)) as f64).ln()
                        - (((c[0] + i + 1) * (c[3] + i + 1)) as f64).ln();
                    lprob = addlog(lprob, llike);
                }
                let min = cmp::min(c[0], c[3]);
                // Calculate amount required to increase delta by decreasing leading diagonal elements
                let adjust = (2.0 * delta).ceil() as usize;
                if adjust <= min {
                    c[0] -= adjust;
                    c[3] -= adjust;
                    c[1] += adjust;
                    c[2] += adjust;
                    llike = konst
                        - self.lfact(c[0])
                        - self.lfact(c[1])
                        - self.lfact(c[2])
                        - self.lfact(c[3]);
                    lprob = addlog(lprob, llike);
                    for i in 0..min - adjust {
                        llike += (((c[0] - i) * (c[3] - i)) as f64).ln()
                            - (((c[1] + i + 1) * (c[2] + i + 1)) as f64).ln();
                        lprob = addlog(lprob, llike);
                    }
                }
            } else {
                // Decrease leading diagonal elements until zero (this will increase delta)
                let min = cmp::min(c[0], c[3]);
                for i in 0..min {
                    llike += (((c[0] - i) * (c[3] - i)) as f64).ln()
                        - (((c[1] + i + 1) * (c[2] + i + 1)) as f64).ln();
                    lprob = addlog(lprob, llike);
                }
                let min = cmp::min(c[1], c[2]);
                let adjust = cmp::max((-2.0 * delta).ceil() as usize, 1);
                if adjust <= min {
                    c[0] += adjust;
                    c[3] += adjust;
                    c[1] -= adjust;
                    c[2] -= adjust;
                    llike = konst
                        - self.lfact(c[0])
                        - self.lfact(c[1])
                        - self.lfact(c[2])
                        - self.lfact(c[3]);
                    lprob = addlog(lprob, llike);

                    for i in 0..min - adjust {
                        llike += (((c[1] - i) * (c[2] - i)) as f64).ln()
                            - (((c[0] + i + 1) * (c[3] + i + 1)) as f64).ln();
                        lprob = addlog(lprob, llike);
                    }
                }
            }
            lprob.exp()
        }
    }
}
