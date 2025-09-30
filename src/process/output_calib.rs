use std::{
    fs::File,
    io::{BufWriter, Write},
};

use super::ProcWork;
use crate::{
    cli::Config,
    context::{Ctxt5, N_CTXT},
};

pub(super) fn output_calibration_data(cfg: &Config, pw: &ProcWork) -> anyhow::Result<()> {
    if cfg.output_qual_calib() {
        let calc_q = |q: &[usize]| {
            let p = (q[1] + 1) as f64 / ((q[0] + q[1] + 2) as f64);
            -10.0 * p.log10()
        };
        let qual_cal_output = format!("{}_qcal.txt", cfg.output_prefix());
        let mut wrt = BufWriter::new(File::create(&qual_cal_output)?);
        write!(
            wrt,
            "Qual\tEmp_Qual\tEmp_Qual_Del\tMatch\tMismatch\tMatch_Del\tMismatch_Del"
        )?;
        for i in 0..N_CTXT {
            write!(wrt, "\t{:#}\t\t\t\t\t", Ctxt5((i as u16) << 2))?;
        }
        writeln!(wrt)?;
        for (q, qc) in pw.qual_hist.iter().enumerate() {
            if qc[0] + qc[1] + qc[2] + qc[3] > 0 {
                let empirical_q = calc_q(&qc[..2]);
                let empirical_qd = calc_q(&qc[2..]);
                write!(
                    wrt,
                    "{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}",
                    q, empirical_q, empirical_qd, qc[0], qc[1], qc[2], qc[3]
                )?;
                for cc in pw.ctxt_hist.iter() {
                    let qc1 = &cc[q];
                    let empirical_q = calc_q(&qc1[..2]);
                    let empirical_qd = calc_q(&qc1[2..]);
                    write!(
                        wrt,
                        "\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}",
                        empirical_q, empirical_qd, qc1[0], qc1[1], qc1[2], qc1[3]
                    )?;
                }
                writeln!(wrt)?;
            }
        }
    }
    Ok(())
}
