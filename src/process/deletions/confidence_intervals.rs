use super::em_max::EmParam;

/// 0.5 * chisq(p=0.95, df=1)
const THRESH: f64 = 1.920729;

pub fn get_confidence_intervals(em_mle: &EmParam) -> anyhow::Result<(f64, f64)> {
    let fq_mle = em_mle.fq().ok_or(anyhow!("Invalid frequency estimates"))?;
    let max_lk = em_mle.lk().ok_or(anyhow!("Invalid likelihood"))?;

    let mut emp = vec![EmParam::new(em_mle.dels()), EmParam::new(em_mle.dels())];

    let nd = fq_mle.len();

    // Working storage
    let mut tmp_dels: Vec<usize> = Vec::with_capacity(nd - 1);
    let mut cts = vec![0.0f64; nd];

    let lower = get_root(&mut emp, em_mle, max_lk, &mut cts, &mut tmp_dels, false)?;
    let upper = get_root(&mut emp, em_mle, max_lk, &mut cts, &mut tmp_dels, true)?;
    Ok((lower, upper))
}

fn get_root<'a>(
    emp: &mut [EmParam<'a>],
    em_mle: &EmParam<'a>,
    max_lk: f64,
    cts: &mut [f64],
    tmp_dels: &mut Vec<usize>,
    upper: bool,
) -> anyhow::Result<f64> {
    emp[0].copy_into(em_mle);
    emp[1].copy_into(em_mle);

    let (i, inc) = if upper { (1, 1.0) } else { (0, 0.0) };
    if find_lim(cts, tmp_dels, &mut emp[i], max_lk, inc).is_ok() {
        bisect(
            emp,
            |e, x| {
                e.rescale(x);
                e.em_max(cts, tmp_dels, true)
            },
            |l| max_lk - l - THRESH,
        )
    } else {
        Ok(inc)
    }
}

/// On entry, emp[0] and emp[1] have the frequency estimates for the two ends of the starting interval
fn bisect<F, G>(emp: &mut [EmParam], mut f: F, g: G) -> anyhow::Result<f64>
where
    F: FnMut(&mut EmParam, f64) -> anyhow::Result<f64>,
    G: Fn(f64) -> f64,
{
    let l = g(emp[0].lk().unwrap());
    let lmid = g(emp[1].lk().unwrap());
    assert!(l * lmid <= 0.0);
    let x1 = emp[0].wt_fq().unwrap();
    let x2 = emp[1].wt_fq().unwrap();

    let (mut rtb, mut dx, mut ix) = if l < 0.0 {
        (x1, x2 - x1, 0)
    } else {
        (x2, x1 - x2, 1)
    };

    for _ in 0..1000 {
        let xmid = rtb + dx;
        dx *= 0.5;
        let lmid = g(f(&mut emp[1 - ix], xmid)?);
        if lmid <= 0.0 {
            rtb = xmid;
            ix = 1 - ix;
        }
        if dx.abs() < 1.0e-6 || lmid == 0.0 {
            break;
        }
    }
    Ok(rtb)
}

fn find_lim(
    cts: &mut [f64],
    tmp_dels: &mut Vec<usize>,
    e: &mut EmParam,
    max_lk: f64,
    inc: f64,
) -> anyhow::Result<f64> {
    let mut q = e.wt_fq().ok_or(anyhow!("Invalid freq estimates"))?;

    let lk = loop {
        let new_q = (q + inc) * 0.5;
        e.rescale(new_q);
        let lk = e.em_max(cts, tmp_dels, true)?;
        assert!(lk <= max_lk);
        if max_lk - lk > THRESH {
            break lk;
        }
        if (q - new_q).abs() < 1.0e-6 {
            return Err(anyhow!("Failed to find limit"));
        }
        q = new_q;
    };
    Ok(lk)
}
