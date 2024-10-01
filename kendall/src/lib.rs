mod avl;

use avl::*;
use std::{cmp::Ordering, collections::HashMap, hash::Hash};

const SMALL_NUMBER: f64 = 1e-9;

fn safe_div(a: f64, b: f64) -> f64 {
    a / (b + SMALL_NUMBER)
}

#[derive(PartialEq, Eq, Hash)]
struct FloatHash {
    value: u64
}

impl FloatHash {
    pub fn new(f: f64) -> Self { Self { value: f.to_bits() } }
}

fn kendall_core(x: &[f64], y: &[f64], t: &mut HashMap<FloatHash, i128>, u: &mut HashMap<FloatHash, i128>) -> Option<(i128, i128, i128)> {
    let n = x.len();
    if n != y.len() || n < 2 { return None; }
    let (mut kp, mut kn) = (0i128, 0i128);
    for i in 0..n {
        *t.entry(FloatHash::new(x[i])).or_insert(0i128) += 1;
        *u.entry(FloatHash::new(y[i])).or_insert(0i128) += 1;
        for j in i+1..n {
            let xy = (x[i] - x[j]) * (y[i] - y[j]);
            if xy > 0.0 { kp += 1; }
            else if xy < 0.0 { kn += 1; }
        }
    }

    Some((n as i128, kp, kn))
}

pub fn tau(x: &[f64], y: &[f64]) -> Option<f64> {
    let (mut t, mut u) = (HashMap::new(), HashMap::new());
    match kendall_core(x, y, &mut t, &mut u) {
        Some((n, kp, kn)) => {
            let n0 = n * (n - 1) / 2;
            let n1 = t.values().map(|&ti| ti * (ti - 1)).sum::<i128>() / 2;
            let n2 = u.values().map(|&ui| ui * (ui - 1)).sum::<i128>() / 2;
            let num = (kp - kn) as f64;
            let den = (((n0 - n1) * (n0 - n2)) as f64).sqrt();
            
            Some(safe_div(num, den))
        }
        None => None,
    }
}

pub fn zb(x: &[f64], y: &[f64]) -> Option<f64> {
    let (mut t, mut u) = (HashMap::new(), HashMap::new());
    match kendall_core(x, y, &mut t, &mut u) {
        Some((n, kp, kn)) => {
            let v0 = n * (n - 1) * (2 * n + 5);
            let vt = t.values().map(|&ti| ti * (ti - 1) * (2 * ti + 5)).sum::<i128>();
            let vu = u.values().map(|&uj| uj * (uj - 1) * (2 * uj + 5)).sum::<i128>();
            let v1 = (
                t.values().map(|&ti| ti * (ti - 1)).sum::<i128>()
                * u.values().map(|&uj| uj * (uj - 1)).sum::<i128>()
            ) as f64 / (2 * n * (n - 1)) as f64;
            let v2 = safe_div(
                (
                    t.values().map(|&ti| ti * (ti - 1) * (ti - 2)).sum::<i128>()
                    * u.values().map(|&uj| uj * (uj - 1) * (uj - 2)).sum::<i128>()
                ) as f64,
                (9 * n  * (n - 1) * (n - 2)) as f64);
            let v = (v0 - vt - vu) as f64 / 18.0 + v1 + v2;

            Some(safe_div((kp - kn) as f64, v.sqrt()))
        }
        None => None,
    }
}

pub fn tau_avltree(x: &[f64], y: &[f64]) -> Option<f64> {
    let n = x.len();
    if n != y.len() || n < 2 { return None; }
    let mut xy = x.iter().zip(y.iter()).map(|(&xi, &yi)| (xi, yi)).collect::<Vec<(f64, f64)>>();
    xy.sort_by(|&(xi, _), &(xj, _)| xi.partial_cmp(&xj).unwrap_or(Ordering::Greater));
    let mut tree = AVLTree::new();
    let c = xy.iter().fold(0, |m, &(_, yi)| m + tree.insert(yi) as i128);
    let num = 4.0 * c as f64;
    let den = (n as i128 * (n as i128 - 1)) as f64;

    Some(safe_div(num, den) - 1.0)
}

pub fn zb_avltree(x: &[f64], y: &[f64]) -> Option<f64> {
    let n = x.len() as i128;
    match tau_avltree(x, y) {
        Some(tau) => {
            let num = (9 * n * (n - 1)) as f64;
            let den = (2 * (2 * n + 5)) as f64;
            Some(tau * safe_div(num, den).sqrt())
        }
        None => None
    }
}

pub fn binary_core(x: &[bool], y: &[bool]) -> Option<(i128, i128, i128, i128, i128)> {
    let n = x.len();
    if n != y.len() || n < 2 { return None; }
    let (mut a, mut s, mut dx, mut dy) = (0i128, 0i128, 0i128, 0i128);
    let (mut sx, mut sy) = (0i128, 0i128);
    for i in 0..n {
        match (x[i], y[i]) {
            (true, true) => { a += 1; sx += 1; sy += 1; }
            (true, false) => { dx += 1; sx += 1; }
            (false, true) => { dy += 1; sy += 1; }
            (false, false) => { s += 1; }
        }
    }
    let kp = a * s;
    let kn = dx * dy;

    Some((n as i128, kp, kn, sx, sy))
}

pub fn tau_binary(x: &[bool], y: &[bool]) -> Option<f64> {
    match binary_core(x, y) {
        Some((n, kp, kn, sx, sy)) => {
            let n0 = n * (n - 1) / 2;
            let n1 = (sx * (sx - 1) + (n - sx) * (n - sx - 1)) / 2;
            let n2 = (sy * (sy - 1) + (n - sy) * (n - sy - 1)) / 2;
            let num = (kp - kn) as f64;
            let den = (((n0 - n1) * (n0 - n2)) as f64).sqrt();

            Some(safe_div(num, den))
        }
        None => None,
    }
}

pub fn zb_binary(x: &[bool], y: &[bool]) -> Option<f64> {
    match binary_core(x, y) {
        Some((n, kp, kn, sx, sy)) => {
            let (t, u) = ([n - sx, sx], [n - sy, sy]);
            let v0 = n * (n - 1) * (2 * n + 5);
            let vt = t.iter().map(|ti| ti * (ti - 1) * (2 * ti + 5)).sum::<i128>();
            let vu = u.iter().map(|ui| ui * (ui - 1) * (2 * ui + 5)).sum::<i128>();
            let v1 = (
                t.iter().map(|&ti| ti * (ti - 1)).sum::<i128>()
                * u.iter().map(|&uj| uj * (uj - 1)).sum::<i128>()
            ) as f64 / (2 * n * (n - 1)) as f64;
            let v2 = safe_div(
                (
                    t.iter().map(|&ti| ti * (ti - 1) * (ti - 2)).sum::<i128>()
                    * u.iter().map(|&uj| uj * (uj - 1) * (uj - 2)).sum::<i128>()
                ) as f64,
                (9 * n  * (n - 1) * (n - 2)) as f64);
            let v = (v0 - vt - vu) as f64 / 18.0 + v1 + v2;

            Some(safe_div((kp - kn) as f64, v.sqrt()))
        }
        None => None,
    }
}

pub fn zero_core(x: &[f64], y: &[f64]) -> Option<(i128, i128, i128, i128, i128)> {
    let n = x.len();
    if n != y.len() || n < 2 { return None; }
    let (mut xy, mut xz, mut zy, mut nzz) = (Vec::new(), Vec::new(), Vec::new(), 0i128);
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        if xi > 0.0 {
            if yi > 0.0 { xy.push((xi, yi)); }
            else { xz.push(xi); }
        } else {
            if yi > 0.0 { zy.push(yi); }
            else { nzz += 1; }
        }
    }
    let (nxy, nxz, nzy) = (xy.len() as i128, xz.len() as i128, zy.len() as i128);
    let t0 = nzz + nzy;
    let u0 = nzz + nxz;
    let (mut kp, mut kn) = (nxy * nzz, nxz * nzy);
    xy.sort_by(|&(xi, _), &(xj, _)| xi.partial_cmp(&xj).unwrap_or(Ordering::Greater));
    let mut tree = AVLTree::new();
    let c = xy.iter().fold(0, |m, &(_, yi)| m + tree.insert(yi) as i128);
    kp += c;
    kn += nxy * (nxy - 1) / 2 - c;

    xz.sort_by(|&xi, xj| xi.partial_cmp(xj).unwrap_or(Ordering::Greater));
    zy.sort_by(|&yi, yj| yi.partial_cmp(yj).unwrap_or(Ordering::Greater));
    let mut count = 0;
    for &xi in xz.iter() {
        loop {
            if count >= nxy || xi <= xy[count as usize].0 { break; }
            count += 1;
        }
        kp += nxy - count;
        kn += count;
    }
    xy.sort_by(|&(_, yi), &(_, yj)| yi.partial_cmp(&yj).unwrap_or(Ordering::Greater));
    let mut count = 0;
    for &yi in zy.iter() {
        loop {
            if count >= nxy || yi <= xy[count as usize].1 { break; }
            count += 1;
        }
        kp += nxy - count;
        kn += count;
    }
    
    Some((n as i128, kp, kn, t0, u0))
}

pub fn tau_zero(x: &[f64], y: &[f64]) -> Option<f64> {
    match zero_core(x, y) {
        Some((n, kp, kn, t0, u0)) => {
            let n0 = n * (n - 1) / 2;
            let n1 = t0 * (t0 - 1) / 2;
            let n2 = u0 * (u0 - 1) / 2;
            let num = (kp - kn) as f64;
            let den = (((n0 - n1) * (n0 - n2)) as f64).sqrt();
            
            Some(safe_div(num, den))
        }
        None => None,
    }
}

pub fn zb_zero(x: &[f64], y: &[f64]) -> Option<f64> {
    match zero_core(x, y) {
        Some((n, kp, kn, t0, u0)) => {
            let v0 = n * (n - 1) * (2 * n + 5);
            let vt = t0 * (t0 - 1) * (2 * t0 + 5);
            let vu = u0 * (u0 - 1) * (2 * u0 + 5);
            let v1 = (t0 * (t0 - 1) * u0 * (u0 - 1)) as f64 / (2 * n * (n - 1)) as f64;
            let v2 = safe_div(
                (t0 * (t0 - 1) * (t0 - 2) * u0 * (u0 - 1) * (u0 - 2)) as f64,
                (9 * n * (n - 1) * (n - 2)) as f64);
            let v = (v0 - vt - vu) as f64 / 18.0 + v1 + v2;

            Some(safe_div((kp - kn) as f64, v.sqrt()))
        }
        None => None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    #[test]
    fn same_array() {
        let mut rng = thread_rng();
        let x = (0..100).map(|_| rng.gen()).collect::<Vec<f64>>();
        assert!((tau(&x, &x).unwrap() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn parallel_array() {
        let x = (0..100).map(|_| 1.0).collect::<Vec<f64>>();
        let y = (0..100).map(|_| 0.0).collect::<Vec<f64>>();
        assert!((tau(&x, &y).unwrap() - 0.0).abs() < 1e-3);
    }

    #[test]
    fn tau_avltree_test() {
        let mut rng = thread_rng();
        let x = (0..100).map(|_| rng.gen()).collect::<Vec<f64>>();
        let y = (0..100).map(|_| rng.gen()).collect::<Vec<f64>>();
        assert!((tau_avltree(&x, &y).unwrap() - tau(&x, &y).unwrap()).abs() < 1e-3);
    }

    #[test]
    fn zb_avltree_test() {
        let mut rng = thread_rng();
        let x = (0..100).map(|_| rng.gen()).collect::<Vec<f64>>();
        let y = (0..100).map(|_| rng.gen()).collect::<Vec<f64>>();
        assert!((zb_avltree(&x, &y).unwrap() - zb(&x, &y).unwrap()).abs() < 1e-3);
    }

    #[test]
    fn tau_binary_test() {
        let mut rng = thread_rng();
        let xb = (0..100).map(|_| rng.gen_bool(0.5)).collect::<Vec<bool>>();
        let yb = (0..100).map(|_| rng.gen_bool(0.5)).collect::<Vec<bool>>();
        let xf = xb.iter().map(|&xi| if xi { 1.0 } else { 0.0 }).collect::<Vec<f64>>();
        let yf = yb.iter().map(|&yi| if yi { 1.0 } else { 0.0 }).collect::<Vec<f64>>();
        assert!((tau_binary(&xb, &yb).unwrap() - tau(&xf, &yf).unwrap()).abs() < 1e-3);
    }

    #[test]
    fn zb_binary_test() {
        let mut rng = thread_rng();
        let xb = (0..100).map(|_| rng.gen_bool(0.5)).collect::<Vec<bool>>();
        let yb = (0..100).map(|_| rng.gen_bool(0.5)).collect::<Vec<bool>>();
        let xf = xb.iter().map(|&xi| if xi { 1.0 } else { 0.0 }).collect::<Vec<f64>>();
        let yf = yb.iter().map(|&yi| if yi { 1.0 } else { 0.0 }).collect::<Vec<f64>>();
        assert!((zb_binary(&xb, &yb).unwrap() - zb(&xf, &yf).unwrap()).abs() < 1e-3);
    }

    #[test]
    fn tau_zero_test() {
        let mut rng = thread_rng();
        let x = (0..100).map(|_| if rng.gen_bool(0.5) { rng.gen() } else { 0.0 }).collect::<Vec<f64>>();
        let y = (0..100).map(|_| if rng.gen_bool(0.5) { rng.gen() } else { 0.0 }).collect::<Vec<f64>>();
        assert!((tau_zero(&x, &y).unwrap() - tau(&x, &y).unwrap()).abs() < 1e-3)
    }

    #[test]
    fn zb_zero_test() {
        let mut rng = thread_rng();
        let x = (0..100).map(|_| if rng.gen_bool(0.5) { rng.gen() } else { 0.0 }).collect::<Vec<f64>>();
        let y = (0..100).map(|_| if rng.gen_bool(0.5) { rng.gen() } else { 0.0 }).collect::<Vec<f64>>();
        assert!((zb_zero(&x, &y).unwrap() - zb(&x, &y).unwrap()).abs() < 1e-3)
    }
}