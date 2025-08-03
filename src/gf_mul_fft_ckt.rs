use crate::{
    builder::Circuit,
    gf_ckt::Gf,
    gf_interpolate_ckt::{emit_gf9_interpolate_fft, pow_table},
    gf9_ckt::{emit_gf9_mul, emit_gf9_square},
    gf9_eval_ckt::emit_poly_eval_domain,
    gf9_ref::{gf9ref_mul, gf9ref_pow},
};

/// Generate 52 orbit representatives modulo 2⁹-1=511, skipping multiples of 73
fn representatives() -> Vec<usize> {
    let mut seen = [false; 511];
    let mut reps = Vec::with_capacity(52);
    for k in 1..511 {
        if seen[k] {
            continue;
        }
        reps.push(k);
        let mut t = k;
        let rounds = if t % 73 == 0 { 3 } else { 9 };
        for _ in 0..rounds {
            seen[t] = true;
            t = (t * 2) % 511;
        }
    }
    reps
}

fn extend_domain(rep: &[usize], xs: &[u16]) -> Vec<u16> {
    let mut exts = Vec::with_capacity(511);
    exts.extend_from_slice(&xs[0..1]);
    for i in 1..xs.len() {
        let mut t = xs[i];
        let rounds = if rep[i - 1] % 73 == 0 { 3 } else { 9 };
        for _ in 0..rounds {
            exts.push(t);
            t = gf9ref_mul(t, t);
        }
    }
    exts
}

/// Build evaluation points: 0,1 and 52 orbits of length 9 each
fn domain(reps: &[usize]) -> Vec<u16> {
    let mut xs = Vec::with_capacity(465);
    xs.push(1);
    const G: u16 = 0x02;
    for &r in reps.iter() {
        let t = gf9ref_pow(G, r as u32);
        xs.push(t);
    }
    xs
}

fn extend_evaluation<T: Circuit>(
    bld: &mut T,
    reps: &[usize],
    rep_evs: &[[usize; 9]],
) -> Vec<[usize; 9]> {
    // first two evs are 0, 1
    // rest should be for representatives
    let mut ys: Vec<[usize; 9]> = Vec::with_capacity(465);
    ys.push(rep_evs[0]);

    for i in 1..rep_evs.len() {
        let mut t = rep_evs[i];
        let rounds = if reps[i - 1] % 73 == 0 { 3 } else { 9 };
        for _ in 0..rounds {
            ys.push(t);
            t = emit_gf9_square(bld, t);
        }
    }
    ys
}

/// FFT-based multiplication: BigUint → BigUint
pub(crate) fn emit_gf_mul<T: Circuit>(bld: &mut T, a: &Gf, b: &Gf) -> Gf {
    // Evaluate, multiply, interpolate
    let rep = representatives();

    let xs = domain(&rep);

    let ya = emit_poly_eval_domain(bld, a, &xs);
    let yb = emit_poly_eval_domain(bld, b, &xs);
    let yc: Vec<[usize; 9]> = ya
        .iter()
        .zip(yb.iter())
        .map(|(&u, &v)| emit_gf9_mul(bld, u, v))
        .collect();
    assert_eq!(yc.len(), 59);
    let ycs = extend_evaluation(bld, &rep, &yc);

    let xs = extend_domain(&rep, &xs); // todo: replace with cheap extend_domain
    assert_eq!(xs.len(), ycs.len());

    // reorder domain and evals in the fft domain order
    let pys = {
        let pxs = pow_table();

        let mut xs_sorted = xs.to_vec();
        xs_sorted.sort_unstable();

        let mut pxs_sorted = pow_table(); // returns Vec<T>
        pxs_sorted.sort_unstable();

        assert_eq!(
            xs_sorted, pxs_sorted,
            "pow_table() is NOT a permutation of xs"
        );

        let mut pos_map = vec![0; pxs.len()]; // pxs-index → xs-index

        for (xs_idx, &x) in xs.iter().enumerate() {
            if let Some(pxs_idx) = pxs.iter().position(|&p| p == x) {
                pos_map[pxs_idx] = xs_idx;
            }
        }

        let mut pys = vec![[0; 9]; pxs.len()];
        for (i, &xs_idx) in pos_map.iter().enumerate() {
            pys[i] = ycs[xs_idx];
        }
        pys
    };

    let pys: [[usize; 9]; 511] = pys.try_into().unwrap();

    // Convert result back to BigUint
    emit_gf9_interpolate_fft(bld, &pys)
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use crate::{
        builder::CktBuilder,
        gf_ref::{gfref_mul, gfref_to_bits},
    };

    use super::*;
    use num_bigint::{BigUint, RandomBits};
    use rand::{Rng, SeedableRng, rngs::StdRng};

    #[test]
    fn test_emit_gf_mul() {
        let mut rng = StdRng::seed_from_u64(0xC0FFEE);

        let a: u64 = rng.r#gen::<u64>() & 0x1FF; // random 9-bit ints
        let b: u64 = rng.r#gen::<u64>() & 0x1FF;
        let size = 64;

        let mut bld = CktBuilder::default();
        let a_bits: Gf = bld.fresh();
        let b_bits: Gf = bld.fresh();
        let r = Instant::now();
        let out_bits = emit_gf_mul(&mut bld, &a_bits, &b_bits);
        let el = r.elapsed();
        println!("emit_gf_mul takes {:?} seconds", el.as_secs());

        let mut witness = Vec::<bool>::with_capacity(233 * 2);
        witness.extend((0..size).map(|i| (a >> i) & 1 != 0));
        witness.extend((size..233).map(|_| false));
        witness.extend((0..size).map(|i| (b >> i) & 1 != 0));
        witness.extend((size..233).map(|_| false));

        let wires = bld.eval_gates(&witness);

        let hw: BigUint = out_bits
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });
        let chw = gfref_mul(&BigUint::from(a), &BigUint::from(b));
        bld.show_gate_counts();
        assert_eq!(hw, chw);
    }

    #[test]
    fn test_emit_gf_mul2() {
        let mut rng = StdRng::seed_from_u64(0xC0FFEE);

        let a: BigUint = rng.sample(RandomBits::new(232));
        let b: BigUint = rng.sample(RandomBits::new(232));

        let mut bld = CktBuilder::default();
        let a_bits: Gf = bld.fresh();
        let b_bits: Gf = bld.fresh();
        let out_bits = emit_gf_mul(&mut bld, &a_bits, &b_bits);

        let mut witness = Vec::<bool>::with_capacity(233 * 2);
        witness.extend(gfref_to_bits(&a));
        witness.extend(gfref_to_bits(&b));
        assert_eq!(witness.len(), 233 * 2);
        let wires = bld.eval_gates(&witness);

        let hw: BigUint = out_bits
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });
        let chw = gfref_mul(&a, &b);
        assert_eq!(hw, chw);

        bld.show_gate_counts()
    }
}
