use std::str::FromStr;

use mcircuit::Operation;
use num_bigint::BigUint;

use crate::{
    builder::{Circuit, CktBuilder, GateOperation, Template},
    gf_ckt::{GF_LEN, Gf, emit_gf_add, emit_gf_square},
    gf_mul_fft_ckt::emit_gf_mul,
};

/// A point on the xsk233 curve in internal representation
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub(crate) struct CurvePoint {
    pub x: Gf,
    pub s: Gf,
    pub z: Gf,
    pub t: Gf,
}

pub(crate) const COMPRESSED_POINT_LEN: usize = 30;
pub(crate) type CompressedCurvePoint = [[usize; 8]; COMPRESSED_POINT_LEN];
pub(crate) type CompressedCurvePointRef = [u8; COMPRESSED_POINT_LEN];

fn gf_zero<T: Circuit>(b: &mut T) -> Gf {
    let z = b.zero();
    [z; GF_LEN]
}

fn gf_to_bits(n: &BigUint) -> [u8; GF_LEN] {
    let bytes = n.to_bytes_le();
    let mut bits = [0u8; GF_LEN];
    for i in 0..GF_LEN {
        let byte = if i / 8 < bytes.len() { bytes[i / 8] } else { 0 };
        bits[i] = (byte >> (i % 8)) & 1;
    }
    bits
}

fn gf_one<T: Circuit>(b: &mut T) -> Gf {
    let z = b.zero();
    let o = b.one();
    let mut zs = [z; GF_LEN];
    zs[0] = o; // lsb at 0 bit
    zs
}

pub(crate) fn emit_point_equals<T: Circuit>(
    bld: &mut T,
    p1: &CurvePoint,
    p2: &CurvePoint,
) -> usize {
    // tmp1 = S₁·T₂
    let tmp1: Gf = emit_gf_mul(bld, &p1.s, &p2.t);

    // tmp2 = S₂·T₁
    let tmp2: Gf = emit_gf_mul(bld, &p2.s, &p1.t);

    let mut or_res = bld.zero();
    tmp1.iter().zip(tmp2).for_each(|(a, b)| {
        let a_xor_b = bld.xor_wire(*a, b);
        or_res = bld.or_wire(a_xor_b, or_res);
    });
    let one = bld.one();

    bld.xor_wire(or_res, one)
}

impl CurvePoint {
    /// Returns the identity element (point at infinity) of the curve.
    /// The exact representation depends on the coordinate system.
    /// For projective coordinates (X, S, Z), it's often (0, 1, 0).
    pub(crate) fn identity<T: Circuit>(bld: &mut T) -> Self {
        CurvePoint {
            x: gf_zero(bld),
            s: gf_one(bld), // Or Y=1
            z: gf_one(bld),
            t: gf_zero(bld), // Often X*Z = 0 for identity
        }
    }

    pub(crate) fn generator<T: Circuit>(bld: &mut T) -> Self {
        let x = BigUint::from_str(
            "13283792768796718556929275469989697816663440403339868882741001477299174",
        )
        .unwrap();
        let s = BigUint::from_str(
            "6416386389908495168242210184454780244589215014363767030073322872085145",
        )
        .unwrap();
        let z = BigUint::from_str("1").unwrap();
        let t = BigUint::from_str(
            "13283792768796718556929275469989697816663440403339868882741001477299174",
        )
        .unwrap();

        let x = gf_to_bits(&x);
        let s = gf_to_bits(&s);
        let z = gf_to_bits(&z);
        let t = gf_to_bits(&t);

        let x: Vec<usize> = x
            .iter()
            .map(|xi| if *xi == 1 { bld.one() } else { bld.zero() })
            .collect();
        let s: Vec<usize> = s
            .iter()
            .map(|xi| if *xi == 1 { bld.one() } else { bld.zero() })
            .collect();
        let z: Vec<usize> = z
            .iter()
            .map(|xi| if *xi == 1 { bld.one() } else { bld.zero() })
            .collect();
        let t: Vec<usize> = t
            .iter()
            .map(|xi| if *xi == 1 { bld.one() } else { bld.zero() })
            .collect();

        CurvePoint {
            x: x.try_into().unwrap(),
            s: s.try_into().unwrap(),
            z: z.try_into().unwrap(),
            t: t.try_into().unwrap(),
        }
    }
}

pub(crate) fn emit_point_add<T: Circuit>(
    bld: &mut T,
    p1: &CurvePoint,
    p2: &CurvePoint,
) -> CurvePoint {
    /*
     * x1x2 <- X1*X2
     * s1s2 <- S1*S2
     * z1z2 <- Z1*Z2
     * (suppressed: t1t2 <- T1*T2)
     * d <- (S1 + T1)*(S2 + T2)
     * (suppressed: e <- (a^2)*t1t2)
     * f <- x1x2^2
     * g <- z1z2^2
     * X3 <- d + s1s2
     * S3 <- sqrt(b)*(g*s1s2 + f*d) note: sqrt(b) = 1 for xsk233
     * Z3 <- sqrt(b)*(f + g)
     * T3 <- X3*Z3
     */

    // Step 1: Calculate products
    let x1x2 = emit_gf_mul(bld, &p1.x, &p2.x);
    let s1s2 = emit_gf_mul(bld, &p1.s, &p2.s);
    let z1z2 = emit_gf_mul(bld, &p1.z, &p2.z);

    // Step 2: Calculate d = (S1 + T1)*(S2 + T2)
    let tmp1 = emit_gf_add(bld, &p1.s, &p1.t);
    let tmp2 = emit_gf_add(bld, &p2.s, &p2.t);
    let d = emit_gf_mul(bld, &tmp1, &tmp2);

    // Step 3: Calculate squares
    let f = emit_gf_square(bld, &x1x2);
    let g = emit_gf_square(bld, &z1z2);

    // Step 4: Calculate output coordinates
    let p3_x = emit_gf_add(bld, &d, &s1s2);

    let tmp1 = emit_gf_mul(bld, &s1s2, &g);
    let tmp2 = emit_gf_mul(bld, &d, &f);
    let p3_s = emit_gf_add(bld, &tmp1, &tmp2);

    let p3_z = emit_gf_add(bld, &f, &g);
    let p3_t = emit_gf_mul(bld, &p3_x, &p3_z);

    CurvePoint {
        x: p3_x,
        s: p3_s,
        z: p3_z,
        t: p3_t,
    }
}

/// Apply the Frobenius endomorphism on a point (i.e. square all coordinates).
///
/// Squares all coordinates of a xsk233 curve point.
pub(crate) fn emit_point_frob<T: Circuit>(bld: &mut T, p1: &CurvePoint) -> CurvePoint {
    // Square all coordinates
    let p3_x = emit_gf_square(bld, &p1.x);
    let p3_z = emit_gf_square(bld, &p1.z);
    let p3_s = emit_gf_square(bld, &p1.s);
    let p3_t = emit_gf_square(bld, &p1.t);
    CurvePoint {
        x: p3_x,
        s: p3_s,
        z: p3_z,
        t: p3_t,
    }
}

pub(crate) fn template_emit_point_add() -> Template {
    println!("Initializing template_emit_point_add");
    let mut bld = CktBuilder::default();
    // define const wires
    let const_wire_zero = bld.zero();
    let const_wire_one = bld.one();
    let p1 = CurvePoint {
        x: bld.fresh(),
        s: bld.fresh(),
        z: bld.fresh(),
        t: bld.fresh(),
    };
    let p2 = CurvePoint {
        x: bld.fresh(),
        s: bld.fresh(),
        z: bld.fresh(),
        t: bld.fresh(),
    };
    let mut input_wires = vec![];
    input_wires.extend_from_slice(&p1.x);
    input_wires.extend_from_slice(&p1.s);
    input_wires.extend_from_slice(&p1.z);
    input_wires.extend_from_slice(&p1.t);

    input_wires.extend_from_slice(&p2.x);
    input_wires.extend_from_slice(&p2.s);
    input_wires.extend_from_slice(&p2.z);
    input_wires.extend_from_slice(&p2.t);

    let start_wire_idx = bld.next_wire();
    let res = emit_point_add(&mut bld, &p1, &p2);
    let end_wire_idx = bld.next_wire();

    let mut output_wires = vec![];
    output_wires.extend_from_slice(&res.x);
    output_wires.extend_from_slice(&res.s);
    output_wires.extend_from_slice(&res.z);
    output_wires.extend_from_slice(&res.t);

    let gates = bld.get_gates();

    let (tmp_mul, tmp_add) = {
        let mut temp_mul_gates_count = 0;
        let mut temp_xor_gates_count = 0;
        for g in gates {
            if let GateOperation::Base(bg) = g {
                match bg {
                    Operation::Add(_, _, _) => {
                        temp_xor_gates_count += 1;
                    }
                    Operation::Mul(_, _, _) => {
                        temp_mul_gates_count += 1;
                    }
                    _ => unreachable!(),
                }
            }
        }
        (temp_mul_gates_count, temp_xor_gates_count)
    };

    Template {
        input_wires: input_wires.clone(),
        output_wires,
        gates: gates.clone(),
        start_wire_idx,
        end_wire_idx,
        const_wire_one,
        const_wire_zero,
        stats: (tmp_mul, tmp_add),
    }
}

#[cfg(test)]
mod test {

    use std::{str::FromStr, time::Instant};

    use crate::builder::Circuit;
    use num_bigint::{BigUint, RandomBits};
    use rand::Rng;

    use crate::{
        builder::CktBuilder,
        curve_ref::{CurvePointRef as InnerPointRef, point_add as ref_point_add},
        gf_ref::{gfref_mul, gfref_to_bits},
    };

    use super::{CurvePoint, emit_point_add};

    // Creates a random point ensuring T = X*Z
    fn random_point() -> InnerPointRef {
        let mut rng = rand::thread_rng();
        let max_bit_len = 232;
        let x = rng.sample(RandomBits::new(max_bit_len));
        let s = rng.sample(RandomBits::new(max_bit_len));
        let z = rng.sample(RandomBits::new(max_bit_len));

        let t = gfref_mul(&x, &z);

        InnerPointRef { x, s, z, t }
    }

    #[test]
    fn test_point_add() {
        let pt = InnerPointRef {
            x: BigUint::from_str(
                "13283792768796718556929275469989697816663440403339868882741001477299174",
            )
            .unwrap(),
            s: BigUint::from_str(
                "6416386389908495168242210184454780244589215014363767030073322872085145",
            )
            .unwrap(),
            z: BigUint::from_str("1").unwrap(),
            t: BigUint::from_str(
                "13283792768796718556929275469989697816663440403339868882741001477299174",
            )
            .unwrap(),
        };

        let pt2 = random_point();
        let ptadd = ref_point_add(&pt, &pt2);

        let mut bld = CktBuilder::default();
        let c_pt = CurvePoint {
            x: bld.fresh(),
            s: bld.fresh(),
            z: bld.fresh(),
            t: bld.fresh(),
        };

        let c_pt2 = CurvePoint {
            x: bld.fresh(),
            s: bld.fresh(),
            z: bld.fresh(),
            t: bld.fresh(),
        };

        println!("Emit Pt Add START");

        let st = Instant::now();
        bld.zero();
        bld.one();
        let c_ptadd = emit_point_add(&mut bld, &c_pt, &c_pt2);
        let el = st.elapsed();
        println!("emit_point_add took {} seconds to compile ", el.as_secs());

        println!("Emit Pt Add DONE");
        let mut witness = Vec::<bool>::with_capacity(233 * 8);
        witness.extend(gfref_to_bits(&pt.x));
        witness.extend(gfref_to_bits(&pt.s));
        witness.extend(gfref_to_bits(&pt.z));
        witness.extend(gfref_to_bits(&pt.t));

        witness.extend(gfref_to_bits(&pt2.x));
        witness.extend(gfref_to_bits(&pt2.s));
        witness.extend(gfref_to_bits(&pt2.z));
        witness.extend(gfref_to_bits(&pt2.t));

        let wires = bld.eval_gates(&witness);

        let c_ptadd_x: BigUint = c_ptadd
            .x
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });

        let c_ptadd_s: BigUint = c_ptadd
            .s
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });

        let c_ptadd_z: BigUint = c_ptadd
            .z
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });

        let c_ptadd_t: BigUint = c_ptadd
            .t
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });

        let c_ptadd_val = InnerPointRef {
            x: c_ptadd_x,
            s: c_ptadd_s,
            z: c_ptadd_z,
            t: c_ptadd_t,
        };
        assert_eq!(c_ptadd_val, ptadd);

        bld.show_gate_counts();
    }
}
