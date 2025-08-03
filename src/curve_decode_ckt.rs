// -------------------------------------------------------------------------
//  GF(2^233) – basic constant-time arithmetic (school-book, easy to read).
// -------------------------------------------------------------------------

use crate::{
    builder::Circuit,
    curve_ckt::{CompressedCurvePoint, CurvePoint},
    gf_ckt::{
        GF_LEN, Gf, emit_gf_add, emit_gf_decode, emit_gf_halftrace, emit_gf_inv, emit_gf_square,
        emit_gf_trace,
    },
    gf_mul_fft_ckt::emit_gf_mul,
};

fn is_zero<T: Circuit>(bld: &mut T, w: Gf) -> usize {
    let mut acc = bld.zero();
    for x in w {
        acc = bld.or_wire(acc, x);
    }
    let one = bld.one();

    bld.xor_wire(acc, one)
}

pub(crate) fn emit_xsk233_decode<T: Circuit>(
    bld: &mut T,
    src: &CompressedCurvePoint,
) -> (CurvePoint, usize) {
    let mut w = [0; GF_LEN];
    let ve = emit_gf_decode(bld, &mut w, src);
    let wz = is_zero(bld, w);

    let d = emit_gf_square(bld, &w);
    let d = emit_gf_add(bld, &d, &w);
    let d_is_zero = is_zero(bld, d);
    let one = bld.one();
    let rp = bld.xor_wire(d_is_zero, one);

    let e = emit_gf_square(bld, &d);
    let e = emit_gf_inv(bld, e);
    let etr = emit_gf_trace(bld, &e);
    let rp_tr = bld.xor_wire(etr, one);
    let rp = bld.and_wire(rp, rp_tr);

    let f = emit_gf_halftrace(bld, &e);

    let x = emit_gf_mul(bld, &d, &f);
    let xtr = emit_gf_trace(bld, &x);
    let rp_tr = bld.xor_wire(xtr, one);
    let rp = bld.and_wire(rp, rp_tr);

    let g = emit_gf_halftrace(bld, &x);
    ///////

    let g = emit_gf_add(bld, &g, &w);
    let g = emit_gf_mul(bld, &g, &x);
    let g_tr = emit_gf_trace(bld, &g);
    let masked_a = {
        let mut tmp = [0; GF_LEN];
        for i in 0..GF_LEN {
            tmp[i] = bld.and_wire(x[i], g_tr);
        }
        tmp
    };
    let x = emit_gf_add(bld, &d, &masked_a);

    let s = emit_gf_square(bld, &w);
    let s = emit_gf_mul(bld, &s, &x);

    // condset
    let rp_or_wz = bld.or_wire(rp, wz);
    let r = bld.and_wire(ve, rp_or_wz);

    let mut one_233 = [bld.zero(); GF_LEN];
    one_233[0] = one;
    (
        CurvePoint {
            x,
            s,
            z: one_233,
            t: x,
        },
        r,
    )
}

#[cfg(test)]
mod tests {
    use std::os::raw::c_void;

    use xs233_sys::xsk233_generator;

    use crate::{
        builder::CktBuilder, curve_ckt::COMPRESSED_POINT_LEN,
        curve_ref::CurvePointRef as InnerPointRef,
    };
    use num_bigint::BigUint;

    use super::*;

    #[test]
    fn test_point_decode_roundtrip() {
        let mut bld = CktBuilder::default();

        let dst = {
            let mut bytes = Vec::new();
            for _ in 0..COMPRESSED_POINT_LEN {
                let byte: [usize; 8] = bld.fresh();
                bytes.push(byte);
            }
            let bytes: CompressedCurvePoint = bytes.try_into().unwrap();
            bytes
        };
        let (c_ptadd, out_ok_label) = emit_xsk233_decode(&mut bld, &dst);

        let witness = {
            unsafe {
                let pt = xsk233_generator;
                let mut dst = [0u8; COMPRESSED_POINT_LEN];
                xs233_sys::xsk233_encode(dst.as_mut_ptr() as *mut c_void, &pt);
                let mut wit = Vec::new();
                for d in dst {
                    let mut vs: Vec<bool> = (0..8).map(|i| (d >> i) & 1 != 0).collect();
                    wit.append(&mut vs);
                }
                wit
            }
        };

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

        assert_eq!(c_ptadd_val, InnerPointRef::generator());

        let out_label = wires[out_ok_label];
        println!("out_label {}", out_label);
        assert!(out_label, "should be 1 for valid input");
        bld.show_gate_counts();
    }
}
