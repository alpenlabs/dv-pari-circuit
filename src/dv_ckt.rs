//! Binary circuit implementation of DV Verifier Program
//!
use std::str::FromStr;

use crate::{
    blake3_ckt,
    builder::{Circuit, CktBuilder},
    curve_ckt::{
        CompressedCurvePoint, CompressedCurvePointRef, CurvePoint, emit_point_add,
        emit_point_equals, emit_xsk233_decode,
    },
    curve_scalar_mul_ckt::point_scalar_mul::emit_mul_windowed_tau,
    fr_ckt::{
        FR_LEN, Fr, const_mod_n, emit_fr_add as fr_add, emit_fr_mul as fr_mul,
        emit_fr_sub as fr_sub, ge_unsigned,
    },
    fr_ref::{FrRef, frref_to_bits},
};

const PUBLIC_INPUT_LEN: usize = 64;

#[derive(Debug)]
/// VerifierPayloadRef
pub struct VerifierPayloadRef {
    /// proof
    pub proof: ProofRef,
    /// public_input
    pub public_input: RawPublicInputsRef,
    /// trapdoor
    pub trapdoor: TrapdoorRef,
}

#[derive(Debug)]
/// ProofRef
pub struct ProofRef {
    /// commit_p
    pub commit_p: CompressedCurvePointRef, // commitment to witness folding & quotient
    /// kzg_k
    pub kzg_k: CompressedCurvePointRef, // combined KZG evaluation proof
    /// a0
    pub a0: FrRef,
    /// b0
    pub b0: FrRef,
}

const PROOF_BIT_LEN: usize = 944;
const PUBINP_BIT_LEN: usize = 64;
const TRAPDOOR_BIT_LEN: usize = 696;

impl ProofRef {
    // (2 * 30 + 2 * 29)*8 = 944
    /// Serialize ProofRef
    pub fn to_bits(&self) -> [bool; PROOF_BIT_LEN] {
        let mut commit_p: Vec<bool> = self
            .commit_p
            .iter()
            .flat_map(|x| u8_to_bits_le(*x).to_vec())
            .collect();
        let mut kzg_k: Vec<bool> = self
            .kzg_k
            .iter()
            .flat_map(|x| u8_to_bits_le(*x).to_vec())
            .collect();
        let mut a0 = frref_to_bits(&self.a0).to_vec();
        let mut b0 = frref_to_bits(&self.b0).to_vec();

        let mut witness: Vec<bool> = Vec::new();

        witness.append(&mut commit_p);
        witness.append(&mut kzg_k);
        witness.append(&mut a0);
        witness.append(&mut b0);

        let witness: [bool; PROOF_BIT_LEN] = witness.try_into().unwrap();
        witness
    }
}

/// RawPublicInputsRef
#[derive(Debug)]
pub struct RawPublicInputsRef {
    /// deposit_index
    pub deposit_index: u64,
}

impl RawPublicInputsRef {
    // deposit_index is 64 bits
    /// Serialize RawPublicInputsRef
    pub fn to_bits(&self) -> [bool; PUBINP_BIT_LEN] {
        let deposit_index = u64_to_bits_le(self.deposit_index).to_vec();
        let deposit_index: [bool; 64] = deposit_index.try_into().unwrap();
        deposit_index
    }
}

/// TrapdoorRef
#[derive(Debug)]
pub struct TrapdoorRef {
    /// tau
    pub tau: FrRef, // trapdoor (only known to verifier)
    /// delta
    pub delta: FrRef,
    /// epsilon
    pub epsilon: FrRef,
}

impl TrapdoorRef {
    /// Serialize TrapdoorRef
    pub fn to_bits(&self) -> [bool; TRAPDOOR_BIT_LEN] {
        let mut tau = frref_to_bits(&self.tau).to_vec();
        let mut delta = frref_to_bits(&self.delta).to_vec();
        let mut epsilon = frref_to_bits(&self.epsilon).to_vec();

        let mut witness: Vec<bool> = Vec::new();
        witness.append(&mut tau);
        witness.append(&mut delta);
        witness.append(&mut epsilon);

        let witnes: [bool; TRAPDOOR_BIT_LEN] = witness.try_into().unwrap();
        witnes
    }
}

fn u64_to_bits_le(n: u64) -> [bool; 64] {
    let v: Vec<bool> = (0..64).map(|i| (n >> i) & 1 != 0).collect();
    v.try_into().unwrap()
}

fn u8_to_bits_le(n: u8) -> [bool; 8] {
    let v: Vec<bool> = (0..8).map(|i| (n >> i) & 1 != 0).collect();
    v.try_into().unwrap()
}

impl VerifierPayloadRef {
    fn get_labels(bld: &mut CktBuilder) -> (Proof, RawPublicInputs, Trapdoor) {
        let secrets = Trapdoor {
            tau: bld.fresh(),
            delta: bld.fresh(),
            epsilon: bld.fresh(),
        };
        let rpin = RawPublicInputs {
            deposit_index: bld.fresh(),
        };
        let commit_p = {
            let r: [usize; 240] = bld.fresh();
            let r: Vec<[usize; 8]> = r
                .chunks(8)
                .map(|x| {
                    let y: [usize; 8] = x.try_into().unwrap();
                    y
                })
                .collect();
            let r: [[usize; 8]; 30] = r.try_into().unwrap();
            r
        };
        let kzg_k = {
            let r: [usize; 240] = bld.fresh();
            let r: Vec<[usize; 8]> = r
                .chunks(8)
                .map(|x| {
                    let y: [usize; 8] = x.try_into().unwrap();
                    y
                })
                .collect();
            let r: [[usize; 8]; 30] = r.try_into().unwrap();
            r
        };
        let proof = Proof {
            commit_p,
            kzg_k,
            a0: bld.fresh(),
            b0: bld.fresh(),
        };
        (proof, rpin, secrets)
    }

    /// Serialize VerifierPayloadRef
    pub fn to_bits(&self) -> [bool; PROOF_BIT_LEN + PUBINP_BIT_LEN + TRAPDOOR_BIT_LEN] {
        let mut secret_bits = self.trapdoor.to_bits().to_vec();
        let mut deposit_index_bits = self.public_input.to_bits().to_vec();
        let mut proof_bits = self.proof.to_bits().to_vec();

        let mut witness: Vec<bool> = Vec::new();

        witness.append(&mut secret_bits);
        witness.append(&mut deposit_index_bits);
        witness.append(&mut proof_bits);

        let witness: [bool; PROOF_BIT_LEN + PUBINP_BIT_LEN + TRAPDOOR_BIT_LEN] =
            witness.try_into().unwrap();
        witness
    }
}

/// Proof
#[derive(Debug)]
pub struct Proof {
    /// commit_p
    pub commit_p: CompressedCurvePoint, // commitment to witness folding & quotient
    /// kzg_k
    pub kzg_k: CompressedCurvePoint, // combined KZG evaluation proof
    /// a0
    pub a0: Fr,
    /// b0
    pub b0: Fr,
}

/// RawPublicInputs
#[derive(Debug)]
pub struct RawPublicInputs {
    /// deposit_index
    pub deposit_index: [usize; PUBLIC_INPUT_LEN],
}

/// Trapdoor
#[derive(Debug)]
pub struct Trapdoor {
    /// tau
    pub tau: Fr, // trapdoor (only known to verifier)
    /// delta
    pub delta: Fr,
    /// epsilon
    pub epsilon: Fr,
}

/// Label Info
// Fr:232 x 5 + Pt:240 x 2 + Pub:64 = 1704 bits
// 0 -> const_zero
// 1 -> const_one
// (2, 1705) -> input wires
#[derive(Debug)]
pub struct LabelInfo {
    /// input wire inclusive range
    pub input_wire_range: (usize, usize),
    /// const zero wire label
    pub const_zero: usize,
    /// const one wire label
    pub const_one: usize,
    /// output label
    pub output_label: usize,
}

fn u8_arr_to_labels_le<T: Circuit>(bld: &mut T, ns: &[u8]) -> Vec<[usize; 8]> {
    let zero = bld.zero();
    let one = bld.one();
    let mut vs: Vec<[usize; 8]> = Vec::new();
    for n in ns {
        let v: Vec<usize> = (0..8)
            .map(|i| {
                let r = (n >> i) & 1 != 0;
                if r { one } else { zero }
            })
            .collect();
        let v: [usize; 8] = v.try_into().unwrap();
        vs.push(v);
    }
    vs
}

fn get_fs_challenge<T: Circuit>(
    bld: &mut T,
    commit_p: CompressedCurvePoint,
    public_inputs: [Fr; 2],
    srs_bytes: Vec<u8>,
    circuit_info_bytes: Vec<u8>,
) -> Fr {
    // convert Vec<u8> into its usize version
    let witness_commitment_hash = blake3_ckt::hash(bld, commit_p.to_vec());

    let public_inputs_hash = {
        let mut buf = Vec::new();
        for pubin in public_inputs {
            let mut pubin_240 = [bld.zero(); FR_LEN];
            pubin_240.copy_from_slice(&pubin[0..FR_LEN]);
            let mut r: Vec<[usize; 8]> = pubin_240
                .chunks(8)
                .map(|x| {
                    let y: [usize; 8] = x.try_into().unwrap();
                    y
                })
                .collect();
            buf.append(&mut r);
        }

        blake3_ckt::hash(bld, buf)
    };

    let compile_time_hash = {
        let srs_bytes = u8_arr_to_labels_le(bld, &srs_bytes);
        let srs_hash = blake3_ckt::hash(bld, srs_bytes);
        let circuit_info_bytes = u8_arr_to_labels_le(bld, &circuit_info_bytes);
        let circuit_info_hash = blake3_ckt::hash(bld, circuit_info_bytes);
        let mut compile_time_bytes: Vec<[usize; 8]> = Vec::new();
        compile_time_bytes.extend_from_slice(&srs_hash);
        compile_time_bytes.extend_from_slice(&circuit_info_hash);
        blake3_ckt::hash(bld, compile_time_bytes)
    };

    let runtime_hash = {
        let mut runtime_bytes: Vec<[usize; 8]> = Vec::new();
        runtime_bytes.extend_from_slice(&witness_commitment_hash);
        runtime_bytes.extend_from_slice(&public_inputs_hash);

        blake3_ckt::hash(bld, runtime_bytes)
    };

    let mut root_hash = {
        let mut root_bytes: Vec<[usize; 8]> = Vec::new();
        root_bytes.extend_from_slice(&compile_time_hash);
        root_bytes.extend_from_slice(&runtime_hash);

        blake3_ckt::hash(bld, root_bytes)
    };

    let zero_byte = [bld.zero(); 8];
    // truncate msb
    root_hash[28..].copy_from_slice(&[zero_byte; 4]); // mask top 3 bytes, 256-24=232 bits

    // convert to Fr
    let root_hash_flat: Vec<usize> = root_hash.into_iter().flatten().collect();
    let root_fr: Fr = root_hash_flat[0..FR_LEN].try_into().unwrap();
    root_fr
}

/// Convert raw public inputs into scalar field element as done by SP1
fn get_pub_hash_from_raw_pub_inputs<T: Circuit>(bld: &mut T, raw_pub_in: &RawPublicInputs) -> Fr {
    let inps: Vec<[usize; 8]> = raw_pub_in
        .deposit_index
        .chunks(8)
        .map(|x| {
            let y: [usize; 8] = x.try_into().unwrap();
            y
        })
        .collect();
    let mut out_hash = blake3_ckt::hash(bld, inps);
    let zero_byte = [bld.zero(); 8];
    out_hash[0..4].copy_from_slice(&[zero_byte; 4]); // MSB masked assuming BE
    out_hash.reverse();

    let out_hash_flat: Vec<usize> = out_hash.into_iter().flatten().collect();
    let out_fr: Fr = out_hash_flat[0..FR_LEN].try_into().unwrap();
    out_fr
}

fn const_biguint_to_labels<T: Circuit>(bld: &mut T, num: FrRef) -> Fr {
    let num_bits = frref_to_bits(&num);
    let r: Vec<usize> = num_bits
        .iter()
        .map(|xi| if *xi { bld.one() } else { bld.zero() })
        .collect();
    let r: Fr = r.try_into().unwrap();
    r
}

/// Function to compile dvsnark verifier circuit
pub fn compile_verifier(sp1_vk: &str) -> (CktBuilder, LabelInfo) {
    let mut bld = CktBuilder::default();

    let input_wire_start = bld.next_wire();
    let (proof, rpin, secrets) = VerifierPayloadRef::get_labels(&mut bld);
    let input_wire_end = bld.next_wire();

    // Prepare
    let passed_label = verify(&mut bld, proof, rpin, sp1_vk, secrets);
    let label_info = LabelInfo {
        input_wire_range: (input_wire_start, input_wire_end - 1), // -1 because inclusive range
        const_zero: bld.zero(),
        const_one: bld.one(),
        output_label: passed_label,
    };
    (bld, label_info)
}

/// evaluate verifier
pub fn evaluate_verifier(bld: &mut CktBuilder, witness: [bool; 1704], output_label: usize) -> bool {
    let wires = bld.eval_gates(&witness);

    wires[output_label]
}

/// verify
pub(crate) fn verify<T: Circuit>(
    bld: &mut T,
    proof: Proof,
    raw_public_inputs: RawPublicInputs,
    sp1_vk: &str,
    secrets: Trapdoor,
) -> usize {
    let (proof_commit_p, decode_proof_commit_p_success) = emit_xsk233_decode(bld, &proof.commit_p);
    let (proof_kzg_k, decode_proof_kzg_k_success) = emit_xsk233_decode(bld, &proof.kzg_k);

    let is_input_valid = {
        let one_wire = bld.one();
        let fr_modulus = const_mod_n(bld);
        let proof_a0_invalid = ge_unsigned(bld, &proof.a0, &fr_modulus); // a0 should be less than modulus
        let proof_b0_invalid = ge_unsigned(bld, &proof.b0, &fr_modulus);
        let proof_scalars_invalid = bld.or_wire(proof_a0_invalid, proof_b0_invalid); // either invalid
        let proof_scalars_valid = bld.xor_wire(proof_scalars_invalid, one_wire); // both scalars valid
        let is_pts_valid = bld.and_wire(decode_proof_commit_p_success, decode_proof_kzg_k_success); // pts valid
        bld.and_wire(proof_scalars_valid, is_pts_valid) // points and scalars valid
    };

    let public_inputs_1 = get_pub_hash_from_raw_pub_inputs(bld, &raw_public_inputs);
    let public_inputs_0_vk_const = {
        let num = FrRef::from_str(sp1_vk).unwrap(); // vk
        const_biguint_to_labels(bld, num)
    };

    let fs_challenge_alpha = get_fs_challenge(
        bld,
        proof.commit_p,
        [public_inputs_0_vk_const, public_inputs_1],
        vec![],
        vec![],
    );

    let i0 = {
        let t0 = fr_mul(bld, &public_inputs_1, &fs_challenge_alpha);
        fr_add(bld, &t0, &public_inputs_0_vk_const)
    };

    let r0 = {
        //&proof.a0 * &proof.b0 - &proof.i0
        let t0 = fr_mul(bld, &proof.a0, &proof.b0);

        fr_sub(bld, &t0, &i0)
    };

    // Step 3. Compute u₀ and v₀
    let delta2 = fr_mul(bld, &secrets.delta, &secrets.delta);
    let u0 = {
        //(&proof.a0 + &secrets.delta * &proof.b0 + &delta2 * &r0) * &secrets.epsilon;
        let db0 = fr_mul(bld, &secrets.delta, &proof.b0);
        let a0_p_db0 = fr_add(bld, &proof.a0, &db0);
        let d2_r0 = fr_mul(bld, &delta2, &r0);
        let t1 = fr_add(bld, &a0_p_db0, &d2_r0);

        fr_mul(bld, &t1, &secrets.epsilon)
    };
    let tmp0 = fr_sub(bld, &secrets.tau, &fs_challenge_alpha);
    let v0 = fr_mul(bld, &tmp0, &secrets.epsilon);

    let w = 5;
    let generator = CurvePoint::generator(bld);
    let v0_k = emit_mul_windowed_tau(bld, &v0, &proof_kzg_k, w);
    let u0_g = emit_mul_windowed_tau(bld, &u0, &generator, w);
    let lhs = emit_point_add(bld, &v0_k, &u0_g);
    let rhs: CurvePoint = proof_commit_p;

    let verify_success = emit_point_equals(bld, &lhs, &rhs);

    bld.and_wire(verify_success, is_input_valid)
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use num_bigint::BigUint;

    use crate::{
        builder::{Circuit, CktBuilder},
        dv_ckt::{
            PUBLIC_INPUT_LEN, ProofRef, RawPublicInputsRef, TrapdoorRef, VerifierPayloadRef,
            compile_verifier, evaluate_verifier, u8_to_bits_le, u64_to_bits_le,
        },
        dv_ref::{self},
        fr_ckt::Fr,
        fr_ref::{FrRef, frref_to_bits},
    };

    use super::{RawPublicInputs, get_fs_challenge, get_pub_hash_from_raw_pub_inputs};

    #[test]
    #[ignore] // ignore because of being long running
    fn test_verify_over_mock_inputs() {
        const SP1_VK_BLAKE3_FIBO: &str =
            "3211415145189105019978395454939297393438090639215215871422959911100063";
        let (mut bld, label_info) = compile_verifier(SP1_VK_BLAKE3_FIBO);

        // Prepare VerifierPayloadRef
        let tau: FrRef = BigUint::from_str(
            "2472308663339583895498147954222995510858962633570970238431638506807949",
        )
        .unwrap();
        let delta = BigUint::from_str(
            "2194316856053929337106370775922152179496555179813841848311939628788959",
        )
        .unwrap();
        let epsilon = BigUint::from_str(
            "1154785560216858119874588837659951154401760642599649999302917233356517",
        )
        .unwrap();
        let deposit_index = 2838366633794829684;
        let commit_p: [u8; 30] = [
            145, 195, 86, 210, 230, 219, 176, 179, 148, 236, 194, 133, 166, 240, 60, 111, 152, 154,
            62, 190, 248, 224, 197, 250, 131, 57, 145, 237, 213, 0,
        ];
        let kzg_k: [u8; 30] = [
            239, 6, 89, 163, 169, 250, 184, 159, 153, 181, 70, 47, 167, 56, 153, 92, 52, 197, 196,
            244, 10, 197, 235, 26, 46, 57, 18, 194, 56, 0,
        ];
        let a0 = BigUint::from_str(
            "3042729463975785669077695901360320813980996043134603468597671969223884",
        )
        .unwrap();
        let b0 = BigUint::from_str(
            "3099898550361810144312020021372781014260137110595159844100103797269587",
        )
        .unwrap();

        let verifier_payload: VerifierPayloadRef = VerifierPayloadRef {
            proof: ProofRef {
                commit_p,
                kzg_k,
                a0,
                b0,
            },
            public_input: RawPublicInputsRef { deposit_index },
            trapdoor: TrapdoorRef {
                tau,
                delta,
                epsilon,
            },
        };
        let witness = verifier_payload.to_bits();

        bld.show_gate_counts();

        println!("label_info {:?}", label_info);
        let passed_val = evaluate_verifier(&mut bld, witness, label_info.output_label);
        assert!(passed_val, "verification failed");
    }

    #[test]
    #[ignore] // ignore because of being long running
    fn test_invalid_proof_over_mock_inputs() {
        const SP1_VK_BLAKE3_FIBO: &str =
            "3211415145189105019978395454939297393438090639215215871422959911100063";
        let (mut bld, label_info) = compile_verifier(SP1_VK_BLAKE3_FIBO);

        // Prepare VerifierPayloadRef
        let tau: FrRef = BigUint::from_str(
            "2472308663339583895498147954222995510858962633570970238431638506807949",
        )
        .unwrap();
        let delta = BigUint::from_str(
            "2194316856053929337106370775922152179496555179813841848311939628788959",
        )
        .unwrap();
        let epsilon = BigUint::from_str(
            "1154785560216858119874588837659951154401760642599649999302917233356517",
        )
        .unwrap();
        let deposit_index = 2838366633794829684;
        let commit_p: [u8; 30] = [
            145, 195, 86, 210, 230, 219, 176, 179, 148, 236, 194, 133, 166, 240, 60, 111, 152, 154,
            62, 190, 248, 224, 197, 250, 131, 57, 145, 237, 213, 7,
        ];
        let kzg_k: [u8; 30] = [
            239, 6, 89, 163, 169, 250, 184, 159, 153, 181, 70, 47, 167, 56, 153, 92, 52, 197, 196,
            244, 10, 197, 235, 26, 46, 57, 18, 194, 56, 0,
        ];
        let a0 = BigUint::from_str(
            "3042729463975785669077695901360320813980996043134603468597671969223884",
        )
        .unwrap();
        let b0 = BigUint::from_str(
            "3099898550361810144312020021372781014260137110595159844100103797269587",
        )
        .unwrap();

        let verifier_payload: VerifierPayloadRef = VerifierPayloadRef {
            proof: ProofRef {
                commit_p,
                kzg_k,
                a0,
                b0,
            },
            public_input: RawPublicInputsRef { deposit_index },
            trapdoor: TrapdoorRef {
                tau,
                delta,
                epsilon,
            },
        };
        let witness = verifier_payload.to_bits();

        bld.show_gate_counts();

        println!("label_info {:?}", label_info);
        let passed_val = evaluate_verifier(&mut bld, witness, label_info.output_label);
        assert!(!passed_val, "verification should have failed but passed");
    }

    #[test]
    fn test_get_pub_hash_from_raw_pub_inputs() {
        let mut bld = CktBuilder::default();

        let deposit_index_labels: [usize; PUBLIC_INPUT_LEN] = bld.fresh();
        let raw_pub_in = &RawPublicInputs {
            deposit_index: deposit_index_labels,
        };
        let out_labels = get_pub_hash_from_raw_pub_inputs(&mut bld, raw_pub_in);

        bld.show_gate_counts();

        let deposit_index_raw_bytes = u64::from_le_bytes([55, 0, 0, 0, 89, 0, 0, 0]);
        let deposit_index = u64_to_bits_le(deposit_index_raw_bytes).to_vec();
        let wires = bld.eval_gates(&deposit_index);

        let out_val: BigUint = out_labels
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });

        let out_val_ref = dv_ref::get_pub_hash_from_raw_pub_inputs(&RawPublicInputsRef {
            deposit_index: u64::from_le_bytes([55, 0, 0, 0, 89, 0, 0, 0]),
        });
        assert_eq!(out_val, out_val_ref);
    }

    #[test]
    fn test_get_fs_challenge() {
        let mut bld = CktBuilder::default();

        let commit_p = {
            let r: [usize; 240] = bld.fresh();
            let r: Vec<[usize; 8]> = r
                .chunks(8)
                .map(|x| {
                    let y: [usize; 8] = x.try_into().unwrap();
                    y
                })
                .collect();
            let r: [[usize; 8]; 30] = r.try_into().unwrap();
            r
        };

        let pub0: Fr = bld.fresh();
        let pub1: Fr = bld.fresh();
        let pubs = [pub0, pub1];

        let challenge_labels = get_fs_challenge(&mut bld, commit_p, pubs, vec![], vec![]);

        let mut witness = Vec::new();
        let commit_p: [u8; 30] = [
            149, 102, 73, 129, 207, 1, 170, 225, 187, 192, 126, 126, 208, 3, 54, 148, 170, 148,
            114, 143, 39, 215, 251, 62, 10, 32, 20, 146, 207, 0,
        ];
        let mut commit_p: Vec<bool> = commit_p
            .iter()
            .flat_map(|x| u8_to_bits_le(*x).to_vec())
            .collect();
        let mut pub0 = frref_to_bits(
            &BigUint::from_str(
                "7527402554317099476086310993202889463751940730940407143885949231928",
            )
            .unwrap(),
        )
        .to_vec();
        let mut pub1 = frref_to_bits(
            &BigUint::from_str(
                "19542051593079647282099705468191403958371264520862632234952945594121",
            )
            .unwrap(),
        )
        .to_vec();

        witness.append(&mut commit_p);
        witness.append(&mut pub0);
        witness.append(&mut pub1);

        let wires = bld.eval_gates(&witness);
        let challenge_val: BigUint = challenge_labels
            .iter()
            .enumerate()
            .fold(BigUint::ZERO, |acc, (i, &w_id)| {
                acc + (BigUint::from(wires[w_id] as u16) << i)
            });

        let challenge_val_ref = {
            let commit_p = [
                149, 102, 73, 129, 207, 1, 170, 225, 187, 192, 126, 126, 208, 3, 54, 148, 170, 148,
                114, 143, 39, 215, 251, 62, 10, 32, 20, 146, 207, 0,
            ];
            let public_inputs = vec![
                BigUint::from_str(
                    "7527402554317099476086310993202889463751940730940407143885949231928",
                )
                .unwrap(),
                BigUint::from_str(
                    "19542051593079647282099705468191403958371264520862632234952945594121",
                )
                .unwrap(),
            ];
            dv_ref::get_fs_challenge(commit_p, public_inputs, vec![], vec![])
        };

        assert_eq!(challenge_val, challenge_val_ref);
    }
}
