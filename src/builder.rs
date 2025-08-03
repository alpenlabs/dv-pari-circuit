use std::collections::HashMap;

use mcircuit::Operation;
use mcircuit::exporters::bool_circuit_to_bristol;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelExtend,
    ParallelIterator,
};

#[derive(Debug, Clone)]
pub struct CustomGateParams {
    gate_type: CustomGateType,
    input_wire_labels: Vec<usize>,
    diff_start_wire: usize,
}

#[derive(Debug, Clone)]
pub enum GateOperation {
    Base(Operation<bool>),
    Custom(CustomGateParams),
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum CustomGateType {
    PointAdd,
    // add more type here e.g fr_mul when needed for memory balance
}

pub(crate) trait Circuit {
    fn fresh_one(&mut self) -> usize;

    fn fresh<const N: usize>(&mut self) -> [usize; N];

    fn zero(&mut self) -> usize;

    fn one(&mut self) -> usize;

    fn xor_wire(&mut self, x: usize, y: usize) -> usize;

    fn or_wire(&mut self, x: usize, y: usize) -> usize;

    fn and_wire(&mut self, x: usize, y: usize) -> usize;

    fn push_custom_gate(&mut self, params: CustomGateParams, new_wire_idx: usize);

    fn get_gates(&self) -> &Vec<GateOperation>;

    fn next_wire(&self) -> usize;

    fn push_template_for_gate(&mut self, templ_type: CustomGateType) -> &Template;

    fn get_template(&self, templ_type: CustomGateType) -> Option<&Template>;
}

impl Circuit for CktBuilder {
    fn fresh_one(&mut self) -> usize {
        let w = self.next_wire;
        self.next_wire += 1;
        w
    }

    fn fresh<const N: usize>(&mut self) -> [usize; N] {
        let mut out = [0; N];
        for slot in &mut out {
            *slot = self.fresh_one();
        }
        out
    }

    fn zero(&mut self) -> usize {
        /* constant “0” (one per circuit) */
        if let Some(z) = self.zero {
            return z;
        }
        let w = self.fresh_one();
        self.zero = Some(w);
        w
    }

    fn one(&mut self) -> usize {
        /* constant “1” (one per circuit) */
        if let Some(z) = self.one {
            return z;
        }
        let w = self.fresh_one();
        self.one = Some(w);
        w
    }

    fn xor_wire(&mut self, x: usize, y: usize) -> usize {
        let w = self.fresh_one();
        self.gates
            .push(GateOperation::Base(Operation::Add(w, x, y)));
        w
    }

    fn or_wire(&mut self, x: usize, y: usize) -> usize {
        let xpy = self.fresh_one();
        let xmy = self.fresh_one();
        let z = self.fresh_one();
        self.gates
            .push(GateOperation::Base(Operation::Add(xpy, x, y)));
        self.gates
            .push(GateOperation::Base(Operation::Mul(xmy, x, y)));
        self.gates
            .push(GateOperation::Base(Operation::Add(z, xpy, xmy)));
        z
    }

    fn and_wire(&mut self, x: usize, y: usize) -> usize {
        let w = self.fresh_one();
        self.gates
            .push(GateOperation::Base(Operation::Mul(w, x, y)));
        w
    }

    fn push_custom_gate(&mut self, params: CustomGateParams, new_wire_idx: usize) {
        self.gates.push(GateOperation::Custom(params));
        self.next_wire = new_wire_idx;
    }

    fn get_gates(&self) -> &Vec<GateOperation> {
        &self.gates
    }

    fn next_wire(&self) -> usize {
        self.next_wire
    }

    fn push_template_for_gate(&mut self, templ_type: CustomGateType) -> &Template {
        match templ_type {
            CustomGateType::PointAdd => {
                if self.templates.ptadd_template.is_none() {
                    let add_templ = template_emit_point_add();
                    self.templates.ptadd_template = Some(add_templ);
                }
            }
        }
        self.templates.ptadd_template.as_ref().unwrap()
    }

    fn get_template(&self, templ_type: CustomGateType) -> Option<&Template> {
        match templ_type {
            CustomGateType::PointAdd => self.templates.ptadd_template.as_ref(),
        }
    }
}

#[derive(Default, Debug)]
pub struct Templates {
    pub(crate) ptadd_template: Option<Template>,
}

/// CktBuilder
#[derive(Debug)]
pub struct CktBuilder {
    pub gates: Vec<GateOperation>,
    pub next_wire: usize,
    pub zero: Option<usize>,
    pub one: Option<usize>,
    pub templates: Templates,
}

impl Default for CktBuilder {
    fn default() -> Self {
        CktBuilder {
            gates: vec![],
            next_wire: 2,
            zero: Some(0),
            one: Some(1),
            templates: Templates {
                ptadd_template: None,
            },
        }
    }
}

pub(crate) fn xor_three<T: Circuit>(b: &mut T, x: usize, y: usize, z: usize) -> usize {
    let x_xor_y = b.xor_wire(x, y);
    b.xor_wire(x_xor_y, z)
}

/* ───────────── 1 · tiny “template” wrapper (inputs only, auto-mapping) ─── */

pub(crate) fn xor_vec<T: Circuit>(b: &mut T, a: &[usize], b_: &[usize]) -> Vec<usize> {
    let len = a.len().max(b_.len());
    (0..len)
        .map(|i| match (a.get(i), b_.get(i)) {
            (Some(&x), Some(&y)) => b.xor_wire(x, y),
            (Some(&x), None) => x,
            (None, Some(&y)) => y,
            _ => unreachable!(),
        })
        .collect()
}

pub(crate) fn add_shifted<T: Circuit>(
    b: &mut T,
    dst: &mut [Option<usize>],
    shift: usize,
    src: &[usize],
) {
    for (i, &w) in src.iter().enumerate() {
        let idx = i + shift;
        dst[idx] = match dst[idx] {
            None => Some(w),
            Some(prev) => Some(b.xor_wire(prev, w)),
        };
    }
}
pub(crate) fn finish_bits<T: Circuit>(b: &mut T, bits: Vec<Option<usize>>) -> Vec<usize> {
    bits.into_iter()
        .map(|opt| opt.unwrap_or_else(|| b.zero()))
        .collect()
}

/* ────────────────────  GF(2⁹) modular reduction  ────────────────────── */
/*  polynomial:  m(x) = x⁹ + x⁴ + 1    => 17-bit → 9-bit linear map       */

/* XOR helper for an arbitrary list */
pub(crate) fn xor_many<T: Circuit>(b: &mut T, items: impl IntoIterator<Item = usize>) -> usize {
    let mut it = items.into_iter();
    let first = match it.next() {
        Some(w) => w,
        None => return b.zero(),
    };
    it.fold(first, |acc, w| b.xor_wire(acc, w))
}

#[derive(Default, Clone, Debug)]
pub(crate) struct Template {
    pub input_wires: Vec<usize>,
    pub gates: Vec<GateOperation>,
    pub output_wires: Vec<usize>,
    pub start_wire_idx: usize,
    pub end_wire_idx: usize,
    pub const_wire_zero: usize,
    pub const_wire_one: usize,
    pub stats: (usize, usize),
}

impl Template {
    pub(crate) fn emit_point_add_custom<T: Circuit>(
        bld: &mut T,
        p1: &CurvePoint,
        p2: &CurvePoint,
    ) -> CurvePoint {
        let mut input_wires = vec![];
        input_wires.extend_from_slice(&p1.x);
        input_wires.extend_from_slice(&p1.s);
        input_wires.extend_from_slice(&p1.z);
        input_wires.extend_from_slice(&p1.t);

        input_wires.extend_from_slice(&p2.x);
        input_wires.extend_from_slice(&p2.s);
        input_wires.extend_from_slice(&p2.z);
        input_wires.extend_from_slice(&p2.t);

        let output_wires = Self::emit_custom(bld, input_wires, CustomGateType::PointAdd);

        CurvePoint {
            x: output_wires[0..GF_LEN].try_into().unwrap(),
            s: output_wires[GF_LEN..GF_LEN * 2].try_into().unwrap(),
            z: output_wires[GF_LEN * 2..GF_LEN * 3].try_into().unwrap(),
            t: output_wires[GF_LEN * 3..GF_LEN * 4].try_into().unwrap(),
        }
    }

    fn emit_custom<T: Circuit>(
        bld: &mut T,
        input_wires: Vec<usize>,
        gate_type: CustomGateType,
    ) -> Vec<usize> {
        if bld.get_template(gate_type).is_none() {
            bld.push_template_for_gate(gate_type);
        }
        let tmpl = bld.get_template(gate_type).unwrap();

        assert_eq!(tmpl.input_wires.len(), input_wires.len());

        let diff_start_wire = bld.next_wire() - tmpl.start_wire_idx;

        let ref_output_set: Vec<usize> = tmpl
            .output_wires
            .iter()
            .map(|x| x + diff_start_wire)
            .collect();

        let templ_params = CustomGateParams {
            gate_type,
            input_wire_labels: input_wires,
            diff_start_wire,
        };
        let next_wire = bld.next_wire();
        bld.push_custom_gate(
            templ_params,
            next_wire + (tmpl.end_wire_idx - tmpl.start_wire_idx),
        );

        ref_output_set
    }

    pub(crate) fn unroll_custom_gate(
        &self,
        zero: Option<usize>,
        one: Option<usize>,
        diff_start_wire: usize,
        input_wires: &[usize],
    ) -> Vec<Operation<bool>> {
        assert_eq!(self.input_wires.len(), input_wires.len());

        let mut wire_map: HashMap<usize, usize> = self
            .input_wires
            .par_iter()
            .copied()
            .zip(input_wires.to_owned())
            .collect();
        wire_map.par_extend(
            (self.start_wire_idx..self.end_wire_idx)
                .into_par_iter()
                .map(|old_idx| (old_idx, old_idx + diff_start_wire)),
        );
        if let Some(zero) = zero {
            wire_map.insert(self.const_wire_zero, zero);
        }
        if let Some(one) = one {
            wire_map.insert(self.const_wire_one, one);
        }

        let ops: Vec<Operation<bool>> = self
            .gates
            .par_iter()
            .map(|h| match h {
                GateOperation::Base(g) => {
                    let ret = match *g {
                        Operation::Add(d, x, y) => {
                            let nd = wire_map.get(&d).unwrap();
                            let nx = wire_map.get(&x).unwrap();
                            let ny = wire_map.get(&y).unwrap();
                            Operation::Add(*nd, *nx, *ny)
                        }
                        Operation::Mul(d, x, y) => {
                            let nd = wire_map.get(&d).unwrap();
                            let nx = wire_map.get(&x).unwrap();
                            let ny = wire_map.get(&y).unwrap();
                            Operation::Mul(*nd, *nx, *ny)
                        }
                        Operation::Const(d, v) => {
                            let nd = wire_map.get(&d).unwrap();
                            Operation::Const(*nd, v)
                        }
                        _ => unreachable!(),
                    };
                    vec![ret]
                }
                GateOperation::Custom(g) => {
                    self.unroll_custom_gate(zero, one, g.diff_start_wire, &g.input_wire_labels)
                }
            })
            .flatten()
            .collect();

        ops
    }
}

use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::curve_ckt::{CurvePoint, template_emit_point_add};
use crate::gf_ckt::GF_LEN;

impl CktBuilder {
    pub fn write_bristol_periodic(
        &mut self,
        // bool_witness: &[bool],
        path: &str,
    ) -> io::Result<()> {
        const CHUNK_GATES: usize = 10_000; // how many gates to render at once
        const FLUSH_EVERY: usize = 4 << 20; // flush after ~4 MiB have been written
        let mut writer = BufWriter::with_capacity(4 << 20, File::create(path)?);

        let mut bytes_since_flush = 0;

        let gates = &self.gates;
        let zero_gate = self.zero;
        let one_gate = self.one;

        let mut flushable_gates = vec![];
        for i in 0..gates.len() {
            let gate = &gates[i];
            let mut inner_gates = match gate {
                GateOperation::Base(g) => {
                    vec![*g]
                }
                GateOperation::Custom(params) => {
                    let custom_gate_template = match params.gate_type {
                        CustomGateType::PointAdd => {
                            self.templates.ptadd_template.as_ref().unwrap() // assumed it exists
                        }
                    };

                    custom_gate_template.unroll_custom_gate(
                        zero_gate,
                        one_gate,
                        params.diff_start_wire,
                        &params.input_wire_labels,
                    )
                }
            };

            flushable_gates.append(&mut inner_gates);
            if flushable_gates.len() > CHUNK_GATES || i == gates.len() - 1 {
                for chunk in flushable_gates.chunks(CHUNK_GATES) {
                    // Render the whole chunk at once with the original helper.
                    let rendered: String = bool_circuit_to_bristol(chunk, &[]);

                    // Stream it to the file.
                    writer.write_all(rendered.as_bytes())?;
                    bytes_since_flush += rendered.len();

                    // Force a flush every FLUSH_EVERY bytes so data hits the OS regularly.
                    if bytes_since_flush >= FLUSH_EVERY {
                        writer.flush()?;
                        bytes_since_flush = 0;
                    }
                }
                flushable_gates.clear();
            }
        }

        writer.flush()?; // final flush
        Ok(())
    }

    pub fn show_gate_counts(&mut self) {
        let mut and_gates_count = 0;
        let mut xor_gates_count = 0;
        let mut const_gates_count = 0;
        let mut custom_gates_count = 0;

        let (tmp_mul, tmp_add) = {
            if let Some(tmpl) = self.get_template(CustomGateType::PointAdd) {
                tmpl.stats
            } else {
                (0, 0)
            }
        };

        for h in self.get_gates() {
            match h {
                GateOperation::Base(g) => match g {
                    Operation::Add(_, _, _) => {
                        xor_gates_count += 1;
                    }
                    Operation::Mul(_, _, _) => {
                        and_gates_count += 1;
                    }
                    Operation::Const(_, _) => {
                        const_gates_count += 1;
                    }
                    _ => unreachable!(),
                },
                GateOperation::Custom(_) => {
                    custom_gates_count += 1;
                    xor_gates_count += tmp_add;
                    and_gates_count += tmp_mul;
                }
            }
        }
        println!("and_gates_count {}", and_gates_count);
        println!("xor_gates_count {}", xor_gates_count);
        println!("const_gates_count {}", const_gates_count);
        println!("custom_gates_count {}", custom_gates_count);
        println!("total {}", self.get_gates().len());
    }

    pub fn eval_gates(&self, witness: &[bool]) -> Vec<bool> {
        let n_wires = self.next_wire;
        let mut w = vec![None; n_wires];
        let gates = &self.gates;
        let zero_gate = self.zero;
        let one_gate = self.one;
        assert_eq!(zero_gate, Some(0));
        w[zero_gate.unwrap()] = Some(false);
        assert_eq!(one_gate, Some(1));
        w[one_gate.unwrap()] = Some(true);

        const WIRE_OFFSET: usize = 2; // due to 0 and 1 at the first two index
        for (id, &bit) in witness.iter().enumerate() {
            w[WIRE_OFFSET + id] = Some(bit);
        }

        for h in gates {
            match h {
                GateOperation::Base(g) => {
                    match *g {
                        Operation::Add(d, x, y) => w[d] = Some(w[x].unwrap() ^ w[y].unwrap()),
                        Operation::Mul(d, x, y) => w[d] = Some(w[x].unwrap() & w[y].unwrap()),
                        Operation::Const(d, v) => w[d] = Some(v),
                        _ => unreachable!(), // no other variants used
                    }
                }
                GateOperation::Custom(params) => {
                    let custom_gate_template = match params.gate_type {
                        CustomGateType::PointAdd => self.templates.ptadd_template.as_ref().unwrap(),
                    };
                    let gates = custom_gate_template.unroll_custom_gate(
                        zero_gate,
                        one_gate,
                        params.diff_start_wire,
                        &params.input_wire_labels,
                    );
                    for g in gates {
                        match g {
                            Operation::Add(d, x, y) => w[d] = Some(w[x].unwrap() ^ w[y].unwrap()),
                            Operation::Mul(d, x, y) => w[d] = Some(w[x].unwrap() & w[y].unwrap()),
                            Operation::Const(d, v) => w[d] = Some(v),
                            _ => unreachable!(), // no other variants used
                        }
                    }
                }
            }
        }

        w.iter().map(|x| x.unwrap()).collect()
    }
}
