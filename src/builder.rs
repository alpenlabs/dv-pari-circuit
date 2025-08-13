//! Module provides a way to build binary circuits
use std::collections::HashMap;

use mcircuit::Operation;
use mcircuit::exporters::bool_circuit_to_bristol;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelExtend,
    ParallelIterator,
};

/// Boolean Operation
#[derive(Debug, Clone)]
pub enum GateOperation {
    /// Basic Boolean Operations: AND, XOR
    Base(Operation<bool>),
    /// Custom Boolean Circuit: Point Add, Scalar Mul, etc
    Custom(CustomGateParams),
}

/// Custom Gate Type
#[derive(Debug, Clone, Copy)]
pub(crate) enum CustomGateType {
    PointAdd,
    // add more type here e.g fr_mul when needed for memory balance
}

/// Parameters that specify an instance of boolean circuit
///
/// An instance is uniquely identified by its type (which dictates circuit configuration),
/// what wires it is connected to (i.e. input wires to this circuit) and
/// a unique assignment for its internal wires. These are given by `gate_type`, `input_wire_labels`
/// and `internal_wire_label_start_index`
#[derive(Debug, Clone)]
pub struct CustomGateParams {
    /// type of custom gate: e.g Point Add
    gate_type: CustomGateType,
    /// wire labels that the custom gate takes as input
    input_wire_labels: Vec<usize>,
    /// offset for internal wire labels
    // We need to assign unique labels to internal wire
    // `internal_wire_label_start_index` is the starting index
    // all internal wires and output wires are increments from this value
    internal_wire_label_start_index: usize,
}

/// Circuit Trait specifies how you represent the entire binary circuit.
/// It provides a way to prepare wire labels (input and constant gates),
/// assemble different type of gates (boolean or custom) and inspect them.
///
/// How wire labels are assigned and how the gates are represented in memory
/// depends upon the implementation.
///
/// Wire lables are represented by 'usize' right now.
/// TODO: do not tie the trait with concrete type 'usize',
/// this way you can use u32 to reduce memory usage if there are less wire counts.
pub(crate) trait Circuit {
    /// get fresh wire label
    fn fresh_one(&mut self) -> usize;

    /// get fresh wire labels
    fn fresh<const N: usize>(&mut self) -> [usize; N];

    /// get wire label for constant zero wire
    fn zero(&mut self) -> usize;

    /// get wire label for constant one wire
    fn one(&mut self) -> usize;

    /// XOR two input wires and return wire label of output
    fn xor_wire(&mut self, x: usize, y: usize) -> usize;

    /// OR two input wires and return wire label of output
    fn or_wire(&mut self, x: usize, y: usize) -> usize;

    // AND two input wires and return wire label of output
    fn and_wire(&mut self, x: usize, y: usize) -> usize;

    /// push an instance of custom gate. `params` specifies the current instance of Custom Gate.
    /// `new_wire_idx` represents the new next wire label
    fn push_custom_gate(&mut self, params: CustomGateParams, new_wire_idx: usize);

    /// get all gates
    fn get_gates(&self) -> &Vec<GateOperation>;

    /// next wire label
    fn next_wire(&self) -> usize;

    /// initialize cicruit configuration for custom gate
    /// `Template` specifies circuit configuration.
    // This is known only by executing the operation represented by circuit type.
    // For example to determine circuit configuration of a Point Addition,
    // you execute circuit compiler for that corresponding module
    fn init_circuit_config_for_custom_gate(&mut self, templ_type: CustomGateType) -> &Template;

    /// get template (circuit configuration) of `templ_type` if it exists
    fn get_template(&self, templ_type: CustomGateType) -> Option<&Template>;
}

/// CktBuilder: implementation of Circuit Trait.
/// Implements method that walks through the list of gates to
/// export them, evaluate circuit or obtain statistics (gate types and counts)
///
/// Wire labels are unique integer values and incremented by one each time.
/// The first two wire labels are hardcoded to const wires 0 and 1.
/// All input wires should then be initialized by a call to CktBuilder::fresh().
/// Successive wires are internal wires.
/// Relevant function that compiles circuit will return output wire labels for that task.
#[derive(Debug)]
pub struct CktBuilder {
    // all gates
    pub gates: Vec<GateOperation>,
    // next wire label
    pub next_wire: usize,
    // wire label of const zero wire
    pub zero: Option<usize>,
    // wire label of const one wire
    pub one: Option<usize>,
    // custom gate configurations
    pub templates: Templates,
}

#[derive(Default, Debug)]
pub struct Templates {
    // circuit configuration for point addition
    pub(crate) ptadd_template: Option<Template>,
}

impl Default for CktBuilder {
    fn default() -> Self {
        CktBuilder {
            gates: vec![],
            next_wire: 2,
            zero: Some(0),
            one: Some(1),
            templates: Templates {
                // initialized only when it's needed because
                // not all circuits may require this template
                ptadd_template: None,
            },
        }
    }
}

impl Circuit for CktBuilder {
    fn fresh_one(&mut self) -> usize {
        let w = self.next_wire;
        self.next_wire += 1;
        // increment wire label by 1 to have new wire label be unique
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

    fn init_circuit_config_for_custom_gate(&mut self, templ_type: CustomGateType) -> &Template {
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

// few helper modules
pub(crate) fn xor_three<T: Circuit>(b: &mut T, x: usize, y: usize, z: usize) -> usize {
    let x_xor_y = b.xor_wire(x, y);
    b.xor_wire(x_xor_y, z)
}

pub(crate) fn xor_vec<T: Circuit>(bld: &mut T, a: &[usize], b: &[usize]) -> Vec<usize> {
    let len = a.len().max(b.len());
    (0..len)
        .map(|i| match (a.get(i), b.get(i)) {
            (Some(&x), Some(&y)) => bld.xor_wire(x, y),
            (Some(&x), None) => x,
            (None, Some(&y)) => y,
            _ => unreachable!(),
        })
        .collect()
}

// XOR helper for an arbitrary list
pub(crate) fn xor_many<T: Circuit>(bld: &mut T, items: impl IntoIterator<Item = usize>) -> usize {
    let mut it = items.into_iter();
    let first = match it.next() {
        Some(w) => w,
        None => return bld.zero(),
    };
    it.fold(first, |acc, w| bld.xor_wire(acc, w))
}

/// Custom Gate Configuration
///
/// Configuration of a circuit is defined by its input wires, gates and output wires.
/// We include `start_wire_idx` and `end_wire_idx` so that a different instance of
/// a binary circuit of a configuration can be instantiated. These fields help uniquely
/// label internal wires, while the other fields specify how these wire labels are connected.
///
/// The `stats` field is a metric to obtain gate counts for benchmarks
#[derive(Default, Clone, Debug)]
pub(crate) struct Template {
    /// input wire labels to this circuit
    pub input_wires: Vec<usize>,
    /// logic gates in this circuit
    pub gates: Vec<GateOperation>,
    /// output wire labels from this circuit
    pub output_wires: Vec<usize>,
    /// wire labels corresponding to constant value zero
    pub const_wire_zero: usize,
    // wire labels corresponding to constant value one
    pub const_wire_one: usize,

    /// starting wire label of internal wires
    pub start_wire_idx: usize,
    /// final wire label of internal wires plus one
    // "plus one" here because end_wire_idx saves value of CktBuilder::next_wire
    pub end_wire_idx: usize,

    // count of AND, XOR gates in this circuit
    pub stats: (usize, usize),
}

impl Template {
    /// Generate binary circuit for point addition using cached configuration 'template' for this circuit.
    /// Input is the same as it would normally be for point addition, which is two CurvePoints
    pub(crate) fn emit_point_add_custom<T: Circuit>(
        bld: &mut T,
        p1: &CurvePoint,
        p2: &CurvePoint,
    ) -> CurvePoint {
        // serialize wire labels in the same order accepted by the circuit configuration
        let mut input_wires = vec![];
        input_wires.extend_from_slice(&p1.x);
        input_wires.extend_from_slice(&p1.s);
        input_wires.extend_from_slice(&p1.z);
        input_wires.extend_from_slice(&p1.t);

        input_wires.extend_from_slice(&p2.x);
        input_wires.extend_from_slice(&p2.s);
        input_wires.extend_from_slice(&p2.z);
        input_wires.extend_from_slice(&p2.t);

        // Generate an instance of this custom gate using circuit configuration for the template
        // and return output wire labels
        let output_wires = Self::emit_custom(bld, input_wires, CustomGateType::PointAdd);

        // Deserialize wire labels into expected data structure
        CurvePoint {
            x: output_wires[0..GF_LEN].try_into().unwrap(),
            s: output_wires[GF_LEN..GF_LEN * 2].try_into().unwrap(),
            z: output_wires[GF_LEN * 2..GF_LEN * 3].try_into().unwrap(),
            t: output_wires[GF_LEN * 3..GF_LEN * 4].try_into().unwrap(),
        }
    }

    /// Instantiate binary circuit of known configuration.
    /// Configuration is known by `gate_type`. Each instance of the circuit has its own set of inputs,
    /// which is defined by `input_wires`
    fn emit_custom<T: Circuit>(
        bld: &mut T,
        input_wires: Vec<usize>,
        gate_type: CustomGateType,
    ) -> Vec<usize> {
        // Build Configuration if it doesn't already exist
        if bld.get_template(gate_type).is_none() {
            bld.init_circuit_config_for_custom_gate(gate_type);
        }
        let tmpl = bld.get_template(gate_type).unwrap();

        assert_eq!(tmpl.input_wires.len(), input_wires.len());

        let internal_wire_starting_index = bld.next_wire() - tmpl.start_wire_idx;

        // output wires are also labelled as any internal wire i.e. assigned a unique value
        let ref_output_set: Vec<usize> = tmpl
            .output_wires
            .iter()
            .map(|x| x + internal_wire_starting_index)
            .collect();

        let next_wire = bld.next_wire();
        bld.push_custom_gate(
            CustomGateParams {
                gate_type,
                input_wire_labels: input_wires,
                internal_wire_label_start_index: internal_wire_starting_index,
            },
            // new wire label is obtained by adding the next wire label by
            // the total number of internal wires defined within this configuration.
            // like any batch wire assignment
            next_wire + (tmpl.end_wire_idx - tmpl.start_wire_idx),
        );

        ref_output_set
    }

    /// Convert custom gate to basic gates
    pub(crate) fn unroll_custom_gate(
        &self,
        zero: Option<usize>,
        one: Option<usize>,
        internal_wire_start_offset: usize,
        input_wires: &[usize],
    ) -> Vec<Operation<bool>> {
        assert_eq!(self.input_wires.len(), input_wires.len());

        // wire_map maps template wire labels to instance specific wire labels
        let mut wire_map: HashMap<usize, usize> = self
            .input_wires
            .par_iter()
            .copied()
            .zip(input_wires.to_owned())
            .collect();

        // iterate through internal wire labels of the template
        // and offset by a value `internal_wire_start_offset` specific to this instance to generate new wire labels
        wire_map.par_extend(
            (self.start_wire_idx..self.end_wire_idx)
                .into_par_iter()
                .map(|old_idx| (old_idx, old_idx + internal_wire_start_offset)),
        );
        if let Some(zero) = zero {
            wire_map.insert(self.const_wire_zero, zero);
        }
        if let Some(one) = one {
            wire_map.insert(self.const_wire_one, one);
        }

        // Iterate through each of the gates and map wire labels from the ones in template
        // to wire lables unique to this instant
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
                GateOperation::Custom(g) => self.unroll_custom_gate(
                    zero,
                    one,
                    g.internal_wire_label_start_index,
                    &g.input_wire_labels,
                ),
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
    /// Dump gates to bristol format
    // Iterate through gates, unroll custom gates when needed,
    // and write them to file in bristol format
    pub fn write_bristol_periodic(&mut self, path: &str) -> io::Result<()> {
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
                        params.internal_wire_label_start_index,
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

    /// show gate counts in the current CktBuilder
    pub fn show_gate_counts(&mut self) {
        let mut and_gates_count = 0;
        let mut xor_gates_count = 0;
        let mut custom_gates_count = 0;

        let (tmp_mul, tmp_add) = {
            if let Some(tmpl) = self.get_template(CustomGateType::PointAdd) {
                (Some(tmpl.stats.0), Some(tmpl.stats.1))
            } else {
                // if custom gate is not initialized return 0
                (None, None)
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
                    _ => unreachable!(),
                },
                GateOperation::Custom(_) => {
                    custom_gates_count += 1;
                    if let Some(tmp_add) = tmp_add {
                        xor_gates_count += tmp_add;
                    }
                    if let Some(tmp_mul) = tmp_mul {
                        and_gates_count += tmp_mul;
                    }
                }
            }
        }
        println!("and_gates_count {}", and_gates_count);
        println!("xor_gates_count {}", xor_gates_count);
        println!("custom_gates_count {}", custom_gates_count);
        println!("total_gates_count {}", self.get_gates().len());
    }

    /// Evaluate a binary circuit given `witness` as input wire values.
    /// Assumes that these `witness` values correspond to wire labels from 2 to 2+num_input_wires.
    /// as the first two wire labels are always constant wire labels 0 and 1.
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
                        params.internal_wire_label_start_index,
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
