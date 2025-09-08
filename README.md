# DV-Pari Circuit

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache-blue.svg)](https://opensource.org/licenses/apache-2-0)
[![ci](https://github.com/alpenlabs/dv-pari-circuit/actions/workflows/lint.yml/badge.svg?event=push)](https://github.com/alpenlabs/dv-pari-circuit/actions)

This repo includes binary circuit representation of DV-Pari verifier program designed for low number of AND gates -- currently around 11 million.

A verifier program is built from the following functions represented as circuits:

- Base Field Multiplication
- Scalar Field Multiplication
- Curve Operations
- Blake3 Hash Operation

User operations:

- At setup time, compile the circuit and export it to a format (like bristol)
- At runtime, evaluate the circuit for a given bit stream

To generate the circuit, run the following:

```bash
cargo test --release --package dv-pari-circuit --lib -- dv_ckt::test::test_verify_over_mock_inputs --exact --show-output --ignored
```

Compilation and evaluation currently takes around 13 minutes with 9 GB of peak memory.

## Details of the Circuit

- SHA256 hash: 2b290772ddce07a2503571df147ef88755b6db417c59197b5b8846248ffaab7a
- Total Gates: 3284035290
- AND Gates: 11337170
- XOR Gates: 3272698120
- Total wires: 3286566319
- Primary inputs: 1706
- Intermediate wires: 3283320345
- Primary outputs: 714945
- Important Output Wire: 3284036995
- Missing/unused wires: 0

Check out [docs](docs/) for more.
