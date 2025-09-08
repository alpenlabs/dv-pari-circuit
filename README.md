# DV-Pari Circuit

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache-blue.svg)](https://opensource.org/licenses/apache-2-0)
[![ci](https://github.com/alpenlabs/rust-template/actions/workflows/lint.yml/badge.svg?event=push)](https://github.com/alpenlabs/rust-template/actions)
[![docs](https://img.shields.io/badge/docs-docs.rs-orange)](https://docs.rs/rust-template)

This repo includes binary circuit representation of DV-Pari verifier program designed for low number of AND gates -- currently around 11 million.

A verifier program is built from the following functions represented as circuits:
- Base Field Multiplication
- Scalar Field Multiplication
- Curve Operations
- Blake3 Hash Operation

User operations:
- At setup time, compile the circuit and export it to a format (like bristol)
- At runtime, evaluate the circuit for a given bit stream

Compilation and evaluation currently takes around 13 minutes with 9 GB of peak memory.

Check out [docs](https://github.com/alpenlabs/dv-pari-circuit/tree/main/docs) for more.
