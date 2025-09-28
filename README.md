# DreidTyper

**DreidTyper** is a high-performance, foundational software library for the automated assignment of DREIDING force field atom types and the perception of molecular topologies. It provides a modern, robust solution for translating simple chemical connectivity (a molecular graph) into a complete, engine-agnostic topological description essential for molecular simulations. This library is engineered from the ground up in **Rust** for exceptional performance, memory safety, and strict adherence to the principles of modular software design.

The core mission of DreidTyper is to provide a reliable, predictable, and easy-to-integrate tool for developers and researchers building the next generation of simulation tools for general chemistry, materials science, and drug discovery.

## Features

- **DREIDING Atom Typing**: Assigns canonical DREIDING atom types from molecular connectivity.
- **Full Topology Perception**: Identifies bonds, angles, and proper/improper dihedrals.
- **Memory Safe & Fast**: Built in Rust for guaranteed memory safety and high performance.
- **Rule-Based Engine**: Uses a clear TOML-based rule system for atom typing logic.
- **Engine-Agnostic**: Produces a pure topological representation independent of any MD engine.

## Tech Stack

- **Core Language**: Rust
- **Build System**: Cargo
- **Rule Format**: TOML

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
