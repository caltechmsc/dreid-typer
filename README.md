# DreidTyper

**DreidTyper** is a high-performance, foundational software library for the automated assignment of DREIDING force field atom types and the perception of molecular topologies. It provides a modern, robust solution for translating simple chemical connectivity (a molecular graph) into a complete, engine-agnostic topological description essential for molecular simulations. This library is engineered from the ground up in **Rust** for exceptional performance, memory safety, and strict adherence to the principles of modular software design.

The core mission of DreidTyper is to provide a reliable, predictable, and easy-to-integrate tool for developers and researchers building the next generation of simulation tools for general chemistry, materials science, and drug discovery.

## Features

- **DREIDING Atom Typing**: Assigns canonical DREIDING atom types from molecular connectivity.
- **Full Topology Perception**: Identifies bonds, angles, and proper/improper dihedrals.
- **Memory Safe & Fast**: Built in Rust for guaranteed memory safety and high performance.
- **Rule-Based Engine**: Uses a clear TOML-based rule system for atom typing logic.
- **Engine-Agnostic**: Produces a pure topological representation independent of any MD engine.

## Getting Started

To get started with DreidTyper, add it as a dependency in your `Cargo.toml`:

```toml
[dependencies]
dreid-typer = "0.2.0"
```

Then, you can use it in your Rust code as follows:

```rust
use dreid_typer::{
    assign_topology, MolecularGraph, MolecularTopology,
    Element, BondOrder,
};

// 1. Define the molecule's connectivity using a `MolecularGraph`.
let mut graph = MolecularGraph::new();
let c1 = graph.add_atom(Element::C); // CH3
let c2 = graph.add_atom(Element::C); // CH2
let o = graph.add_atom(Element::O);
let h_c1_1 = graph.add_atom(Element::H);
let h_c1_2 = graph.add_atom(Element::H);
let h_c1_3 = graph.add_atom(Element::H);
let h_c2_1 = graph.add_atom(Element::H);
let h_c2_2 = graph.add_atom(Element::H);
let h_o = graph.add_atom(Element::H);

graph.add_bond(c1, c2, BondOrder::Single).unwrap();
graph.add_bond(c2, o, BondOrder::Single).unwrap();
graph.add_bond(c1, h_c1_1, BondOrder::Single).unwrap();
graph.add_bond(c1, h_c1_2, BondOrder::Single).unwrap();
graph.add_bond(c1, h_c1_3, BondOrder::Single).unwrap();
graph.add_bond(c2, h_c2_1, BondOrder::Single).unwrap();
graph.add_bond(c2, h_c2_2, BondOrder::Single).unwrap();
graph.add_bond(o, h_o, BondOrder::Single).unwrap();

// 2. Call the main function to perceive the topology.
let topology: MolecularTopology = assign_topology(&graph).unwrap();

// 3. Inspect the results.
assert_eq!(topology.atoms.len(), 9);
assert_eq!(topology.bonds.len(), 8);
assert_eq!(topology.angles.len(), 13);
assert_eq!(topology.proper_dihedrals.len(), 12);

// Check the assigned DREIDING atom types.
assert_eq!(topology.atoms[c1].atom_type, "C_3");   // sp3 Carbon
assert_eq!(topology.atoms[c2].atom_type, "C_3");   // sp3 Carbon
assert_eq!(topology.atoms[o].atom_type, "O_3");    // sp3 Oxygen
assert_eq!(topology.atoms[h_o].atom_type, "H_HB"); // Hydrogen-bonding Hydrogen
assert_eq!(topology.atoms[h_c1_1].atom_type, "H_"); // Standard Hydrogen
```

> **Note**: This is a simplified example. For more complex molecules and edge cases, please refer to the [API Documentation](https://docs.rs/dreid-typer).

## Documentation

- [API Documentation](https://docs.rs/dreid-typer) - Comprehensive reference for all public types and functions.
- [Architecture Documents](docs/ARCHITECTURE.md) - In-depth design and implementation details.

## Tech Stack

- **Core Language**: Rust
- **Build System**: Cargo
- **Rule Format**: TOML

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
