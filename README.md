# DreidTyper

**DreidTyper** is a Rust library that turns a minimal `MolecularGraph` (atoms + bonds) into a fully typed, DREIDING-compatible topology. The pipeline is deterministic, aggressively validated, and designed for integrators who need trustworthy chemistry without shipping their own perception code.

At a high level the library walks through:

1. **Perception:** six ordered passes (rings → Kekulé expansion → electron bookkeeping → aromaticity → resonance → hybridization) that upgrade raw connectivity into a rich `AnnotatedMolecule`.
2. **Typing:** an iterative, priority-sorted rule engine that resolves the final DREIDING atom label for every atom.
3. **Building:** a pure graph traversal that emits canonical bonds, angles, torsions, and inversions as a `MolecularTopology`.

## Features

- **Chemically faithful perception:** built-in algorithms cover SSSR ring search, strict Kekulé expansion, charge/lone pair templates for heteroatoms, aromaticity categorization (including anti-aromatic detection), resonance propagation, and hybridization inference.
- **Deterministic typing engine:** TOML rules are sorted by priority and evaluated until a fixed point, making neighbor-dependent rules (e.g., `H_HB`) converge without guesswork.
- **Engine-agnostic topology:** outputs canonicalized bonds, angles, torsions, and inversions ready for any simulator that consumes DREIDING-style terms.
- **Extensible ruleset:** ship with curated defaults (`resources/default.rules.toml`) and load or merge custom rule files at runtime.
- **Rust-first ergonomics:** zero `unsafe`, comprehensive unit/integration tests, and precise error variants for validation, perception, and typing failures.

## Getting Started

Add the crate to your `Cargo.toml`:

```toml
[dependencies]
dreid-typer = "0.5.0"
```

Run the full pipeline from connectivity to topology:

```rust
use dreid_typer::{assign_topology, Element, GraphBondOrder, MolecularGraph, MolecularTopology};

let mut graph = MolecularGraph::new();
let c1 = graph.add_atom(Element::C);
let c2 = graph.add_atom(Element::C);
let o = graph.add_atom(Element::O);
let h_o = graph.add_atom(Element::H);
let h_atoms: Vec<_> = (0..5).map(|_| graph.add_atom(Element::H)).collect();

graph.add_bond(c1, c2, GraphBondOrder::Single).unwrap();
graph.add_bond(c2, o, GraphBondOrder::Single).unwrap();
graph.add_bond(o, h_o, GraphBondOrder::Single).unwrap();
for (carbon, chunk) in [(c1, &h_atoms[0..3]), (c2, &h_atoms[3..5])].into_iter() {
    for &hydrogen in chunk {
        graph.add_bond(carbon, hydrogen, GraphBondOrder::Single).unwrap();
    }
}

let topology: MolecularTopology = assign_topology(&graph).expect("perception + typing succeed");

assert_eq!(topology.atoms[c1].atom_type, "C_3");
assert_eq!(topology.atoms[c2].atom_type, "C_3");
assert_eq!(topology.atoms[o].atom_type, "O_3");
assert_eq!(topology.atoms[h_o].atom_type, "H_HB");
```

Need custom chemistry? Parse a TOML file and extend the default rules:

```rust
use dreid_typer::{assign_topology_with_rules, rules::{get_default_rules, parse_rules}, MolecularGraph};

// Start with the default DREIDING rules
let mut all_rules = get_default_rules().to_vec();

// Parse and append custom rules from a TOML file
let extra_toml = std::fs::read_to_string("my_metals.rules.toml")?;
all_rules.extend(parse_rules(&extra_toml)?);

// Run the pipeline with extended rules
let topology = assign_topology_with_rules(&graph, &all_rules)?;
```

## Documentation

- [API Documentation](https://docs.rs/dreid-typer) - Comprehensive reference for all public types and functions.
- [Architecture Documents](docs/ARCHITECTURE.md) - In-depth design and implementation details.

## Tech Stack

- **Core Language**: Rust
- **Build System**: Cargo
- **Rule Format**: TOML

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
