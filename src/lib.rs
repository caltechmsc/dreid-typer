//! A high-performance Rust library for DREIDING force field atom typing
//! and molecular topology perception.
//!
//! # DreidTyper
//!
//! **DreidTyper** is a foundational software library for the automated assignment
//! of DREIDING atom types from molecular connectivity. It provides a modern, robust
//! solution for translating a simple chemical graph (`MolecularGraph`) into a
//! complete, simulation-ready `MolecularTopology`.
//!
//! The library is engineered in Rust for performance, memory safety, and reliability.
//! It employs a deterministic, three-phase pipeline (**Perceive-Type-Build**) to
//! ensure accurate and reproducible results.
//!
//! # Quickstart
//!
//! The primary entry point of the library is the [`assign_topology`] function.
//! Here's how to build a simple ethanol molecule and perceive its topology:
//!
//! ```
//! use dreid_typer::{
//!     assign_topology, MolecularGraph, MolecularTopology,
//!     Element, GraphBondOrder,
//! };
//!
//! // 1. Define the molecule's connectivity using a `MolecularGraph`.
//! let mut graph = MolecularGraph::new();
//! let c1 = graph.add_atom(Element::C); // Atom for the CH3 group
//! let c2 = graph.add_atom(Element::C); // Atom for the CH2 group
//! let o = graph.add_atom(Element::O);
//! let h_c1_1 = graph.add_atom(Element::H);
//! let h_c1_2 = graph.add_atom(Element::H);
//! let h_c1_3 = graph.add_atom(Element::H);
//! let h_c2_1 = graph.add_atom(Element::H);
//! let h_c2_2 = graph.add_atom(Element::H);
//! let h_o = graph.add_atom(Element::H);
//!
//! graph.add_bond(c1, c2, GraphBondOrder::Single).unwrap();
//! graph.add_bond(c2, o, GraphBondOrder::Single).unwrap();
//! graph.add_bond(c1, h_c1_1, GraphBondOrder::Single).unwrap();
//! graph.add_bond(c1, h_c1_2, GraphBondOrder::Single).unwrap();
//! graph.add_bond(c1, h_c1_3, GraphBondOrder::Single).unwrap();
//! graph.add_bond(c2, h_c2_1, GraphBondOrder::Single).unwrap();
//! graph.add_bond(c2, h_c2_2, GraphBondOrder::Single).unwrap();
//! graph.add_bond(o, h_o, GraphBondOrder::Single).unwrap();
//!
//! // 2. Call the main function to perceive the topology using default rules.
//! let topology: MolecularTopology = assign_topology(&graph).unwrap();
//!
//! // 3. Inspect the results.
//! assert_eq!(topology.atoms.len(), 9);
//! assert_eq!(topology.bonds.len(), 8);
//! assert_eq!(topology.angles.len(), 13);
//! assert_eq!(topology.torsions.len(), 12);
//!
//! // Check the assigned DREIDING atom types.
//! assert_eq!(topology.atoms[c1].atom_type, "C_3");    // sp3 Carbon
//! assert_eq!(topology.atoms[c2].atom_type, "C_3");    // sp3 Carbon
//! assert_eq!(topology.atoms[o].atom_type, "O_3");     // sp3 Oxygen
//! assert_eq!(topology.atoms[h_o].atom_type, "H_HB");  // Hydrogen-bonding Hydrogen
//! assert_eq!(topology.atoms[h_c1_1].atom_type, "H_"); // Standard Hydrogen
//! ```

mod builder;
mod core;
mod perception;
mod typing;

pub use crate::core::error::{AssignmentError, GraphValidationError, PerceptionError, TyperError};
pub use crate::core::graph::{AtomNode, BondEdge, MolecularGraph};
pub use crate::core::properties::{
    Element, GraphBondOrder, Hybridization, ParseBondOrderError, ParseElementError,
    ParseHybridizationError, TopologyBondOrder,
};
pub use crate::core::topology::{
    Angle, Atom, Bond, ImproperDihedral, MolecularTopology, ProperDihedral,
};

/// Rule parsing and customization utilities.
///
/// The core types needed to parse and inspect DREIDING
/// atom-typing rules from TOML configuration files.
pub mod rules {
    pub use crate::typing::rules::{Conditions, Rule, get_default_rules, parse_rules};
}

/// Assigns a full molecular topology using the default embedded DREIDING ruleset.
///
/// This is the primary, high-level entry point for the library. It orchestrates the
/// entire three-phase pipeline: chemical perception, atom typing, and topology
/// construction. It takes a [`MolecularGraph`] representing the chemical
/// connectivity and returns a complete [`MolecularTopology`] ready for use in
/// molecular simulations.
///
/// # Arguments
///
/// * `graph` - A reference to the [`MolecularGraph`] to be processed.
///
/// # Returns
///
/// A `Result` containing the fully perceived [`MolecularTopology`] on success.
///
/// # Errors
///
/// Returns a [`TyperError`] if any stage of the process fails. This can include:
/// * [`GraphValidationError`] if the input graph is inconsistent.
/// * [`PerceptionError`] if the chemical perception logic fails.
/// * [`AssignmentError`] if the rule engine cannot assign a type to one or more atoms.
///
/// # Panics
///
/// Panics if the embedded default rules file is malformed, which indicates a
/// critical library bug.
pub fn assign_topology(graph: &MolecularGraph) -> Result<MolecularTopology, TyperError> {
    let default_rules = typing::rules::get_default_rules();
    assign_topology_internal(graph, default_rules)
}

/// Assigns a full molecular topology using a custom set of typing rules.
///
/// This function provides the same functionality as [`assign_topology`] but allows
/// for customization of the atom typing logic by supplying a custom slice of [`rules::Rule`]s.
/// This is useful for extending the DREIDING force field to new elements or
/// defining special types for specific chemical environments.
///
/// # Arguments
///
/// * `graph` - A reference to the [`MolecularGraph`] to be processed.
/// * `rules` - A slice of [`rules::Rule`] structs that the typing engine will use.
///
/// # Errors
///
/// Returns a [`TyperError`] under the same conditions as [`assign_topology`].
pub fn assign_topology_with_rules(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    assign_topology_internal(graph, rules)
}

/// Internal core function that executes the perception, typing, and building pipeline.
fn assign_topology_internal(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    let annotated_molecule = perception::perceive(graph)?;

    let atom_types = typing::engine::assign_types(&annotated_molecule, rules)
        .map_err(TyperError::AssignmentFailed)?;

    let topology = builder::build_topology(&annotated_molecule, &atom_types);

    Ok(topology)
}
