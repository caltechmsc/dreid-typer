//! A high-performance Rust library for DREIDING force field atom typing
//! and molecular topology perception.
//!
//! # DreidTyper
//!
//! **DreidTyper** is a foundational software library for the automated assignment
//! of DREIDING force field atom types and the perception of molecular topologies.
//! It provides a modern, robust solution for translating simple chemical connectivity
//! (a `MolecularGraph`) into a complete, engine-agnostic topological description
//! (`MolecularTopology`) essential for molecular simulations.
//!
//! The core mission of DreidTyper is to provide a reliable, predictable, and
//! easy-to-integrate tool for developers and researchers building the next
//! generation of simulation tools for general chemistry, materials science,
//! and drug discovery.
//!
//! # Quickstart
//!
//! The primary entry point of the library is the [`assign_topology`] function.
//! Here's how to build a simple ethanol molecule and perceive its topology:
//!
//! ```
//! use dreid_typer::{
//!     assign_topology, MolecularGraph, MolecularTopology,
//!     Element, BondOrder,
//! };
//!
//! // 1. Define the molecule's connectivity using a `MolecularGraph`.
//! let mut graph = MolecularGraph::new();
//! let c1 = graph.add_atom(Element::C); // CH3
//! let c2 = graph.add_atom(Element::C); // CH2
//! let o = graph.add_atom(Element::O);
//! let h_c1_1 = graph.add_atom(Element::H);
//! let h_c1_2 = graph.add_atom(Element::H);
//! let h_c1_3 = graph.add_atom(Element::H);
//! let h_c2_1 = graph.add_atom(Element::H);
//! let h_c2_2 = graph.add_atom(Element::H);
//! let h_o = graph.add_atom(Element::H);
//!
//! graph.add_bond(c1, c2, BondOrder::Single).unwrap();
//! graph.add_bond(c2, o, BondOrder::Single).unwrap();
//! graph.add_bond(c1, h_c1_1, BondOrder::Single).unwrap();
//! graph.add_bond(c1, h_c1_2, BondOrder::Single).unwrap();
//! graph.add_bond(c1, h_c1_3, BondOrder::Single).unwrap();
//! graph.add_bond(c2, h_c2_1, BondOrder::Single).unwrap();
//! graph.add_bond(c2, h_c2_2, BondOrder::Single).unwrap();
//! graph.add_bond(o, h_o, BondOrder::Single).unwrap();
//!
//! // 2. Call the main function to perceive the topology.
//! let topology: MolecularTopology = assign_topology(&graph).unwrap();
//!
//! // 3. Inspect the results.
//! assert_eq!(topology.atoms.len(), 9);
//! assert_eq!(topology.bonds.len(), 8);
//! assert_eq!(topology.angles.len(), 13);
//! assert_eq!(topology.proper_dihedrals.len(), 12);
//!
//! // Check the assigned DREIDING atom types.
//! assert_eq!(topology.atoms[c1].atom_type, "C_3");   // sp3 Carbon
//! assert_eq!(topology.atoms[c2].atom_type, "C_3");   // sp3 Carbon
//! assert_eq!(topology.atoms[o].atom_type, "O_3");    // sp3 Oxygen
//! assert_eq!(topology.atoms[h_o].atom_type, "H_HB"); // Hydrogen-bonding Hydrogen
//! assert_eq!(topology.atoms[h_c1_1].atom_type, "H_"); // Standard Hydrogen
//! ```

// Internal modules responsible for the pipeline stages.
mod builder;
mod core;
mod processor;

// Re-export core data structures for convenient access at the crate root.
pub use crate::core::graph::{
    Angle, Atom, Bond, ImproperDihedral, MolecularGraph, MolecularTopology, ProperDihedral,
};
pub use crate::core::{BondOrder, Element, Hybridization};

// Re-export error types for comprehensive error handling.
pub use crate::core::error::{AnnotationError, AssignmentError, GraphValidationError, TyperError};

/// Provides functionality for parsing and managing DREIDING atom typing rules.
///
/// While the library can be used out-of-the-box with its default ruleset,
/// this module allows advanced users to define and load custom rules from
/// TOML-formatted strings using the [`rules::parse_rules`] function.
pub mod rules;

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
/// * [`GraphValidationError`] if the input graph is inconsistent (e.g., dangling bonds).
/// * [`AnnotationError`] if the chemical perception logic fails (e.g., cannot determine hybridization).
/// * [`AssignmentError`] if the rule engine cannot assign a type to one or more atoms.
///
/// # Examples
///
/// See the [crate-level documentation](crate) for a detailed example.
pub fn assign_topology(graph: &MolecularGraph) -> Result<MolecularTopology, TyperError> {
    // Retrieve the statically compiled default rules.
    let default_rules = rules::get_default_rules()?;
    assign_topology_internal(graph, default_rules)
}

/// Assigns a full molecular topology using a user-provided set of rules.
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
///
/// # Examples
///
/// ```
/// use dreid_typer::{
///     assign_topology_with_rules, rules, MolecularGraph, Element, BondOrder
/// };
///
/// // Define a simple molecule: a single Carbon atom.
/// let mut graph = MolecularGraph::new();
/// let c = graph.add_atom(Element::C);
///
/// // Define a custom, high-priority rule for a lone carbon atom.
/// let custom_rules_toml = r#"
///     [[rule]]
///     name = "Lone_Carbon"
///     priority = 1000
///     type = "C_LONE"
///     conditions = { element = "C", degree = 0 }
/// "#;
///
/// let my_rules = rules::parse_rules(custom_rules_toml).unwrap();
///
/// // Assign topology using the custom rules.
/// let topology = assign_topology_with_rules(&graph, &my_rules).unwrap();
///
/// assert_eq!(topology.atoms[c].atom_type, "C_LONE");
/// ```
pub fn assign_topology_with_rules(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    assign_topology_internal(graph, rules)
}

/// Internal core function that executes the perception and typing pipeline.
///
/// This function is the shared implementation for both public-facing `assign_topology`
/// functions. It encapsulates the three-phase process:
/// 1. Perceive chemical properties to create a `ProcessingGraph`.
/// 2. Assign atom types using the provided ruleset.
/// 3. Build the final `MolecularTopology`.
fn assign_topology_internal(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    // Phase 1: Perceive chemical features from the input graph.
    let perception = processor::perceive(graph)?;
    let processing_graph = perception.processing_graph;

    // Phase 2: Run the iterative typing engine with the given rules.
    let atom_types =
        processor::assign_types(&processing_graph, rules).map_err(TyperError::AssignmentFailed)?;

    // Phase 3: Construct the final topology from the annotated graph and types.
    let topology = builder::build_topology(graph, &processing_graph, &atom_types)?;

    Ok(topology)
}
