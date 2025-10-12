#![doc = include_str!("../README.md")]

mod builder;
mod core;
mod processor;

pub use crate::core::graph::{
    Angle, Atom, Bond, ImproperDihedral, MolecularGraph, MolecularTopology, ProperDihedral,
};
pub use crate::core::{BondOrder, Element, Hybridization};

pub use crate::core::error::{AnnotationError, AssignmentError, GraphValidationError, TyperError};

pub mod rules;

pub fn assign_topology(graph: &MolecularGraph) -> Result<MolecularTopology, TyperError> {
    let default_rules = rules::get_default_rules()?;
    assign_topology_internal(graph, default_rules)
}

pub fn assign_topology_with_rules(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    assign_topology_internal(graph, rules)
}

fn assign_topology_internal(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    let perception = processor::perceive(graph)?;
    let processing_graph = perception.processing_graph;

    let atom_types =
        processor::assign_types(&processing_graph, rules).map_err(TyperError::AssignmentFailed)?;

    let topology = builder::build_topology(graph, &processing_graph, &atom_types)?;

    Ok(topology)
}
