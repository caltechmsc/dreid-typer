mod builder;
mod core;
mod perception;
mod typing;

pub use crate::core::error::{AssignmentError, GraphValidationError, PerceptionError, TyperError};
pub use crate::core::graph::{
    Angle, Atom, Bond, ImproperDihedral, MolecularGraph, MolecularTopology, ProperDihedral,
};
pub use crate::core::properties::{BondOrder, Element, Hybridization};

pub mod rules {
    pub use crate::typing::rules::{Rule, parse_rules};
}

pub fn assign_topology(graph: &MolecularGraph) -> Result<MolecularTopology, TyperError> {
    let default_rules = typing::rules::get_default_rules();
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
    let annotated_molecule = perception::perceive(graph)?;

    let atom_types = typing::engine::assign_types(&annotated_molecule, rules)
        .map_err(TyperError::AssignmentFailed)?;

    let topology = builder::build_topology(&annotated_molecule, &atom_types);

    Ok(topology)
}
