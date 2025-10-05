use crate::core::Hybridization;
use crate::core::error::TyperError;
use crate::core::graph::{
    Angle, Atom, Bond, ImproperDihedral, MolecularGraph, MolecularTopology, ProperDihedral,
};
use crate::processor::ProcessingGraph;
use std::collections::HashSet;

pub(crate) fn build_topology(
    initial_graph: &MolecularGraph,
    processing_graph: &ProcessingGraph,
    atom_types: &[String],
) -> Result<MolecularTopology, TyperError> {
    let atoms = build_atoms(processing_graph, atom_types);
    let bonds = build_bonds(initial_graph);
    let angles = build_angles(processing_graph);
    let proper_dihedrals = build_proper_dihedrals(processing_graph);
    let improper_dihedrals = build_improper_dihedrals(processing_graph);

    Ok(MolecularTopology {
        atoms,
        bonds: bonds.into_iter().collect(),
        angles: angles.into_iter().collect(),
        proper_dihedrals: proper_dihedrals.into_iter().collect(),
        improper_dihedrals: improper_dihedrals.into_iter().collect(),
    })
}
