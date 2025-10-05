use super::graph::{AtomView, ProcessingGraph, RingInfo};
use super::hybridization::{ProvisionalHybridization, calculate_provisional_hybridization};
use crate::core::Element;
use crate::core::error::TyperError;
use std::collections::HashSet;

pub(crate) fn perceive_aromaticity(
    graph: &mut ProcessingGraph,
    ring_info: &RingInfo,
) -> Result<(), TyperError> {
    let mut aromatic_atoms = HashSet::new();

    for ring_atom_ids_vec in &ring_info.0 {
        if is_ring_aromatic(ring_atom_ids_vec, graph) {
            for &atom_id in ring_atom_ids_vec {
                aromatic_atoms.insert(atom_id);
            }
        }
    }

    for atom_id in aromatic_atoms {
        graph.atoms[atom_id].is_aromatic = true;
    }

    Ok(())
}
