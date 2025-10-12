use super::graph::{PerceptionSource, ProcessingGraph, RingInfo};
use crate::core::error::{AnnotationError, TyperError};
use crate::core::graph::MolecularGraph;
use crate::core::{BondOrder, Element, Hybridization};
use std::collections::{HashSet, VecDeque};

pub(crate) fn perceive_electron_counts(
    molecular_graph: &MolecularGraph,
) -> Result<ProcessingGraph, TyperError> {
    let mut graph = ProcessingGraph::new(molecular_graph).map_err(TyperError::InvalidInputGraph)?;

    for atom in &mut graph.atoms {
        let valence = get_valence_electrons(atom.element).unwrap_or(0);

        atom.valence_electrons = valence;

        let bonding = graph.adjacency[atom.id]
            .iter()
            .map(|(_, order)| bond_order_contribution(*order))
            .sum::<u8>();
        atom.bonding_electrons = bonding;

        let available = valence as i16 - bonding as i16 - atom.formal_charge as i16;
        let adjusted = available.max(0);
        let lone_pairs = (adjusted / 2) as u8;
        atom.lone_pairs = lone_pairs;
        atom.steric_number = 0;
        atom.hybridization = Hybridization::Unknown;
        atom.is_aromatic = false;
        atom.is_in_ring = false;
        atom.smallest_ring_size = None;
        atom.perception_source = None;
    }

    Ok(graph)
}

pub(crate) fn perceive_rings(graph: &ProcessingGraph) -> RingInfo {
    if graph.atoms.is_empty() {
        return RingInfo::default();
    }

    let mut finder = JohnsonCycleFinder::new(graph);
    let sorted_vec_cycles = finder.find_cycles_internal();
    RingInfo(sorted_vec_cycles)
}

pub(crate) fn apply_ring_annotations(graph: &mut ProcessingGraph, ring_info: &RingInfo) {
    let mut atom_ring_sizes: Vec<Vec<u8>> = vec![vec![]; graph.atoms.len()];
    for ring in &ring_info.0 {
        let ring_len = ring.len() as u8;
        for &atom_id in ring {
            atom_ring_sizes[atom_id].push(ring_len);
        }
    }

    for (atom_id, atom) in graph.atoms.iter_mut().enumerate() {
        if !atom_ring_sizes[atom_id].is_empty() {
            atom.is_in_ring = true;
            atom.smallest_ring_size = atom_ring_sizes[atom_id].iter().min().copied();
        }
    }
}

pub(crate) fn perceive_generic_properties(
    graph: &mut ProcessingGraph,
    ring_info: &RingInfo,
) -> Result<(), AnnotationError> {
    perceive_generic_aromaticity(graph, ring_info)?;
    perceive_generic_hybridization(graph)?;
    Ok(())
}
