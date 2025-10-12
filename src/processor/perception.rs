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
