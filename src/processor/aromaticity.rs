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

fn is_ring_aromatic(ring_atom_ids: &[usize], graph: &ProcessingGraph) -> bool {
    let ring_size = ring_atom_ids.len();
    if !(5..=7).contains(&ring_size) {
        return false;
    }

    let all_can_participate = ring_atom_ids.iter().all(|&id| {
        let prov_hyb = calculate_provisional_hybridization(&graph.atoms[id], graph);
        matches!(
            prov_hyb,
            ProvisionalHybridization::SP | ProvisionalHybridization::SP2
        )
    });

    if !all_can_participate {
        return false;
    }

    let pi_electron_count = count_pi_electrons(ring_atom_ids, graph);

    pi_electron_count >= 2 && (pi_electron_count - 2) % 4 == 0
}

fn count_pi_electrons(ring_atom_ids: &[usize], graph: &ProcessingGraph) -> u8 {
    ring_atom_ids
        .iter()
        .map(|&id| count_atom_pi_contribution(&graph.atoms[id], graph))
        .sum()
}

fn count_atom_pi_contribution(atom: &AtomView, graph: &ProcessingGraph) -> u8 {
    let has_pi_bond = graph.adjacency[atom.id]
        .iter()
        .any(|(_, order)| order.is_pi_bond());

    if has_pi_bond {
        return 1;
    }

    match atom.element {
        Element::N | Element::P | Element::As if atom.degree == 3 => 2,
        Element::O | Element::S | Element::Se if atom.degree == 2 => 2,
        Element::B if atom.degree == 3 => 0,
        _ => 0,
    }
}

use crate::core::BondOrder;
impl BondOrder {
    fn is_pi_bond(&self) -> bool {
        matches!(
            self,
            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
        )
    }
}
