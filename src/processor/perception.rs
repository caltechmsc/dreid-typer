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

pub(crate) fn perceive_generic_aromaticity(
    graph: &mut ProcessingGraph,
    ring_info: &RingInfo,
) -> Result<(), AnnotationError> {
    let mut aromatic_atoms = HashSet::new();

    for ring_atom_ids in &ring_info.0 {
        if is_ring_aromatic(ring_atom_ids, graph) {
            aromatic_atoms.extend(ring_atom_ids.iter().copied());
        }
    }

    for atom_id in aromatic_atoms {
        if graph.atoms[atom_id].perception_source == Some(PerceptionSource::Template) {
            continue;
        }
        graph.atoms[atom_id].is_aromatic = true;
    }

    Ok(())
}

pub(crate) fn perceive_generic_hybridization(
    graph: &mut ProcessingGraph,
) -> Result<(), AnnotationError> {
    let mut initial_hybs = Vec::with_capacity(graph.atoms.len());
    for i in 0..graph.atoms.len() {
        if graph.atoms[i].perception_source == Some(PerceptionSource::Template) {
            initial_hybs.push(graph.atoms[i].hybridization);
            continue;
        }

        let atom = &mut graph.atoms[i];
        let steric = atom.degree + atom.lone_pairs;
        atom.steric_number = steric;

        let hyb = if atom.is_aromatic {
            Hybridization::Resonant
        } else if is_special_non_hybridized(atom.element) || atom.degree == 0 {
            Hybridization::None
        } else {
            match steric {
                4 => Hybridization::SP3,
                3 => Hybridization::SP2,
                2 => Hybridization::SP,
                _ => Hybridization::Unknown,
            }
        };
        initial_hybs.push(hyb);
    }

    for i in 0..graph.atoms.len() {
        if graph.atoms[i].perception_source == Some(PerceptionSource::Template) {
            continue;
        }

        let mut final_hyb = initial_hybs[i];

        if final_hyb == Hybridization::SP3 && graph.atoms[i].lone_pairs > 0 {
            let atom_view = &graph.atoms[i];
            let hydrogen_neighbors = graph.adjacency[atom_view.id]
                .iter()
                .filter(|(neighbor_id, _)| graph.atoms[*neighbor_id].element == Element::H)
                .count();
            let has_protective_hydrogen = atom_view.element == Element::O && hydrogen_neighbors > 0;
            let is_multihydrogen = hydrogen_neighbors >= 2;
            let is_anion = atom_view.formal_charge < 0;

            if !has_protective_hydrogen && !is_anion && !is_multihydrogen {
                let has_pi_system_neighbor =
                    graph.adjacency[atom_view.id]
                        .iter()
                        .any(|(neighbor_id, _)| {
                            let neighbor_hyb = initial_hybs[*neighbor_id];
                            matches!(
                                neighbor_hyb,
                                Hybridization::SP2 | Hybridization::SP | Hybridization::Resonant
                            )
                        });

                if has_pi_system_neighbor {
                    final_hyb = Hybridization::SP2;
                }
            }
        }

        if final_hyb == Hybridization::Unknown {
            return Err(AnnotationError::HybridizationInference {
                atom_id: graph.atoms[i].id,
            });
        }

        graph.atoms[i].hybridization = final_hyb;
        graph.atoms[i].perception_source = Some(PerceptionSource::Generic);

        graph.atoms[i].steric_number = match final_hyb {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            Hybridization::None => graph.atoms[i].degree + graph.atoms[i].lone_pairs,
            Hybridization::Unknown => graph.atoms[i].steric_number,
        };
    }

    Ok(())
}

fn is_ring_aromatic(ring_atom_ids: &[usize], graph: &ProcessingGraph) -> bool {
    if ring_atom_ids.len() < 3 {
        return false;
    }

    let ring_atoms: HashSet<usize> = ring_atom_ids.iter().copied().collect();

    let all_bonds_aromatic = ring_atom_ids.iter().all(|&atom_id| {
        graph.adjacency[atom_id].iter().all(|(neighbor, order)| {
            if ring_atoms.contains(neighbor) {
                matches!(order, BondOrder::Aromatic)
            } else {
                true
            }
        })
    });

    if all_bonds_aromatic {
        return true;
    }

    let mut total_pi_electrons = 0u8;

    for &atom_id in ring_atom_ids {
        match count_atom_pi_contribution(atom_id, &ring_atoms, graph) {
            Some(contribution) => total_pi_electrons += contribution,
            None => return false,
        }
    }

    if total_pi_electrons < 2 {
        return false;
    }

    (total_pi_electrons - 2) % 4 == 0
}

fn count_atom_pi_contribution(
    atom_id: usize,
    ring_atoms: &HashSet<usize>,
    graph: &ProcessingGraph,
) -> Option<u8> {
    let atom = &graph.atoms[atom_id];
    let steric = match atom.perception_source {
        Some(PerceptionSource::Template) => atom.steric_number,
        _ => atom.degree + atom.lone_pairs,
    };

    if steric >= 4 && atom.lone_pairs == 0 {
        return None;
    }

    let has_exocyclic_pi_bond = graph.adjacency[atom_id].iter().any(|(neighbor, order)| {
        !ring_atoms.contains(neighbor) && matches!(order, BondOrder::Double | BondOrder::Triple)
    });

    if has_exocyclic_pi_bond {
        return Some(1);
    }

    match steric {
        2 | 3 => Some(1),
        4 if atom.lone_pairs > 0 => Some(2),
        _ => None,
    }
}

fn bond_order_contribution(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 1,
    }
}

fn is_special_non_hybridized(element: Element) -> bool {
    matches!(
        element,
        Element::H
            | Element::F
            | Element::Cl
            | Element::Br
            | Element::I
            | Element::He
            | Element::Ne
            | Element::Ar
            | Element::Kr
            | Element::Xe
            | Element::Rn
            | Element::Li
            | Element::Na
            | Element::K
            | Element::Rb
            | Element::Cs
            | Element::Fr
            | Element::Be
            | Element::Mg
            | Element::Ca
            | Element::Sr
            | Element::Ba
            | Element::Ra
            | Element::Fe
            | Element::Zn
    )
}
