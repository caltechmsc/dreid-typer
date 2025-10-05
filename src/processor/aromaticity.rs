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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element};

    fn setup(
        edges: &[(usize, usize)],
        elements: &[(usize, Element)],
    ) -> (ProcessingGraph, RingInfo) {
        let mut mg = MolecularGraph::new();
        let max_node = edges
            .iter()
            .flat_map(|(u, v)| vec![*u, *v])
            .max()
            .unwrap_or(0)
            .max(elements.iter().map(|(i, _)| *i).max().unwrap_or(0));

        for _ in 0..=max_node {
            mg.add_atom(Element::C);
        }
        for (i, el) in elements {
            mg.atoms[*i].element = *el;
        }

        for (u, v) in edges {
            mg.add_bond(*u, *v, BondOrder::Aromatic).unwrap();
        }

        let mut pg = ProcessingGraph::new(&mg).unwrap();
        let ring_info = crate::processor::rings::perceive_rings(&pg);
        for (id, atom) in pg.atoms.iter_mut().enumerate() {
            if ring_info.0.iter().any(|r| r.contains(&id)) {
                atom.is_in_ring = true;
            }
        }
        (pg, ring_info)
    }

    #[test]
    fn test_benzene_is_aromatic() {
        let edges = &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        let (mut pg, ring_info) = setup(edges, &[]);

        perceive_aromaticity(&mut pg, &ring_info).unwrap();

        for atom in pg.atoms {
            assert!(
                atom.is_aromatic,
                "Atom {} in benzene should be aromatic",
                atom.id
            );
        }
    }

    #[test]
    fn test_cyclohexane_is_not_aromatic() {
        let mut mg = MolecularGraph::new();
        for _ in 0..6 {
            mg.add_atom(Element::C);
        }
        for i in 0..6 {
            mg.add_bond(i, (i + 1) % 6, BondOrder::Single).unwrap();
        }

        let mut pg = ProcessingGraph::new(&mg).unwrap();
        let ring_info = crate::processor::rings::perceive_rings(&pg);
        for (id, atom) in pg.atoms.iter_mut().enumerate() {
            if ring_info.0.iter().any(|r| r.contains(&id)) {
                atom.is_in_ring = true;
            }
        }

        perceive_aromaticity(&mut pg, &ring_info).unwrap();

        for atom in pg.atoms {
            assert!(
                !atom.is_aromatic,
                "Atom {} in cyclohexane should not be aromatic",
                atom.id
            );
        }
    }

    #[test]
    fn test_pyridine_is_aromatic() {
        let edges = &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        let elements = &[(0, Element::N)];
        let (mut pg, ring_info) = setup(edges, elements);

        perceive_aromaticity(&mut pg, &ring_info).unwrap();

        assert!(pg.atoms.iter().all(|a| a.is_aromatic));
    }

    #[test]
    fn test_pyrrole_is_aromatic() {
        let mut mg = MolecularGraph::new();
        mg.add_atom(Element::N);
        mg.add_atom(Element::C);
        mg.add_atom(Element::C);
        mg.add_atom(Element::C);
        mg.add_atom(Element::C);
        mg.add_atom(Element::H);
        mg.add_atom(Element::H);
        mg.add_atom(Element::H);
        mg.add_atom(Element::H);
        mg.add_atom(Element::H);
        mg.add_bond(0, 1, BondOrder::Single).unwrap();
        mg.add_bond(1, 2, BondOrder::Double).unwrap();
        mg.add_bond(2, 3, BondOrder::Single).unwrap();
        mg.add_bond(3, 4, BondOrder::Double).unwrap();
        mg.add_bond(4, 0, BondOrder::Single).unwrap();
        mg.add_bond(0, 5, BondOrder::Single).unwrap();
        mg.add_bond(1, 6, BondOrder::Single).unwrap();
        mg.add_bond(2, 7, BondOrder::Single).unwrap();
        mg.add_bond(3, 8, BondOrder::Single).unwrap();
        mg.add_bond(4, 9, BondOrder::Single).unwrap();

        let mut pg = ProcessingGraph::new(&mg).unwrap();
        let ring_info = crate::processor::rings::perceive_rings(&pg);
        for (id, atom) in pg.atoms.iter_mut().enumerate() {
            if ring_info.0.iter().any(|r| r.contains(&id)) {
                atom.is_in_ring = true;
            }
        }

        perceive_aromaticity(&mut pg, &ring_info).unwrap();

        for atom in pg.atoms.iter().take(5) {
            assert!(
                atom.is_aromatic,
                "Atom {} in pyrrole ring should be aromatic",
                atom.id
            );
        }
        for atom in pg.atoms.iter().skip(5) {
            assert!(
                !atom.is_aromatic,
                "Atom {} (H) in pyrrole should not be aromatic",
                atom.id
            );
        }
    }
}
