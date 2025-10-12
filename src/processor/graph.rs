use crate::core::error::GraphValidationError;
use crate::core::graph::MolecularGraph;
use crate::core::{BondOrder, Element, Hybridization};
use std::collections::HashSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerceptionSource {
    Generic,
    Template,
}

#[derive(Debug, Clone)]
pub struct AtomView {
    pub id: usize,
    pub element: Element,
    pub formal_charge: i8,
    pub degree: u8,
    pub valence_electrons: u8,
    pub bonding_electrons: u8,
    pub lone_pairs: u8,
    pub steric_number: u8,
    pub hybridization: Hybridization,
    pub is_in_ring: bool,
    pub smallest_ring_size: Option<u8>,
    pub is_aromatic: bool,
    pub perception_source: Option<PerceptionSource>,
}

#[derive(Debug, Clone)]
pub struct ProcessingGraph {
    pub atoms: Vec<AtomView>,
    pub adjacency: Vec<Vec<(usize, BondOrder)>>,
}

#[derive(Debug, Clone, Default)]
pub struct RingInfo(pub HashSet<Vec<usize>>);

impl ProcessingGraph {
    pub fn new(graph: &MolecularGraph) -> Result<Self, GraphValidationError> {
        let num_atoms = graph.atoms.len();
        let mut adjacency = vec![vec![]; num_atoms];

        for bond in &graph.bonds {
            let (u, v) = bond.atom_ids;
            if u >= num_atoms || v >= num_atoms {
                let invalid_id = if u >= num_atoms { u } else { v };
                return Err(GraphValidationError::MissingAtom { id: invalid_id });
            }

            adjacency[u].push((v, bond.order));
            adjacency[v].push((u, bond.order));
        }

        let atoms = graph
            .atoms
            .iter()
            .map(|atom_node| {
                let degree = adjacency[atom_node.id].len() as u8;

                AtomView {
                    id: atom_node.id,
                    element: atom_node.element,
                    formal_charge: atom_node.formal_charge,
                    degree,
                    valence_electrons: 0,
                    bonding_electrons: 0,
                    lone_pairs: 0,
                    steric_number: 0,
                    hybridization: Hybridization::Unknown,
                    is_in_ring: false,
                    smallest_ring_size: None,
                    is_aromatic: false,
                    perception_source: None,
                }
            })
            .collect();

        Ok(Self { atoms, adjacency })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element};

    #[test]
    fn new_processing_graph_for_methane_is_correct() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C);
        let h1 = mg.add_atom(Element::H);
        let h2 = mg.add_atom(Element::H);
        let h3 = mg.add_atom(Element::H);
        let h4 = mg.add_atom(Element::H);
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();
        mg.add_bond(c1, h2, BondOrder::Single).unwrap();
        mg.add_bond(c1, h3, BondOrder::Single).unwrap();
        mg.add_bond(c1, h4, BondOrder::Single).unwrap();

        let pg = ProcessingGraph::new(&mg).unwrap();

        assert_eq!(pg.atoms.len(), 5);
        assert_eq!(pg.adjacency.len(), 5);

        let carbon = &pg.atoms[c1];
        assert_eq!(carbon.id, c1);
        assert_eq!(carbon.element, Element::C);
        assert_eq!(carbon.degree, 4);
        assert_eq!(carbon.hybridization, Hybridization::Unknown);
        assert!(!carbon.is_in_ring);
        assert!(!carbon.is_aromatic);
        assert_eq!(carbon.smallest_ring_size, None);

        let hydrogen = &pg.atoms[h1];
        assert_eq!(hydrogen.degree, 1);

        assert_eq!(pg.adjacency[c1].len(), 4);
        assert_eq!(pg.adjacency[h1].len(), 1);
        assert_eq!(pg.adjacency[h1][0], (c1, BondOrder::Single));
    }

    #[test]
    fn new_processing_graph_for_ethene_is_correct() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C);
        let c2 = mg.add_atom(Element::C);
        mg.add_bond(c1, c2, BondOrder::Double).unwrap();

        let pg = ProcessingGraph::new(&mg).unwrap();
        assert_eq!(pg.atoms.len(), 2);

        assert_eq!(pg.atoms[c1].degree, 1);
        assert_eq!(pg.atoms[c2].degree, 1);

        assert_eq!(pg.adjacency[c1][0], (c2, BondOrder::Double));
        assert_eq!(pg.adjacency[c2][0], (c1, BondOrder::Double));
    }

    #[test]
    fn new_processing_graph_with_invalid_bond_returns_error() {
        let mut mg = MolecularGraph::new();
        mg.add_atom(Element::C);
        mg.add_bond(0, 1, BondOrder::Single).unwrap_err();

        let bad_mg = MolecularGraph {
            atoms: vec![crate::core::graph::AtomNode {
                id: 0,
                element: Element::C,
            }],
            bonds: vec![crate::core::graph::BondEdge {
                id: 0,
                atom_ids: (0, 1),
                order: BondOrder::Single,
            }],
        };

        let result = ProcessingGraph::new(&bad_mg);
        assert!(result.is_err());
        match result.unwrap_err() {
            GraphValidationError::MissingAtom { id } => assert_eq!(id, 1),
            _ => panic!("Expected MissingAtom error"),
        }
    }
}
