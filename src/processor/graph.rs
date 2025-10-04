use crate::core::graph::MolecularGraph;
use crate::core::{BondOrder, Element, Hybridization};
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct AtomView {
    pub id: usize,
    pub element: Element,
    pub degree: u8,
    pub hybridization: Hybridization,
    pub is_in_ring: bool,
    pub smallest_ring_size: Option<u8>,
    pub is_aromatic: bool,
}

#[derive(Debug, Clone)]
pub struct ProcessingGraph {
    pub atoms: Vec<AtomView>,
    pub adjacency: Vec<Vec<(usize, BondOrder)>>,
}

pub struct RingInfo(pub HashSet<HashSet<usize>>);

impl ProcessingGraph {
    pub fn new(graph: &MolecularGraph) -> Result<Self, &'static str> {
        let num_atoms = graph.atoms.len();
        let mut adjacency = vec![vec![]; num_atoms];
        for bond in &graph.bonds {
            let (u, v) = bond.atom_ids;
            if u >= num_atoms || v >= num_atoms {
                return Err("Bond references an out-of-bounds atom ID.");
            }
            adjacency[u].push((v, bond.order));
            adjacency[v].push((u, bond.order));
        }

        let atoms = graph
            .atoms
            .iter()
            .map(|atom_node| AtomView {
                id: atom_node.id,
                element: atom_node.element,
                degree: adjacency[atom_node.id].len() as u8,
                hybridization: Hybridization::Unknown,
                is_in_ring: false,
                smallest_ring_size: None,
                is_aromatic: false,
            })
            .collect();

        Ok(Self { atoms, adjacency })
    }
}
