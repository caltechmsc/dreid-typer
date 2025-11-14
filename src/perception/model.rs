use crate::core::error::GraphValidationError;
use crate::core::graph::{BondEdge, MolecularGraph};
use crate::core::properties::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone)]
pub struct AnnotatedAtom {
    pub id: usize,
    pub element: Element,

    pub formal_charge: i8,
    pub lone_pairs: u8,
    pub degree: u8,

    pub is_in_ring: bool,
    pub smallest_ring_size: Option<u8>,

    pub is_aromatic: bool,
    pub is_anti_aromatic: bool,

    pub is_in_conjugated_system: bool,
    pub is_resonant: bool,
    pub steric_number: u8,
    pub hybridization: Hybridization,
}

#[derive(Debug, Clone)]
pub struct AnnotatedMolecule {
    pub atoms: Vec<AnnotatedAtom>,
    pub bonds: Vec<BondEdge>,
    pub adjacency: Vec<Vec<(usize, BondOrder)>>,
}

impl AnnotatedMolecule {
    pub fn new(graph: &MolecularGraph) -> Result<Self, GraphValidationError> {
        let mut adjacency = vec![vec![]; graph.atoms.len()];
        for bond in &graph.bonds {
            let (u, v) = bond.atom_ids;

            if u >= graph.atoms.len() {
                return Err(GraphValidationError::MissingAtom { atom_id: u });
            }
            if v >= graph.atoms.len() {
                return Err(GraphValidationError::MissingAtom { atom_id: v });
            }

            adjacency[u].push((v, bond.order));
            adjacency[v].push((u, bond.order));
        }

        let atoms = graph
            .atoms
            .iter()
            .map(|node| AnnotatedAtom {
                id: node.id,
                element: node.element,
                degree: adjacency[node.id].len() as u8,
                formal_charge: 0,
                lone_pairs: 0,
                is_in_ring: false,
                smallest_ring_size: None,
                is_aromatic: false,
                is_anti_aromatic: false,
                is_in_conjugated_system: false,
                is_resonant: false,
                steric_number: 0,
                hybridization: Hybridization::Unknown,
            })
            .collect();

        Ok(Self {
            atoms,
            bonds: graph.bonds.clone(),
            adjacency,
        })
    }
}
