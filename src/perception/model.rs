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
