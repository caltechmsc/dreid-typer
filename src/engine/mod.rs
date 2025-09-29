pub mod analyze;
pub mod assign;
pub mod perceive;

use crate::core::{BondOrder, Element, Hybridization};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct AnnotatedAtom {
    pub id: usize,
    pub element: Element,

    pub degree: u8,
    pub valence: u8,
    pub hybridization: Hybridization,
    pub is_in_ring: bool,
    pub is_aromatic: bool,
    pub smallest_ring_size: u8,

    pub atom_type: Option<String>,
}

#[derive(Debug, Clone)]
pub struct AnnotatedGraph {
    pub atoms: Vec<AnnotatedAtom>,
    pub adjacency: Vec<Vec<(usize, BondOrder)>>,
}

#[derive(Debug, Clone, Default)]
pub struct GeometricPerception {
    pub atom_ring_info: HashMap<usize, Vec<usize>>,
}

#[derive(Debug, Clone, Default)]
pub struct ChemicalPerception {
    pub aromatic_atoms: HashSet<usize>,
    pub aromatic_bonds: HashSet<(usize, usize)>,
}
