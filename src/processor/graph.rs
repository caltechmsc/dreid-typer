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
