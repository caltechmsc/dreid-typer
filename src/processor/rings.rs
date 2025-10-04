use super::graph::{ProcessingGraph, RingInfo};
use std::collections::{HashMap, HashSet, VecDeque};

pub(crate) fn perceive_rings(graph: &ProcessingGraph) -> RingInfo {
    if graph.atoms.is_empty() {
        return RingInfo::default();
    }
    let mut finder = JohnsonCycleFinder::new(graph);
    let sorted_vec_cycles = finder.find_cycles_internal();
    RingInfo(sorted_vec_cycles)
}
