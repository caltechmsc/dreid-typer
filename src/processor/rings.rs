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

struct JohnsonCycleFinder<'a> {
    graph: &'a ProcessingGraph,
    blocked: Vec<bool>,
    b_sets: Vec<HashSet<usize>>,
    stack: Vec<usize>,
    all_cycles: HashSet<Vec<usize>>,
}

impl<'a> JohnsonCycleFinder<'a> {
    fn new(graph: &'a ProcessingGraph) -> Self {
        let num_atoms = graph.atoms.len();
        Self {
            graph,
            blocked: vec![false; num_atoms],
            b_sets: vec![HashSet::new(); num_atoms],
            stack: Vec::new(),
            all_cycles: HashSet::new(),
        }
    }

    fn find_cycles_internal(&mut self) -> HashSet<Vec<usize>> {
        let num_atoms = self.graph.atoms.len();
        for i in 0..num_atoms {
            self.find_cycles_from_node(i);
        }
        self.all_cycles.clone()
    }

    fn find_cycles_from_node(&mut self, start_node: usize) {
        let num_atoms = self.graph.atoms.len();
        let mut queue = VecDeque::new();
        queue.push_back(vec![start_node]);

        let mut paths = HashMap::new();
        paths.insert(start_node, vec![vec![start_node]]);

        while let Some(path) = queue.pop_front() {
            let last_node = *path.last().unwrap();

            for (neighbor, _) in &self.graph.adjacency[last_node] {
                if *neighbor < start_node {
                    continue;
                }
                if path.contains(neighbor) {
                    if *neighbor == start_node && path.len() > 2 {
                        let mut cycle = path.clone();
                        cycle.sort_unstable();
                        self.all_cycles.insert(cycle);
                    }
                    continue;
                }

                let mut new_path = path.clone();
                new_path.push(*neighbor);
                queue.push_back(new_path);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element};

    fn build_graph(edges: &[(usize, usize)]) -> ProcessingGraph {
        let mut mg = MolecularGraph::new();
        let max_node = edges
            .iter()
            .flat_map(|(u, v)| vec![*u, *v])
            .max()
            .unwrap_or(0);
        for _ in 0..=max_node {
            mg.add_atom(Element::C);
        }
        for (u, v) in edges {
            mg.add_bond(*u, *v, BondOrder::Single).unwrap();
        }
        ProcessingGraph::new(&mg).unwrap()
    }

    #[test]
    fn test_no_rings() {
        let graph = build_graph(&[(0, 1), (1, 2), (2, 3)]);
        let ring_info = perceive_rings(&graph);
        assert!(ring_info.0.is_empty());
    }

    #[test]
    fn test_single_triangle_ring() {
        let graph = build_graph(&[(0, 1), (1, 2), (2, 0)]);
        let ring_info = perceive_rings(&graph);
        let mut expected_ring: Vec<usize> = vec![0, 1, 2];
        expected_ring.sort();
        assert_eq!(ring_info.0.len(), 1);
        assert!(ring_info.0.contains(&expected_ring));
    }

    #[test]
    fn test_single_square_ring() {
        let graph = build_graph(&[(0, 1), (1, 2), (2, 3), (3, 0)]);
        let ring_info = perceive_rings(&graph);
        let mut expected_ring: Vec<usize> = vec![0, 1, 2, 3];
        expected_ring.sort();
        assert_eq!(ring_info.0.len(), 1);
        assert!(ring_info.0.contains(&expected_ring));
    }

    #[test]
    fn test_two_fused_rings_naphthalene_style() {
        let graph = build_graph(&[(0, 1), (1, 2), (2, 0), (1, 3), (3, 2)]);
        let ring_info = perceive_rings(&graph);

        let mut ring1: Vec<usize> = vec![0, 1, 2];
        ring1.sort();
        let mut ring2: Vec<usize> = vec![1, 2, 3];
        ring2.sort();
        let mut ring3: Vec<usize> = vec![0, 1, 2, 3];
        ring3.sort();

        assert_eq!(ring_info.0.len(), 3);
        assert!(ring_info.0.contains(&ring1));
        assert!(ring_info.0.contains(&ring2));
        assert!(ring_info.0.contains(&ring3));
    }
}
