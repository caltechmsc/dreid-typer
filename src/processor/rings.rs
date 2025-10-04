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
