use super::graph::{AtomView, ProcessingGraph};
use crate::core::Element;
use crate::core::error::AssignmentError;
use crate::rules::{Conditions, Rule};
use std::collections::HashMap;

pub(crate) struct TyperEngine<'a> {
    graph: &'a ProcessingGraph,
    rules: Vec<&'a Rule>,
    atom_states: Vec<Option<(String, i32)>>,
    rounds_completed: u32,
}

impl<'a> TyperEngine<'a> {
    pub fn new(graph: &'a ProcessingGraph, rules: &'a [Rule]) -> Self {
        let num_atoms = graph.atoms.len();

        let mut sorted_rules: Vec<&'a Rule> = rules.iter().collect();
        sorted_rules.sort_by(|a, b| {
            b.priority
                .cmp(&a.priority)
                .then_with(|| a.name.cmp(&b.name))
        });

        Self {
            graph,
            rules: sorted_rules,
            atom_states: vec![None; num_atoms],
            rounds_completed: 0,
        }
    }

    pub fn run(mut self) -> Result<Vec<String>, AssignmentError> {
        loop {
            self.rounds_completed += 1;
            let changes_count = self.run_single_round();
            if changes_count == 0 {
                break;
            }
        }

        let mut final_types = Vec::with_capacity(self.graph.atoms.len());
        let mut untyped_ids = vec![];
        for (i, state) in self.atom_states.into_iter().enumerate() {
            if let Some((type_name, _)) = state {
                final_types.push(type_name);
            } else {
                untyped_ids.push(i);
            }
        }

        if !untyped_ids.is_empty() {
            Err(AssignmentError {
                untyped_atom_ids: untyped_ids,
                rounds_completed: self.rounds_completed,
            })
        } else {
            Ok(final_types)
        }
    }

    fn run_single_round(&mut self) -> usize {
        let mut changes_this_round = vec![];

        for atom_view in &self.graph.atoms {
            if let Some(best_rule) = self.find_best_matching_rule(atom_view) {
                let current_priority = self.atom_states[atom_view.id]
                    .as_ref()
                    .map_or(-1, |(_, p)| *p);

                if best_rule.priority > current_priority {
                    changes_this_round.push((
                        atom_view.id,
                        (best_rule.result_type.clone(), best_rule.priority),
                    ));
                }
            }
        }

        let changes_count = changes_this_round.len();
        if changes_count > 0 {
            for (atom_id, new_state) in changes_this_round {
                self.atom_states[atom_id] = Some(new_state);
            }
        }

        changes_count
    }

    fn find_best_matching_rule(&self, atom: &AtomView) -> Option<&'a Rule> {
        self.rules
            .iter()
            .find(|rule| self.match_conditions(atom, &rule.conditions))
            .copied()
    }

    fn match_conditions(&self, atom: &AtomView, conditions: &Conditions) -> bool {
        if let Some(expected) = conditions.element {
            if expected != atom.element {
                return false;
            }
        }
        if let Some(expected) = conditions.degree {
            if expected != atom.degree {
                return false;
            }
        }
        if let Some(expected) = conditions.hybridization {
            if expected != atom.hybridization {
                return false;
            }
        }
        if let Some(expected) = conditions.is_in_ring {
            if expected != atom.is_in_ring {
                return false;
            }
        }
        if let Some(expected) = conditions.is_aromatic {
            if expected != atom.is_aromatic {
                return false;
            }
        }
        if let Some(expected) = conditions.smallest_ring_size {
            if Some(expected) != atom.smallest_ring_size {
                return false;
            }
        }

        if !self.match_neighbor_elements(atom, &conditions.neighbor_elements) {
            return false;
        }
        if !self.match_neighbor_types(atom, &conditions.neighbor_types) {
            return false;
        }

        true
    }

    fn match_neighbor_elements(
        &self,
        atom: &AtomView,
        expected_neighbors: &HashMap<Element, u8>,
    ) -> bool {
        if expected_neighbors.is_empty() {
            return true;
        }

        let mut actual_counts: HashMap<Element, u8> = HashMap::new();
        for (neighbor_id, _) in &self.graph.adjacency[atom.id] {
            let neighbor_element = self.graph.atoms[*neighbor_id].element;
            *actual_counts.entry(neighbor_element).or_insert(0) += 1;
        }

        expected_neighbors
            .iter()
            .all(|(element, &count)| actual_counts.get(element).copied().unwrap_or(0) == count)
    }

    fn match_neighbor_types(&self, atom: &AtomView, expected_types: &HashMap<String, u8>) -> bool {
        if expected_types.is_empty() {
            return true;
        }

        let mut actual_counts: HashMap<&str, u8> = HashMap::new();
        for (neighbor_id, _) in &self.graph.adjacency[atom.id] {
            if let Some((type_name, _)) = &self.atom_states[*neighbor_id] {
                *actual_counts.entry(type_name).or_insert(0) += 1;
            }
        }

        expected_types.iter().all(|(type_name, &count)| {
            actual_counts.get(type_name.as_str()).copied().unwrap_or(0) == count
        })
    }
}
