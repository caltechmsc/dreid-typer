use super::rules::{Conditions, Rule};
use crate::core::error::AssignmentError;
use crate::core::properties::Element;
use crate::perception::{AnnotatedAtom, AnnotatedMolecule};
use std::collections::HashMap;

pub fn assign_types(
    molecule: &AnnotatedMolecule,
    rules: &[Rule],
) -> Result<Vec<String>, AssignmentError> {
    let mut engine = TyperEngine::new(molecule, rules);
    engine.run()
}

struct TyperEngine<'a> {
    molecule: &'a AnnotatedMolecule,
    sorted_rules: Vec<&'a Rule>,
    atom_states: Vec<Option<(String, i32)>>,
}

impl<'a> TyperEngine<'a> {
    fn new(molecule: &'a AnnotatedMolecule, rules: &'a [Rule]) -> Self {
        let mut sorted_rules: Vec<&'a Rule> = rules.iter().collect();
        sorted_rules.sort_by(|a, b| {
            b.priority
                .cmp(&a.priority)
                .then_with(|| a.name.cmp(&b.name))
        });

        Self {
            molecule,
            sorted_rules,
            atom_states: vec![None; molecule.atoms.len()],
        }
    }

    fn run(&mut self) -> Result<Vec<String>, AssignmentError> {
        let mut rounds = 0;
        const MAX_ROUNDS: u32 = 100;

        loop {
            rounds += 1;
            if rounds > MAX_ROUNDS {
                return Err(self.build_error(rounds));
            }

            let changes = self.run_single_round();
            if changes == 0 {
                break;
            }
        }

        let mut final_types = Vec::with_capacity(self.molecule.atoms.len());
        let mut untyped_ids = Vec::new();

        for (i, state) in self.atom_states.iter().enumerate() {
            if let Some((type_name, _)) = state {
                final_types.push(type_name.clone());
            } else {
                untyped_ids.push(i);
            }
        }

        if untyped_ids.is_empty() {
            Ok(final_types)
        } else {
            Err(self.build_error(rounds))
        }
    }

    fn run_single_round(&mut self) -> usize {
        let mut changes_count = 0;

        for atom in &self.molecule.atoms {
            let current_priority = self.atom_states[atom.id].as_ref().map_or(-1, |(_, p)| *p);

            if let Some(best_rule) = self.find_best_matching_rule(atom) {
                if best_rule.priority > current_priority {
                    self.atom_states[atom.id] =
                        Some((best_rule.result_type.clone(), best_rule.priority));
                    changes_count += 1;
                }
            }
        }
        changes_count
    }

    fn find_best_matching_rule(&self, atom: &AnnotatedAtom) -> Option<&'a Rule> {
        self.sorted_rules
            .iter()
            .find(|rule| self.match_conditions(atom, &rule.conditions))
            .copied()
    }

    fn match_conditions(&self, atom: &AnnotatedAtom, conditions: &Conditions) -> bool {
        if let Some(e) = conditions.element {
            if e != atom.element {
                return false;
            }
        }
        if let Some(fc) = conditions.formal_charge {
            if fc != atom.formal_charge {
                return false;
            }
        }
        if let Some(d) = conditions.degree {
            if d != atom.degree {
                return false;
            }
        }
        if let Some(lir) = conditions.is_in_ring {
            if lir != atom.is_in_ring {
                return false;
            }
        }
        if let Some(lp) = conditions.lone_pairs {
            if lp != atom.lone_pairs {
                return false;
            }
        }
        if let Some(h) = conditions.hybridization {
            if h != atom.hybridization {
                return false;
            }
        }
        if let Some(ia) = conditions.is_aromatic {
            if ia != atom.is_aromatic {
                return false;
            }
        }
        if let Some(iaa) = conditions.is_anti_aromatic {
            if iaa != atom.is_anti_aromatic {
                return false;
            }
        }
        if let Some(ir) = conditions.is_resonant {
            if ir != atom.is_in_conjugated_system {
                return false;
            }
        }

        if !conditions.neighbor_elements.is_empty()
            && !self.match_neighbor_elements(atom, &conditions.neighbor_elements)
        {
            return false;
        }
        if !conditions.neighbor_types.is_empty()
            && !self.match_neighbor_types(atom, &conditions.neighbor_types)
        {
            return false;
        }

        true
    }

    fn match_neighbor_elements(
        &self,
        atom: &AnnotatedAtom,
        expected: &HashMap<Element, u8>,
    ) -> bool {
        let mut actual_counts: HashMap<Element, u8> = HashMap::new();
        for &(neighbor_id, _) in &self.molecule.adjacency[atom.id] {
            let neighbor_element = self.molecule.atoms[neighbor_id].element;
            *actual_counts.entry(neighbor_element).or_default() += 1;
        }

        expected
            .iter()
            .all(|(element, &count)| actual_counts.get(element).copied().unwrap_or(0) == count)
    }

    fn match_neighbor_types(&self, atom: &AnnotatedAtom, expected: &HashMap<String, u8>) -> bool {
        let mut actual_counts: HashMap<&str, u8> = HashMap::new();
        for &(neighbor_id, _) in &self.molecule.adjacency[atom.id] {
            if let Some((type_name, _)) = &self.atom_states[neighbor_id] {
                *actual_counts.entry(type_name).or_default() += 1;
            }
        }

        expected.iter().all(|(type_name, &count)| {
            actual_counts.get(type_name.as_str()).copied().unwrap_or(0) == count
        })
    }

    fn build_error(&self, rounds_completed: u32) -> AssignmentError {
        let untyped_atom_ids = self
            .atom_states
            .iter()
            .enumerate()
            .filter(|(_, state)| state.is_none())
            .map(|(i, _)| i)
            .collect();
        AssignmentError {
            untyped_atom_ids,
            rounds_completed,
        }
    }
}
