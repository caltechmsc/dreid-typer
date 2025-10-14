//! Implements the priority-based, iterative rule engine for DREIDING atom type assignment.
//!
//! This module contains the `TyperEngine` which executes the typing phase of the dreid-typer pipeline.
//! It applies a set of rules to assign atom types based on chemical properties, using iteration to
//! handle rules that depend on neighboring atom types.

use super::graph::{AtomView, ProcessingGraph};
use crate::core::Element;
use crate::core::error::AssignmentError;
use crate::rules::{Conditions, Rule};
use std::collections::HashMap;

/// Priority-based iterative rule engine for assigning DREIDING atom types.
///
/// The engine applies rules in order of decreasing priority, iterating until no more changes occur.
/// This handles dependencies between rules where an atom's type depends on its neighbors' types.
/// Priority-based iterative rule engine for assigning DREIDING atom types.
///
/// The engine applies rules in order of decreasing priority, iterating until no more changes occur.
/// This handles dependencies between rules where an atom's type depends on its neighbors' types.
pub(crate) struct TyperEngine<'a> {
    graph: &'a ProcessingGraph,
    rules: Vec<&'a Rule>,
    atom_states: Vec<Option<(String, i32)>>,
    rounds_completed: u32,
}

impl<'a> TyperEngine<'a> {
    /// Creates a new typer engine with sorted rules.
    ///
    /// Rules are sorted by decreasing priority to ensure higher-priority rules are applied first.
    ///
    /// # Arguments
    ///
    /// * `graph` - The processing graph containing atoms to type.
    /// * `rules` - The rules to apply for atom typing.
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

    /// Runs the iterative typing process until convergence.
    ///
    /// Executes rounds of rule application until no atoms change type in a round.
    /// Returns the final atom types or an error if some atoms remain untyped.
    ///
    /// # Returns
    ///
    /// A vector of atom type strings in atom index order.
    ///
    /// # Errors
    ///
    /// Returns `AssignmentError` if any atoms could not be typed after all rounds.
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

    /// Executes one round of rule application across all atoms.
    ///
    /// Attempts to apply the best matching rule to each atom, updating only if a higher-priority
    /// rule is found. Returns the number of atoms that changed type this round.
    ///
    /// # Returns
    ///
    /// The count of atoms that were assigned or reassigned types in this round.
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

    /// Finds the highest-priority rule that matches the given atom.
    ///
    /// Searches through sorted rules to find the first (highest priority) matching rule.
    ///
    /// # Arguments
    ///
    /// * `atom` - The atom to find a matching rule for.
    ///
    /// # Returns
    ///
    /// The best matching rule, or `None` if no rules match.
    fn find_best_matching_rule(&self, atom: &AtomView) -> Option<&'a Rule> {
        self.rules
            .iter()
            .find(|rule| self.match_conditions(atom, &rule.conditions))
            .copied()
    }

    /// Checks if an atom matches all the specified rule conditions.
    ///
    /// Evaluates each condition in the rule against the atom's properties and neighbors.
    ///
    /// # Arguments
    ///
    /// * `atom` - The atom to check against the conditions.
    /// * `conditions` - The conditions that must all be satisfied.
    ///
    /// # Returns
    ///
    /// `true` if all conditions match, `false` otherwise.
    fn match_conditions(&self, atom: &AtomView, conditions: &Conditions) -> bool {
        if let Some(expected) = conditions.element {
            if expected != atom.element {
                return false;
            }
        }
        if let Some(expected) = conditions.formal_charge {
            if expected != atom.formal_charge {
                return false;
            }
        }
        if let Some(expected) = conditions.degree {
            if expected != atom.degree {
                return false;
            }
        }
        if let Some(expected) = conditions.lone_pairs {
            if expected != atom.lone_pairs {
                return false;
            }
        }
        if let Some(expected) = conditions.steric_number {
            if expected != atom.steric_number {
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

        // Check neighbor element counts and neighbor type counts.
        if !self.match_neighbor_elements(atom, &conditions.neighbor_elements) {
            return false;
        }
        if !self.match_neighbor_types(atom, &conditions.neighbor_types) {
            return false;
        }

        true
    }

    /// Verifies that neighbor element counts match the expected distribution.
    ///
    /// Counts the elements of neighboring atoms and compares against expected counts.
    ///
    /// # Arguments
    ///
    /// * `atom` - The central atom whose neighbors to check.
    /// * `expected_neighbors` - Map of element to expected count.
    ///
    /// # Returns
    ///
    /// `true` if neighbor element counts match exactly, `false` otherwise.
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

    /// Verifies that neighbor type counts match the expected distribution.
    ///
    /// Counts the assigned types of neighboring atoms and compares against expected counts.
    /// Only considers neighbors that have already been assigned types.
    ///
    /// # Arguments
    ///
    /// * `atom` - The central atom whose neighbors to check.
    /// * `expected_types` - Map of atom type to expected count.
    ///
    /// # Returns
    ///
    /// `true` if neighbor type counts match exactly, `false` otherwise.
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

/// Assigns DREIDING atom types to all atoms in the processing graph.
///
/// Creates and runs a `TyperEngine` with the given rules to assign types to all atoms.
/// This is the main entry point for the typing phase.
///
/// # Arguments
///
/// * `graph` - The processing graph with annotated atoms.
/// * `rules` - The rules to use for type assignment.
///
/// # Returns
///
/// A vector of atom type strings corresponding to each atom in the graph.
///
/// # Errors
///
/// Returns `AssignmentError` if any atoms could not be assigned types.
pub(crate) fn assign_types(
    graph: &ProcessingGraph,
    rules: &[Rule],
) -> Result<Vec<String>, AssignmentError> {
    TyperEngine::new(graph, rules).run()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element};
    use crate::processor::perceive;
    use crate::rules::{Rule, parse_rules};

    fn run_typing_test(
        mg: &MolecularGraph,
        custom_rules: Option<&[Rule]>,
    ) -> Result<Vec<String>, AssignmentError> {
        let perception = perceive(mg).expect("Annotation failed during test setup");
        let pg = perception.processing_graph;
        let default_rules_storage;
        let rules_to_use = match custom_rules {
            Some(r) => r,
            None => {
                default_rules_storage = crate::rules::get_default_rules().unwrap().to_vec();
                &default_rules_storage
            }
        };
        assign_types(&pg, rules_to_use)
    }

    #[test]
    fn test_methane_typing() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);

        for &hydrogen in &[h1, h2, h3, h4] {
            mg.add_bond(c1, hydrogen, BondOrder::Single).unwrap();
        }

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c1], "C_3");
        for &hydrogen in &[h1, h2, h3, h4] {
            assert_eq!(types[hydrogen], "H_");
        }
    }

    #[test]
    fn test_ethylene_typing() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);
        mg.add_bond(c1, c2, BondOrder::Double).unwrap();
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();
        mg.add_bond(c1, h2, BondOrder::Single).unwrap();
        mg.add_bond(c2, h3, BondOrder::Single).unwrap();
        mg.add_bond(c2, h4, BondOrder::Single).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c1], "C_2");
        assert_eq!(types[c2], "C_2");
        for &hydrogen in &[h1, h2, h3, h4] {
            assert_eq!(types[hydrogen], "H_");
        }
    }

    #[test]
    fn test_acetylene_typing() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        mg.add_bond(c1, c2, BondOrder::Triple).unwrap();
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();
        mg.add_bond(c2, h2, BondOrder::Single).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c1], "C_1");
        assert_eq!(types[c2], "C_1");
        assert_eq!(types[h1], "H_");
        assert_eq!(types[h2], "H_");
    }

    #[test]
    fn test_benzene_aromatic_typing() {
        let mut mg = MolecularGraph::new();
        for _ in 0..6 {
            mg.add_atom(Element::C, 0);
        }
        for i in 0..6 {
            mg.add_bond(i, (i + 1) % 6, BondOrder::Aromatic).unwrap();
        }

        let types = run_typing_test(&mg, None).unwrap();
        for i in 0..6 {
            assert_eq!(types[i], "C_R");
        }
    }

    #[test]
    fn test_halogens_and_ions_typing() {
        let mut mg = MolecularGraph::new();
        let f = mg.add_atom(Element::F, 0);
        let cl = mg.add_atom(Element::Cl, 0);
        let br = mg.add_atom(Element::Br, 0);
        let na = mg.add_atom(Element::Na, 0);

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[f], "F_");
        assert_eq!(types[cl], "Cl_");
        assert_eq!(types[br], "Br_");
        assert_eq!(types[na], "Na");
    }

    #[test]
    fn test_empty_graph() {
        let mg = MolecularGraph::new();
        let types = run_typing_test(&mg, None).unwrap();
        assert!(types.is_empty());
    }

    #[test]
    fn test_single_atom() {
        let mut mg = MolecularGraph::new();
        let c = mg.add_atom(Element::C, 0);
        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c], "C_1");
    }

    #[test]
    fn test_diatomic_oxygen() {
        let mut mg = MolecularGraph::new();
        let o1 = mg.add_atom(Element::O, 0);
        let o2 = mg.add_atom(Element::O, 0);
        mg.add_bond(o1, o2, BondOrder::Double).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[o1], "O_2");
        assert_eq!(types[o2], "O_2");
    }

    #[test]
    fn test_acetic_acid_relaxation_logic() {
        let mut mg = MolecularGraph::new();
        let c_me = mg.add_atom(Element::C, 0);
        let c_co = mg.add_atom(Element::C, 0);
        let o_co = mg.add_atom(Element::O, 0);
        let o_oh = mg.add_atom(Element::O, 0);
        let h_oh = mg.add_atom(Element::H, 0);
        let h_me1 = mg.add_atom(Element::H, 0);
        let h_me2 = mg.add_atom(Element::H, 0);
        let h_me3 = mg.add_atom(Element::H, 0);

        mg.add_bond(c_me, c_co, BondOrder::Single).unwrap();
        mg.add_bond(c_me, h_me1, BondOrder::Single).unwrap();
        mg.add_bond(c_me, h_me2, BondOrder::Single).unwrap();
        mg.add_bond(c_me, h_me3, BondOrder::Single).unwrap();
        mg.add_bond(c_co, o_co, BondOrder::Double).unwrap();
        mg.add_bond(c_co, o_oh, BondOrder::Single).unwrap();
        mg.add_bond(o_oh, h_oh, BondOrder::Single).unwrap();

        let types = run_typing_test(&mg, None).unwrap();

        assert_eq!(types[c_me], "C_3");
        assert_eq!(types[c_co], "C_2");
        assert_eq!(types[o_co], "O_2");
        assert_eq!(types[o_oh], "O_3");
        assert_eq!(types[h_oh], "H_HB");
        for &hydrogen in &[h_me1, h_me2, h_me3] {
            assert_eq!(types[hydrogen], "H_");
        }
    }

    #[test]
    fn test_diborane_special_hydrogen() {
        let mut mg = MolecularGraph::new();
        let b1 = mg.add_atom(Element::B, 0);
        let b2 = mg.add_atom(Element::B, 0);
        let h_bridge1 = mg.add_atom(Element::H, 0);
        let h_bridge2 = mg.add_atom(Element::H, 0);
        let h_term1a = mg.add_atom(Element::H, 0);
        let h_term1b = mg.add_atom(Element::H, 0);
        let h_term2a = mg.add_atom(Element::H, 0);
        let h_term2b = mg.add_atom(Element::H, 0);

        mg.add_bond(b1, h_bridge1, BondOrder::Single).unwrap();
        mg.add_bond(b2, h_bridge1, BondOrder::Single).unwrap();
        mg.add_bond(b1, h_bridge2, BondOrder::Single).unwrap();
        mg.add_bond(b2, h_bridge2, BondOrder::Single).unwrap();

        mg.add_bond(b1, h_term1a, BondOrder::Single).unwrap();
        mg.add_bond(b1, h_term1b, BondOrder::Single).unwrap();
        mg.add_bond(b2, h_term2a, BondOrder::Single).unwrap();
        mg.add_bond(b2, h_term2b, BondOrder::Single).unwrap();

        let types = run_typing_test(&mg, None).unwrap();

        assert_eq!(types[b1], "B_3");
        assert_eq!(types[b2], "B_3");

        assert_eq!(types[h_bridge1], "H_b");
        assert_eq!(types[h_bridge2], "H_b");

        assert_eq!(types[h_term1a], "H_");
        assert_eq!(types[h_term1b], "H_");
        assert_eq!(types[h_term2a], "H_");
        assert_eq!(types[h_term2b], "H_");
    }

    #[test]
    fn test_untypable_element_fails_gracefully() {
        let mut mg = MolecularGraph::new();
        mg.add_atom(Element::Lr, 0);

        let result = run_typing_test(&mg, None);
        assert!(result.is_err());

        if let Err(e) = result {
            assert_eq!(e.untyped_atom_ids, vec![0]);
            assert!(e.rounds_completed > 0);
        } else {
            panic!("Expected an error but got Ok");
        }
    }

    #[test]
    fn test_priority_logic_chooses_higher_priority_rule() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C, 0);
        mg.add_atom(Element::H, 0);
        mg.add_bond(c1, 1, BondOrder::Single).unwrap();

        let custom_rules_toml = r#"
            [[rule]]
            name = "C_Special_High_P"
            priority = 500
            type = "C_SP"
            conditions = { element = "C", degree = 1 }

            [[rule]]
            name = "C_Tetrahedral_Default"
            priority = 100
            type = "C_3"
            conditions = { element = "C", hybridization = "SP3" }

            [[rule]]
            name = "H_Standard_Default"
            priority = 1
            type = "H_"
            conditions = { element = "H" }
        "#;
        let rules = parse_rules(custom_rules_toml).unwrap();

        let types = run_typing_test(&mg, Some(&rules)).unwrap();

        assert_eq!(types[c1], "C_SP");
    }
}
