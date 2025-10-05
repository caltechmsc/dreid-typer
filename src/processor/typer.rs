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
    use crate::processor::process_graph;
    use crate::rules::{Rule, parse_rules};

    fn run_typing_test(
        mg: &MolecularGraph,
        custom_rules: Option<&[Rule]>,
    ) -> Result<Vec<String>, AssignmentError> {
        let pg = process_graph(mg).expect("Annotation failed during test setup");
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
        let c1 = mg.add_atom(Element::C);
        let h1 = mg.add_atom(Element::H);
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c1], "C_3");
        assert_eq!(types[h1], "H_");
    }

    #[test]
    fn test_ethylene_typing() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C);
        let c2 = mg.add_atom(Element::C);
        mg.add_bond(c1, c2, BondOrder::Double).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c1], "C_2");
        assert_eq!(types[c2], "C_2");
    }

    #[test]
    fn test_acetylene_typing() {
        let mut mg = MolecularGraph::new();
        let c1 = mg.add_atom(Element::C);
        let c2 = mg.add_atom(Element::C);
        mg.add_bond(c1, c2, BondOrder::Triple).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[c1], "C_1");
        assert_eq!(types[c2], "C_1");
    }

    #[test]
    fn test_benzene_aromatic_typing() {
        let mut mg = MolecularGraph::new();
        for _ in 0..6 {
            mg.add_atom(Element::C);
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
        let f = mg.add_atom(Element::F);
        let cl = mg.add_atom(Element::Cl);
        let br = mg.add_atom(Element::Br);
        let na = mg.add_atom(Element::Na);

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
        mg.add_atom(Element::C);
        let result = run_typing_test(&mg, None);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().untyped_atom_ids, vec![0]);
    }

    #[test]
    fn test_diatomic_oxygen() {
        let mut mg = MolecularGraph::new();
        let o1 = mg.add_atom(Element::O);
        let o2 = mg.add_atom(Element::O);
        mg.add_bond(o1, o2, BondOrder::Double).unwrap();

        let types = run_typing_test(&mg, None).unwrap();
        assert_eq!(types[o1], "O_2");
        assert_eq!(types[o2], "O_2");
    }

    #[test]
    fn test_acetic_acid_relaxation_logic() {
        let mut mg = MolecularGraph::new();
        let c_me = mg.add_atom(Element::C);
        let c_co = mg.add_atom(Element::C);
        let o_co = mg.add_atom(Element::O);
        let o_oh = mg.add_atom(Element::O);
        let h_oh = mg.add_atom(Element::H);
        let h_me = mg.add_atom(Element::H);

        mg.add_bond(c_me, c_co, BondOrder::Single).unwrap();
        mg.add_bond(c_me, h_me, BondOrder::Single).unwrap();
        mg.add_bond(c_co, o_co, BondOrder::Double).unwrap();
        mg.add_bond(c_co, o_oh, BondOrder::Single).unwrap();
        mg.add_bond(o_oh, h_oh, BondOrder::Single).unwrap();

        let types = run_typing_test(&mg, None).unwrap();

        assert_eq!(types[c_me], "C_3");
        assert_eq!(types[c_co], "C_2");
        assert_eq!(types[o_co], "O_2");
        assert_eq!(types[o_oh], "O_3");
        assert_eq!(types[h_oh], "H_HB");
        assert_eq!(types[h_me], "H_");
    }

    #[test]
    fn test_diborane_special_hydrogen() {
        let mut mg = MolecularGraph::new();
        let b1 = mg.add_atom(Element::B);
        let b2 = mg.add_atom(Element::B);
        let h_bridge1 = mg.add_atom(Element::H);
        let h_bridge2 = mg.add_atom(Element::H);
        let h_term1a = mg.add_atom(Element::H);
        let h_term1b = mg.add_atom(Element::H);
        let h_term2a = mg.add_atom(Element::H);
        let h_term2b = mg.add_atom(Element::H);

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
        mg.add_atom(Element::Lr);

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
        let c1 = mg.add_atom(Element::C);
        mg.add_atom(Element::H);
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
