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

            if let Some(best_rule) = self
                .find_best_matching_rule(atom)
                .filter(|rule| rule.priority > current_priority)
            {
                self.atom_states[atom.id] =
                    Some((best_rule.result_type.clone(), best_rule.priority));
                changes_count += 1;
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
        if conditions.element.is_some_and(|e| e != atom.element) {
            return false;
        }
        if conditions
            .formal_charge
            .is_some_and(|fc| fc != atom.formal_charge)
        {
            return false;
        }
        if conditions.degree.is_some_and(|d| d != atom.degree) {
            return false;
        }
        if conditions
            .is_in_ring
            .is_some_and(|lir| lir != atom.is_in_ring)
        {
            return false;
        }
        if conditions
            .lone_pairs
            .is_some_and(|lp| lp != atom.lone_pairs)
        {
            return false;
        }
        if conditions
            .hybridization
            .is_some_and(|h| h != atom.hybridization)
        {
            return false;
        }
        if conditions
            .is_aromatic
            .is_some_and(|ia| ia != atom.is_aromatic)
        {
            return false;
        }
        if conditions
            .is_anti_aromatic
            .is_some_and(|iaa| iaa != atom.is_anti_aromatic)
        {
            return false;
        }
        if conditions
            .is_resonant
            .is_some_and(|ir| ir != atom.is_in_conjugated_system)
        {
            return false;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::properties::{BondOrder, Element, Hybridization};
    use crate::perception::AnnotatedMolecule;

    fn linear_ethene_like() -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        let c1 = graph.add_atom(Element::C);
        let c2 = graph.add_atom(Element::C);
        let h1 = graph.add_atom(Element::H);
        let h2 = graph.add_atom(Element::H);

        graph
            .add_bond(c1, c2, BondOrder::Double)
            .expect("valid C=C bond");
        graph
            .add_bond(c1, h1, BondOrder::Single)
            .expect("valid C-H bond");
        graph
            .add_bond(c2, h2, BondOrder::Single)
            .expect("valid C-H bond");

        AnnotatedMolecule::new(&graph).expect("graph should be valid")
    }

    fn annotate_sp2_carbons(molecule: &mut AnnotatedMolecule) {
        for atom in &mut molecule.atoms {
            match atom.element {
                Element::C => {
                    atom.hybridization = Hybridization::SP2;
                    atom.degree = 3;
                }
                Element::H => {
                    atom.hybridization = Hybridization::None;
                    atom.degree = 1;
                }
                _ => {}
            }
        }
    }

    fn assign_types_for(
        molecule: &mut AnnotatedMolecule,
        rules: &[Rule],
    ) -> Result<Vec<String>, AssignmentError> {
        annotate_sp2_carbons(molecule);
        assign_types(molecule, rules)
    }

    fn rule(name: &str, priority: i32, result_type: &str, conditions: Conditions) -> Rule {
        Rule {
            name: name.to_string(),
            priority,
            result_type: result_type.to_string(),
            conditions,
        }
    }

    fn condition() -> Conditions {
        Conditions::default()
    }

    #[test]
    fn assigns_simple_sp2_carbons_and_hydrogens() {
        let mut molecule = linear_ethene_like();
        let rules = vec![
            rule(
                "C_SP2",
                10,
                "C_R",
                Conditions {
                    element: Some(Element::C),
                    hybridization: Some(Hybridization::SP2),
                    ..Conditions::default()
                },
            ),
            rule(
                "H_DEFAULT",
                1,
                "H_",
                Conditions {
                    element: Some(Element::H),
                    ..Conditions::default()
                },
            ),
        ];

        let types = assign_types_for(&mut molecule, &rules).expect("typing should succeed");
        assert_eq!(types.len(), 4);
        assert_eq!(types[0], "C_R");
        assert_eq!(types[1], "C_R");
        assert_eq!(types[2], "H_");
        assert_eq!(types[3], "H_");
    }

    #[test]
    fn neighbor_dependent_rules_require_multiple_rounds() {
        let mut molecule = linear_ethene_like();
        let mut neighbor_conditions = condition();
        neighbor_conditions.element = Some(Element::H);

        let mut carbon_neighbor_types = condition();
        carbon_neighbor_types.element = Some(Element::C);
        carbon_neighbor_types.hybridization = Some(Hybridization::SP2);

        let mut hydrogens_require_carbon_type = condition();
        hydrogens_require_carbon_type.element = Some(Element::H);
        hydrogens_require_carbon_type
            .neighbor_types
            .insert("C_R".to_string(), 1);

        let rules = vec![
            rule("Carbon", 5, "C_R", carbon_neighbor_types),
            rule("Hydrogen", 1, "H_", hydrogens_require_carbon_type),
        ];

        let types = assign_types_for(&mut molecule, &rules).expect("typing should converge");
        assert_eq!(types[2], "H_");
        assert_eq!(types[3], "H_");
    }

    #[test]
    fn higher_priority_rule_overrides_lower_one() {
        let mut molecule = linear_ethene_like();
        let mut base_condition = condition();
        base_condition.element = Some(Element::C);

        let mut specific_condition = base_condition.clone();
        specific_condition.hybridization = Some(Hybridization::SP2);

        let rules = vec![
            rule("BaseCarbon", 1, "C_BASE", base_condition),
            rule("SpecificSp2", 10, "C_R", specific_condition),
            rule(
                "Hydrogens",
                1,
                "H_",
                Conditions {
                    element: Some(Element::H),
                    ..Conditions::default()
                },
            ),
        ];

        let types = assign_types_for(&mut molecule, &rules).expect("typing should succeed");
        assert!(types.iter().take(2).all(|t| t == "C_R"));
    }

    #[test]
    fn returns_assignment_error_when_atoms_remain_untyped() {
        let mut molecule = linear_ethene_like();
        let rules = vec![rule(
            "HydrogenOnly",
            1,
            "H_",
            Conditions {
                element: Some(Element::H),
                ..Conditions::default()
            },
        )];

        let err =
            assign_types_for(&mut molecule, &rules).expect_err("carbons should remain untyped");
        assert!(err.untyped_atom_ids.contains(&0));
        assert!(err.untyped_atom_ids.contains(&1));
    }
}
