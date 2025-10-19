//! Functional group template matching for chemical perception.
//!
//! This module implements an expert system of functional group templates that
//! override default perception algorithms for specific chemical environments.
//! Templates are applied during the perception phase to correctly identify
//! hybridization states, aromaticity, and other properties in complex molecules.

use crate::core::{BondOrder, Element, Hybridization};
use crate::processor::graph::{AtomView, PerceptionSource, ProcessingGraph};
use std::collections::HashMap;
use std::sync::LazyLock;

/// Applies functional group templates to override default perception results.
///
/// This function iterates through predefined templates in order of decreasing
/// size to avoid conflicts, finding non-overlapping matches and applying
/// the associated actions to set correct chemical properties.
pub(crate) fn apply_functional_group_templates(
    graph: &mut ProcessingGraph,
) -> Result<(), crate::core::error::AnnotationError> {
    let mut sorted_templates = TEMPLATES.clone();
    sorted_templates.sort_by(|a, b| b.nodes.len().cmp(&a.nodes.len()));

    for template in sorted_templates.iter() {
        let matches = find_non_overlapping_matches(graph, template);
        for a_match in matches {
            apply_actions(graph, &a_match, &template.actions);
        }
    }
    Ok(())
}

/// Actions that can be applied to atoms when a template matches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Action {
    /// Sets the chemical state of an atom to a predefined configuration.
    SetState(ChemicalState),
}

/// Predefined chemical states that templates can assign to atoms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChemicalState {
    /// Aromatic state with resonant hybridization.
    Aromatic,
    /// Trigonal planar geometry with SP2 hybridization.
    TrigonalPlanar,
    /// Tetrahedral geometry with SP3 hybridization.
    Tetrahedral,
}

/// A node in a functional group query pattern.
///
/// Each query node defines a label and a predicate function that tests
/// whether an atom matches the required properties for that position in the template.
#[derive(Clone)]
struct QueryNode {
    /// The label used to reference this node in edges and actions.
    label: &'static str,
    /// Predicate function that tests if an atom matches this node's requirements.
    predicate: fn(&AtomView) -> bool,
}

/// An edge in a functional group query pattern.
///
/// Query edges define connectivity requirements between labeled nodes,
/// including bond order constraints.
#[derive(Clone)]
struct QueryEdge {
    /// The labels of the two nodes connected by this edge.
    labels: (&'static str, &'static str),
    /// Predicate function that tests if the bond order matches requirements.
    predicate: fn(BondOrder) -> bool,
}

/// A complete functional group template for pattern matching.
///
/// Templates define a subgraph pattern to match against molecules,
/// consisting of nodes, edges, and actions to apply when the pattern is found.
#[derive(Clone)]
struct FunctionalGroupTemplate {
    /// Descriptive name of the functional group.
    #[allow(dead_code)]
    name: &'static str,
    /// The nodes (atoms) in the template pattern.
    nodes: Vec<QueryNode>,
    /// The edges (bonds) connecting the nodes.
    edges: Vec<QueryEdge>,
    /// Actions to apply to matched atoms.
    actions: HashMap<&'static str, Action>,
}

/// A mapping from template node labels to graph atom IDs.
type Match = HashMap<&'static str, usize>;

/// Finds all non-overlapping matches of a template in the graph.
///
/// This function ensures that each atom is used in at most one template match
/// and prefers larger templates by processing them first. Atoms already
/// modified by templates are skipped.
fn find_non_overlapping_matches(
    graph: &ProcessingGraph,
    template: &FunctionalGroupTemplate,
) -> Vec<Match> {
    let mut all_matches = Vec::new();
    let mut matched_graph_atoms = vec![false; graph.atoms.len()];

    for i in 0..graph.atoms.len() {
        if matched_graph_atoms[i]
            || matches!(
                graph.atoms[i].perception_source,
                Some(PerceptionSource::Template)
            )
        {
            continue;
        }

        let mut current_match = HashMap::new();

        if find_first_match_recursive(graph, template, &mut current_match, 0, &matched_graph_atoms)
        {
            for &atom_id in current_match.values() {
                matched_graph_atoms[atom_id] = true;
            }
            all_matches.push(current_match);
        }
    }
    all_matches
}

/// Recursively finds the first complete match for a template starting from a given position.
///
/// This backtracking algorithm attempts to match template nodes to graph atoms
/// in order, ensuring all edge constraints are satisfied.
fn find_first_match_recursive(
    graph: &ProcessingGraph,
    template: &FunctionalGroupTemplate,
    current_match: &mut Match,
    query_node_idx: usize,
    used_atoms: &[bool],
) -> bool {
    if query_node_idx == template.nodes.len() {
        return verify_edges(graph, template, current_match);
    }

    let query_node = &template.nodes[query_node_idx];

    for atom in &graph.atoms {
        if used_atoms[atom.id]
            || current_match.values().any(|&id| id == atom.id)
            || matches!(atom.perception_source, Some(PerceptionSource::Template))
        {
            continue;
        }

        if (query_node.predicate)(atom) {
            current_match.insert(query_node.label, atom.id);

            if find_first_match_recursive(
                graph,
                template,
                current_match,
                query_node_idx + 1,
                used_atoms,
            ) {
                return true;
            }

            current_match.remove(query_node.label);
        }
    }

    false
}

/// Verifies that all edge constraints in a template are satisfied by the current match.
///
/// Checks each required edge between matched atoms to ensure bond orders match
/// the template specifications.
fn verify_edges(
    graph: &ProcessingGraph,
    template: &FunctionalGroupTemplate,
    a_match: &Match,
) -> bool {
    template.edges.iter().all(|edge| {
        let id1 = a_match[edge.labels.0];
        let id2 = a_match[edge.labels.1];
        graph.adjacency[id1]
            .iter()
            .any(|(neighbor_id, order)| *neighbor_id == id2 && (edge.predicate)(*order))
    })
}

/// Applies the actions defined in a template to the matched atoms.
///
/// Updates the perception source and applies state changes to each matched atom
/// according to the template's action specifications.
fn apply_actions(
    graph: &mut ProcessingGraph,
    a_match: &Match,
    actions: &HashMap<&'static str, Action>,
) {
    for (label, atom_id) in a_match.iter() {
        if let Some(action) = actions.get(label) {
            let atom = &mut graph.atoms[*atom_id];
            atom.perception_source = Some(PerceptionSource::Template);

            match action {
                Action::SetState(state) => apply_state_change(atom, *state),
            }
        }
    }
}

/// Applies a predefined chemical state to an atom.
///
/// Sets hybridization, steric number, and aromaticity properties according
/// to the specified chemical state configuration.
fn apply_state_change(atom: &mut AtomView, state: ChemicalState) {
    match state {
        ChemicalState::Aromatic => {
            atom.is_aromatic = true;
            atom.hybridization = Hybridization::Resonant;
            atom.steric_number = 3;
        }
        ChemicalState::TrigonalPlanar => {
            atom.is_aromatic = false;
            atom.hybridization = Hybridization::SP2;
            atom.steric_number = 3;
        }
        ChemicalState::Tetrahedral => {
            atom.is_aromatic = false;
            atom.hybridization = Hybridization::SP3;
            atom.steric_number = 4;
        }
    }
}

/// Lazily-loaded collection of all functional group templates.
///
/// Templates are defined at compile time and sorted by size during application
/// to ensure larger, more specific patterns are matched before smaller ones.
static TEMPLATES: LazyLock<Vec<FunctionalGroupTemplate>> = LazyLock::new(define_templates);

/// Defines all functional group templates used in perception.
///
/// This function creates templates for common functional groups that require
/// special handling beyond general perception algorithms, such as resonance
/// systems and charged species.
fn define_templates() -> Vec<FunctionalGroupTemplate> {
    vec![
        // --- 1. Guanidinium ---
        FunctionalGroupTemplate {
            name: "Guanidinium",
            nodes: vec![
                QueryNode {
                    label: "C",
                    predicate: |a| a.element == Element::C && a.degree == 3,
                },
                QueryNode {
                    label: "N1",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "N2",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "N3",
                    predicate: |a| a.element == Element::N,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("C", "N1"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C", "N2"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C", "N3"),
                    predicate: |_| true,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("C", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N1", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N2", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N3", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 2. Amide ---
        FunctionalGroupTemplate {
            name: "Amide",
            nodes: vec![
                QueryNode {
                    label: "C",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "O",
                    predicate: |a| a.element == Element::O,
                },
                QueryNode {
                    label: "N",
                    predicate: |a| a.element == Element::N,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("C", "O"),
                    predicate: |o| o == BondOrder::Double,
                },
                QueryEdge {
                    labels: ("C", "N"),
                    predicate: |o| o == BondOrder::Single,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("C", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 3a. Thioamide (bis-amide) ---
        FunctionalGroupTemplate {
            name: "ThioamideBis",
            nodes: vec![
                QueryNode {
                    label: "C",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "S",
                    predicate: |a| a.element == Element::S,
                },
                QueryNode {
                    label: "N1",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "N2",
                    predicate: |a| a.element == Element::N,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("C", "S"),
                    predicate: |o| o == BondOrder::Double,
                },
                QueryEdge {
                    labels: ("C", "N1"),
                    predicate: |o| o == BondOrder::Single,
                },
                QueryEdge {
                    labels: ("C", "N2"),
                    predicate: |o| o == BondOrder::Single,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("C", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("S", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N1", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N2", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 3b. Thioamide ---
        FunctionalGroupTemplate {
            name: "Thioamide",
            nodes: vec![
                QueryNode {
                    label: "C",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "S",
                    predicate: |a| a.element == Element::S,
                },
                QueryNode {
                    label: "N",
                    predicate: |a| a.element == Element::N,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("C", "S"),
                    predicate: |o| o == BondOrder::Double,
                },
                QueryEdge {
                    labels: ("C", "N"),
                    predicate: |o| o == BondOrder::Single,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("C", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("S", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("N", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 4. Carboxylate ---
        FunctionalGroupTemplate {
            name: "Carboxylate",
            nodes: vec![
                QueryNode {
                    label: "C",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "O1",
                    predicate: |a| a.element == Element::O,
                },
                QueryNode {
                    label: "O2",
                    predicate: |a| a.element == Element::O && a.lone_pairs == 3,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("C", "O1"),
                    predicate: |o| o == BondOrder::Double,
                },
                QueryEdge {
                    labels: ("C", "O2"),
                    predicate: |o| o == BondOrder::Single,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("C", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("O1", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("O2", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 5. Nitro ---
        FunctionalGroupTemplate {
            name: "Nitro",
            nodes: vec![
                QueryNode {
                    label: "N",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "O1",
                    predicate: |a| a.element == Element::O,
                },
                QueryNode {
                    label: "O2",
                    predicate: |a| a.element == Element::O && a.lone_pairs == 3,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("N", "O1"),
                    predicate: |o| o == BondOrder::Double,
                },
                QueryEdge {
                    labels: ("N", "O2"),
                    predicate: |o| o == BondOrder::Single,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("N", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("O1", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("O2", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 6. Phenol/Enol ---
        FunctionalGroupTemplate {
            name: "Phenol/Enol",
            nodes: vec![
                QueryNode {
                    label: "O",
                    predicate: |a| a.element == Element::O && a.degree == 2,
                },
                QueryNode {
                    label: "C1",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "C2",
                    predicate: |a| a.element == Element::C,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("O", "C1"),
                    predicate: |o| o == BondOrder::Single,
                },
                QueryEdge {
                    labels: ("C1", "C2"),
                    predicate: |o| o == BondOrder::Double || o == BondOrder::Aromatic,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("O", Action::SetState(ChemicalState::Aromatic));
                map
            },
        },
        // --- 7. Phosphate ---
        FunctionalGroupTemplate {
            name: "Phosphate",
            nodes: vec![
                QueryNode {
                    label: "P",
                    predicate: |a| a.element == Element::P,
                },
                QueryNode {
                    label: "O_double",
                    predicate: |a| a.element == Element::O && a.formal_charge == 0,
                },
                QueryNode {
                    label: "O_single",
                    predicate: |a| a.element == Element::O && a.formal_charge == -1,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("P", "O_double"),
                    predicate: |o| o == BondOrder::Double,
                },
                QueryEdge {
                    labels: ("P", "O_single"),
                    predicate: |o| o == BondOrder::Single,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                map.insert("P", Action::SetState(ChemicalState::Tetrahedral));
                map.insert("O_double", Action::SetState(ChemicalState::TrigonalPlanar));
                map.insert("O_single", Action::SetState(ChemicalState::TrigonalPlanar));
                map
            },
        },
        // --- 8. Purine Skeleton ---
        FunctionalGroupTemplate {
            name: "PurineSkeleton",
            nodes: vec![
                QueryNode {
                    label: "N1",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "C2",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "N3",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "C4",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "C5",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "C6",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "N7",
                    predicate: |a| a.element == Element::N,
                },
                QueryNode {
                    label: "C8",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "N9",
                    predicate: |a| a.element == Element::N,
                },
            ],
            edges: vec![
                QueryEdge {
                    labels: ("N1", "C2"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C2", "N3"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("N3", "C4"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C4", "C5"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C5", "C6"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C6", "N1"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C4", "N9"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("N9", "C8"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("C8", "N7"),
                    predicate: |_| true,
                },
                QueryEdge {
                    labels: ("N7", "C5"),
                    predicate: |_| true,
                },
            ],
            actions: {
                let mut map = HashMap::new();
                for label in ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"] {
                    map.insert(label, Action::SetState(ChemicalState::Aromatic));
                }
                map
            },
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::super::graph::PerceptionSource;
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element, Hybridization};

    #[test]
    fn guanidinium_template_sets_planar_states_on_all_atoms() {
        let mut mg = MolecularGraph::new();
        let carbon = mg.add_atom(Element::C, 0);
        let n1 = mg.add_atom(Element::N, 1);
        let n2 = mg.add_atom(Element::N, 1);
        let n3 = mg.add_atom(Element::N, 1);

        mg.add_bond(carbon, n1, BondOrder::Single).expect("C-N1");
        mg.add_bond(carbon, n2, BondOrder::Single).expect("C-N2");
        mg.add_bond(carbon, n3, BondOrder::Single).expect("C-N3");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        for &idx in [carbon, n1, n2, n3].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
        }
    }

    #[test]
    fn amide_template_sets_planar_carbon_and_nitrogen() {
        let mut mg = MolecularGraph::new();
        let carbon = mg.add_atom(Element::C, 0);
        let oxygen = mg.add_atom(Element::O, 0);
        let nitrogen = mg.add_atom(Element::N, 0);

        mg.add_bond(carbon, oxygen, BondOrder::Double).expect("C=O");
        mg.add_bond(carbon, nitrogen, BondOrder::Single)
            .expect("C-N");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        for &idx in [carbon, nitrogen].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
        }

        let oxygen_atom = &pg.atoms[oxygen];
        assert_ne!(
            oxygen_atom.perception_source,
            Some(PerceptionSource::Template)
        );
    }

    #[test]
    fn thioamide_bis_template_sets_trigonal_planar_states() {
        let mut mg = MolecularGraph::new();
        let carbon = mg.add_atom(Element::C, 0);
        let sulfur = mg.add_atom(Element::S, 0);
        let nitrogen1 = mg.add_atom(Element::N, 0);
        let nitrogen2 = mg.add_atom(Element::N, 0);

        mg.add_bond(carbon, sulfur, BondOrder::Double).expect("C=S");
        mg.add_bond(carbon, nitrogen1, BondOrder::Single)
            .expect("C-N1");
        mg.add_bond(carbon, nitrogen2, BondOrder::Single)
            .expect("C-N2");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        for &idx in [carbon, sulfur, nitrogen1, nitrogen2].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
        }
    }

    #[test]
    fn thioamide_template_sets_trigonal_planar_states() {
        let mut mg = MolecularGraph::new();
        let carbon = mg.add_atom(Element::C, 0);
        let sulfur = mg.add_atom(Element::S, 0);
        let nitrogen = mg.add_atom(Element::N, 0);

        mg.add_bond(carbon, sulfur, BondOrder::Double).expect("C=S");
        mg.add_bond(carbon, nitrogen, BondOrder::Single)
            .expect("C-N");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        for &idx in [carbon, sulfur, nitrogen].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
        }
    }

    #[test]
    fn carboxylate_template_sets_planar_states_for_carbon_and_oxygens() {
        let mut mg = MolecularGraph::new();
        let carbon = mg.add_atom(Element::C, 0);
        let oxygen_double = mg.add_atom(Element::O, 0);
        let oxygen_single = mg.add_atom(Element::O, -1);

        mg.add_bond(carbon, oxygen_double, BondOrder::Double)
            .expect("C=O");
        mg.add_bond(carbon, oxygen_single, BondOrder::Single)
            .expect("C-O");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        for &idx in [carbon, oxygen_double, oxygen_single].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
        }
    }

    #[test]
    fn phenol_enol_template_marks_oxygen_aromatic() {
        let mut mg = MolecularGraph::new();
        let oxygen = mg.add_atom(Element::O, 0);
        let c1 = mg.add_atom(Element::C, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let hydrogen = mg.add_atom(Element::H, 0);

        mg.add_bond(oxygen, c1, BondOrder::Single).expect("O-C1");
        mg.add_bond(c1, c2, BondOrder::Aromatic).expect("C1=C2");
        mg.add_bond(oxygen, hydrogen, BondOrder::Single)
            .expect("O-H");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        let oxygen_atom = &pg.atoms[oxygen];
        assert!(oxygen_atom.is_aromatic);
        assert_eq!(oxygen_atom.hybridization, Hybridization::Resonant);
        assert_eq!(oxygen_atom.steric_number, 3);
        assert_eq!(
            oxygen_atom.perception_source,
            Some(PerceptionSource::Template)
        );
    }

    #[test]
    fn nitro_template_sets_trigonal_planar_states() {
        let mut mg = MolecularGraph::new();
        let nitrogen = mg.add_atom(Element::N, 1);
        let oxygen_neutral = mg.add_atom(Element::O, 0);
        let oxygen_anion = mg.add_atom(Element::O, -1);
        mg.add_bond(nitrogen, oxygen_neutral, BondOrder::Double)
            .expect("double bond");
        mg.add_bond(nitrogen, oxygen_anion, BondOrder::Single)
            .expect("single bond");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        for &idx in [nitrogen, oxygen_neutral, oxygen_anion].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
        }
    }

    #[test]
    fn phosphate_template_sets_expected_states() {
        let mut mg = MolecularGraph::new();
        let phosphorus = mg.add_atom(Element::P, 0);
        let oxygen_double = mg.add_atom(Element::O, 0);
        let oxygen_single = mg.add_atom(Element::O, -1);

        mg.add_bond(phosphorus, oxygen_double, BondOrder::Double)
            .expect("P=O");
        mg.add_bond(phosphorus, oxygen_single, BondOrder::Single)
            .expect("P-O");

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph creation failed");
        apply_functional_group_templates(&mut pg).expect("template application failed");

        let phosphorus_atom = &pg.atoms[phosphorus];
        assert_eq!(
            phosphorus_atom.perception_source,
            Some(PerceptionSource::Template)
        );
        assert_eq!(phosphorus_atom.hybridization, Hybridization::SP3);
        assert_eq!(phosphorus_atom.steric_number, 4);

        for &idx in [oxygen_double, oxygen_single].iter() {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert_eq!(atom.steric_number, 3);
        }
    }

    #[test]
    fn purine_skeleton_template_marks_fused_ring_atoms_aromatic() {
        let mut mg = MolecularGraph::new();
        let n1 = mg.add_atom(Element::N, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let n3 = mg.add_atom(Element::N, 0);
        let c4 = mg.add_atom(Element::C, 0);
        let c5 = mg.add_atom(Element::C, 0);
        let c6 = mg.add_atom(Element::C, 0);
        let n7 = mg.add_atom(Element::N, 0);
        let c8 = mg.add_atom(Element::C, 0);
        let n9 = mg.add_atom(Element::N, 0);

        mg.add_bond(n1, c2, BondOrder::Aromatic).unwrap();
        mg.add_bond(c2, n3, BondOrder::Aromatic).unwrap();
        mg.add_bond(n3, c4, BondOrder::Aromatic).unwrap();
        mg.add_bond(c4, c5, BondOrder::Aromatic).unwrap();
        mg.add_bond(c5, c6, BondOrder::Aromatic).unwrap();
        mg.add_bond(c6, n1, BondOrder::Aromatic).unwrap();

        mg.add_bond(c4, n9, BondOrder::Aromatic).unwrap();
        mg.add_bond(n9, c8, BondOrder::Aromatic).unwrap();
        mg.add_bond(c8, n7, BondOrder::Aromatic).unwrap();
        mg.add_bond(n7, c5, BondOrder::Aromatic).unwrap();

        let mut pg = ProcessingGraph::new(&mg).expect("processing graph");
        apply_functional_group_templates(&mut pg).expect("apply templates");

        for &idx in &[n1, c2, n3, c4, c5, c6, n7, c8, n9] {
            let atom = &pg.atoms[idx];
            assert_eq!(atom.perception_source, Some(PerceptionSource::Template));
            assert!(atom.is_aromatic, "atom {} not aromatic", idx);
            assert_eq!(atom.hybridization, Hybridization::Resonant);
        }
    }
}
