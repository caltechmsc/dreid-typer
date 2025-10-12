use crate::core::{BondOrder, Element, Hybridization};
use crate::processor::graph::{AtomView, PerceptionSource, ProcessingGraph};
use std::collections::HashMap;
use std::sync::LazyLock;

pub(crate) fn apply_functional_group_templates(
    graph: &mut ProcessingGraph,
) -> Result<(), crate::core::error::AnnotationError> {
    for template in TEMPLATES.iter() {
        let matches = find_non_overlapping_matches(graph, template);
        for a_match in matches {
            apply_actions(graph, &a_match, &template.actions);
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Action {
    SetState(ChemicalState),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChemicalState {
    Aromatic,
    TrigonalPlanar,
    Tetrahedral,
}

#[derive(Clone)]
struct QueryNode {
    label: &'static str,
    predicate: fn(&AtomView) -> bool,
}

#[derive(Clone)]
struct QueryEdge {
    labels: (&'static str, &'static str),
    predicate: fn(BondOrder) -> bool,
}

#[derive(Clone)]
struct FunctionalGroupTemplate {
    name: &'static str,
    nodes: Vec<QueryNode>,
    edges: Vec<QueryEdge>,
    actions: HashMap<&'static str, Action>,
}

type Match = HashMap<&'static str, usize>;

fn find_non_overlapping_matches(
    graph: &ProcessingGraph,
    template: &FunctionalGroupTemplate,
) -> Vec<Match> {
    let mut all_matches = Vec::new();
    let mut matched_graph_atoms = vec![false; graph.atoms.len()];

    for i in 0..graph.atoms.len() {
        if matched_graph_atoms[i] {
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
        if used_atoms[atom.id] || current_match.values().any(|&id| id == atom.id) {
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

static TEMPLATES: LazyLock<Vec<FunctionalGroupTemplate>> = LazyLock::new(define_templates);

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
        // --- 3. Carboxylate ---
        FunctionalGroupTemplate {
            name: "Carboxylate",
            nodes: vec![
                QueryNode {
                    label: "C",
                    predicate: |a| a.element == Element::C,
                },
                QueryNode {
                    label: "O1",
                    predicate: |a| a.element == Element::O && a.formal_charge == 0,
                },
                QueryNode {
                    label: "O2",
                    predicate: |a| a.element == Element::O && a.formal_charge == -1,
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
        // --- 4. Nitro ---
        FunctionalGroupTemplate {
            name: "Nitro",
            nodes: vec![
                QueryNode {
                    label: "N",
                    predicate: |a| a.element == Element::N && a.formal_charge == 1,
                },
                QueryNode {
                    label: "O1",
                    predicate: |a| a.element == Element::O && a.formal_charge == 0,
                },
                QueryNode {
                    label: "O2",
                    predicate: |a| a.element == Element::O && a.formal_charge == -1,
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
        // --- 5. Phenol/Enol ---
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
    ]
}
