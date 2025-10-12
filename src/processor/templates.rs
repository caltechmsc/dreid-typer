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
