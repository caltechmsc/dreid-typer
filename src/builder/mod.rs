use crate::core::Hybridization;
use crate::core::error::TyperError;
use crate::core::graph::{
    Angle, Atom, Bond, ImproperDihedral, MolecularGraph, MolecularTopology, ProperDihedral,
};
use crate::processor::ProcessingGraph;
use std::collections::HashSet;

pub(crate) fn build_topology(
    initial_graph: &MolecularGraph,
    processing_graph: &ProcessingGraph,
    atom_types: &[String],
) -> Result<MolecularTopology, TyperError> {
    let atoms = build_atoms(processing_graph, atom_types);
    let bonds = build_bonds(initial_graph);
    let angles = build_angles(processing_graph);
    let proper_dihedrals = build_proper_dihedrals(processing_graph);
    let improper_dihedrals = build_improper_dihedrals(processing_graph);

    Ok(MolecularTopology {
        atoms,
        bonds: bonds.into_iter().collect(),
        angles: angles.into_iter().collect(),
        proper_dihedrals: proper_dihedrals.into_iter().collect(),
        improper_dihedrals: improper_dihedrals.into_iter().collect(),
    })
}

fn build_atoms(processing_graph: &ProcessingGraph, atom_types: &[String]) -> Vec<Atom> {
    processing_graph
        .atoms
        .iter()
        .map(|atom_view| Atom {
            id: atom_view.id,
            element: atom_view.element,
            atom_type: atom_types[atom_view.id].clone(),
            hybridization: atom_view.hybridization,
        })
        .collect()
}

fn build_bonds(initial_graph: &MolecularGraph) -> HashSet<Bond> {
    initial_graph
        .bonds
        .iter()
        .map(|edge| {
            let (u, v) = edge.atom_ids;
            Bond::new(u, v, edge.order)
        })
        .collect()
}

fn build_angles(graph: &ProcessingGraph) -> HashSet<Angle> {
    let mut angles = HashSet::new();
    for j in 0..graph.atoms.len() {
        let neighbors = &graph.adjacency[j];
        if neighbors.len() < 2 {
            continue;
        }

        for i_idx in 0..neighbors.len() {
            for k_idx in (i_idx + 1)..neighbors.len() {
                let i = neighbors[i_idx].0;
                let k = neighbors[k_idx].0;
                angles.insert(Angle::new(i, j, k));
            }
        }
    }
    angles
}

fn build_proper_dihedrals(graph: &ProcessingGraph) -> HashSet<ProperDihedral> {
    let mut dihedrals = HashSet::new();
    for j in 0..graph.atoms.len() {
        for &(k, _) in &graph.adjacency[j] {
            if j >= k {
                continue;
            }

            for &(i, _) in &graph.adjacency[j] {
                if i == k {
                    continue;
                }

                for &(l, _) in &graph.adjacency[k] {
                    if l == j {
                        continue;
                    }

                    dihedrals.insert(ProperDihedral::new(i, j, k, l));
                }
            }
        }
    }
    dihedrals
}

fn build_improper_dihedrals(graph: &ProcessingGraph) -> HashSet<ImproperDihedral> {
    let mut dihedrals = HashSet::new();
    for i in 0..graph.atoms.len() {
        let atom = &graph.atoms[i];

        if atom.degree == 3 {
            if matches!(
                atom.hybridization,
                Hybridization::SP2 | Hybridization::Resonant
            ) {
                let neighbors = &graph.adjacency[i];
                let j = neighbors[0].0;
                let k = neighbors[1].0;
                let l = neighbors[2].0;

                dihedrals.insert(ImproperDihedral::new(j, k, i, l));
            }
        }
    }
    dihedrals
}
