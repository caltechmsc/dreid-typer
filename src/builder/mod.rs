use crate::core::graph::{Angle, Atom, Bond, ImproperDihedral, MolecularTopology, ProperDihedral};
use crate::core::properties::Hybridization;
use crate::perception::AnnotatedMolecule;
use std::collections::HashSet;

pub fn build_topology(
    annotated_molecule: &AnnotatedMolecule,
    atom_types: &[String],
) -> MolecularTopology {
    let atoms = build_atoms(annotated_molecule, atom_types);
    let bonds = build_bonds(annotated_molecule);
    let angles = build_angles(annotated_molecule);
    let propers = build_propers(annotated_molecule);
    let impropers = build_impropers(annotated_molecule);

    MolecularTopology {
        atoms,
        bonds: bonds.into_iter().collect(),
        angles: angles.into_iter().collect(),
        propers: propers.into_iter().collect(),
        impropers: impropers.into_iter().collect(),
    }
}

fn build_atoms(annotated_molecule: &AnnotatedMolecule, atom_types: &[String]) -> Vec<Atom> {
    annotated_molecule
        .atoms
        .iter()
        .map(|ann_atom| Atom {
            id: ann_atom.id,
            element: ann_atom.element,
            atom_type: atom_types[ann_atom.id].clone(),
            hybridization: ann_atom.hybridization,
        })
        .collect()
}

fn build_bonds(annotated_molecule: &AnnotatedMolecule) -> HashSet<Bond> {
    annotated_molecule
        .bonds
        .iter()
        .map(|edge| Bond::new(edge.atom_ids.0, edge.atom_ids.1, edge.order))
        .collect()
}

fn build_angles(annotated_molecule: &AnnotatedMolecule) -> HashSet<Angle> {
    let mut angles = HashSet::new();
    for j in 0..annotated_molecule.atoms.len() {
        let neighbors = &annotated_molecule.adjacency[j];
        if neighbors.len() < 2 {
            continue;
        }
        for i in 0..neighbors.len() {
            for k in (i + 1)..neighbors.len() {
                let atom_i_id = neighbors[i].0;
                let atom_k_id = neighbors[k].0;
                angles.insert(Angle::new(atom_i_id, j, atom_k_id));
            }
        }
    }
    angles
}

fn build_propers(annotated_molecule: &AnnotatedMolecule) -> HashSet<ProperDihedral> {
    let mut propers = HashSet::new();
    for bond_jk in &annotated_molecule.bonds {
        let (j, k) = bond_jk.atom_ids;

        for &(i, _) in &annotated_molecule.adjacency[j] {
            if i == k {
                continue;
            }
            for &(l, _) in &annotated_molecule.adjacency[k] {
                if l == j || l == i {
                    continue;
                }
                propers.insert(ProperDihedral::new(i, j, k, l));
            }
        }
    }
    propers
}

fn build_impropers(annotated_molecule: &AnnotatedMolecule) -> HashSet<ImproperDihedral> {
    let mut impropers = HashSet::new();
    for atom in &annotated_molecule.atoms {
        if atom.degree == 3 {
            if matches!(
                atom.hybridization,
                Hybridization::SP2 | Hybridization::Resonant
            ) {
                let neighbors = &annotated_molecule.adjacency[atom.id];
                let p1 = neighbors[0].0;
                let p2 = neighbors[1].0;
                let p3 = neighbors[2].0;
                impropers.insert(ImproperDihedral::new(p1, p2, atom.id, p3));
            }
        }
    }
    impropers
}
