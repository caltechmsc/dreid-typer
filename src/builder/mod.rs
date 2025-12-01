//! Converts annotated molecules and assigned atom types into a full molecular topology.
//!
//! The builder stage takes the perception output and typing assignments, emitting atoms, bonds,
//! angles, proper dihedrals, and improper dihedrals expected by downstream force-field tooling.

use crate::core::properties::{GraphBondOrder, Hybridization, TopologyBondOrder};
use crate::core::topology::{
    Angle, Atom, Bond, ImproperDihedral, MolecularTopology, ProperDihedral,
};
use crate::perception::AnnotatedMolecule;
use std::collections::HashSet;

/// Builds the `MolecularTopology` aggregate from perception results and atom-type labels.
///
/// This function effectively serializes the `AnnotatedMolecule` into the graph structures used
/// by force-field consumers by delegating to specialized helpers for each topology term.
///
/// # Arguments
///
/// * `annotated_molecule` - Molecule carrying ring, hybridization, and bonding metadata.
/// * `atom_types` - Slice of final atom-type names aligned with the molecule's atom ordering.
///
/// # Returns
///
/// A populated [`MolecularTopology`] containing atoms, bonds, angles, proper dihedrals, and
/// improper dihedrals.
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

/// Creates the atom list with element, type, and hybridization copies.
///
/// # Arguments
///
/// * `annotated_molecule` - Source molecule whose atoms provide structural metadata.
/// * `atom_types` - Slice of assigned atom-type labels.
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

/// Extracts unique bonds from the annotated molecule.
///
/// This function determines the final `TopologyBondOrder` by checking if a bond belongs to
/// any detected `ResonanceSystem`. If so, it is promoted to `Resonant`. Otherwise, the
/// Kekulized order (`Single`, `Double`, `Triple`) is used.
fn build_bonds(annotated_molecule: &AnnotatedMolecule) -> HashSet<Bond> {
    let resonant_bond_ids: HashSet<usize> = annotated_molecule
        .resonance_systems
        .iter()
        .flat_map(|sys| sys.bond_ids.iter())
        .copied()
        .collect();

    annotated_molecule
        .bonds
        .iter()
        .map(|edge| {
            let topology_order = if resonant_bond_ids.contains(&edge.id) {
                TopologyBondOrder::Resonant
            } else {
                match edge.order {
                    GraphBondOrder::Single => TopologyBondOrder::Single,
                    GraphBondOrder::Double => TopologyBondOrder::Double,
                    GraphBondOrder::Triple => TopologyBondOrder::Triple,
                    GraphBondOrder::Aromatic => TopologyBondOrder::Single, // Fallback; should not occur here
                }
            };

            Bond::new(edge.atom_ids.0, edge.atom_ids.1, topology_order)
        })
        .collect()
}

/// Generates all angle triplets by enumerating neighbor pairs around each atom.
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

/// Builds proper dihedrals by extending each bond to its neighboring atoms.
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

/// Builds improper dihedrals for planar degree-three centers with SP2-like hybridization.
fn build_impropers(annotated_molecule: &AnnotatedMolecule) -> HashSet<ImproperDihedral> {
    let mut impropers = HashSet::new();
    for atom in &annotated_molecule.atoms {
        if atom.degree == 3
            && matches!(
                atom.hybridization,
                Hybridization::SP2 | Hybridization::Resonant
            )
        {
            let neighbors = &annotated_molecule.adjacency[atom.id];
            let p1 = neighbors[0].0;
            let p2 = neighbors[1].0;
            let p3 = neighbors[2].0;
            impropers.insert(ImproperDihedral::new(p1, p2, atom.id, p3));
        }
    }
    impropers
}
