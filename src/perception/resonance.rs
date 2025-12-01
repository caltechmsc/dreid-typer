//! Detects specific, strong resonance systems via strict substructure matching.
//!
//! Unlike generalized conjugation detection, this module uses an allowlist of
//! chemically significant motifs (Carboxylate, Nitro, Guanidinium, Amide).
//! When a motif is found, its atoms are marked `is_resonant`, and the system
//! (atoms + bonds) is recorded to ensure the correct bond order in the topology.

use super::model::{AnnotatedMolecule, ResonanceSystem};
use crate::core::error::PerceptionError;
use crate::core::properties::{Element, GraphBondOrder};

/// Runs strict resonance perception in two phases.
///
/// First, it identifies core resonance systems (aromatics are handled upstream) like
/// carboxylates and amides, registering both their atoms and bonds. Second, it
//  propagates the `is_resonant` flag to adjacent heteroatoms that can participate
//  in conjugation.
///
/// Note: Aromatic resonance is handled in the `aromaticity` module. This module
/// focuses on non-aromatic conjugated systems.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose atoms will be tagged and systems recorded.
///
/// # Returns
///
/// `Ok(())` always, as this process is infallible.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    detect_core_functional_groups(molecule);
    propagate_resonance_to_periphery(molecule);
    Ok(())
}

/// Detects core resonance systems via substructure matching.
fn detect_core_functional_groups(molecule: &mut AnnotatedMolecule) {
    let mut processed = vec![false; molecule.atoms.len()];

    detect_carboxylate_groups(molecule, &mut processed);
    detect_nitro_groups(molecule, &mut processed);
    detect_guanidinium_groups(molecule, &mut processed);
    detect_amide_groups(molecule, &mut processed);
}

/// Propagates resonance flags to peripheral heteroatoms bonded to resonant systems.
fn propagate_resonance_to_periphery(molecule: &mut AnnotatedMolecule) {
    let mut newly_resonant = Vec::new();

    for i in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[i];

        if atom.is_resonant
            || !matches!(atom.element, Element::O | Element::N | Element::S)
            || atom.lone_pairs == 0
        {
            continue;
        }

        let is_bonded_to_resonant_atom = molecule.adjacency[i]
            .iter()
            .any(|&(neighbor_id, _)| molecule.atoms[neighbor_id].is_resonant);

        if is_bonded_to_resonant_atom {
            newly_resonant.push(i);
        }
    }

    for atom_id in newly_resonant {
        molecule.atoms[atom_id].is_resonant = true;
    }
}

/// Helper to register a detected core system.
fn register_core_system(
    molecule: &mut AnnotatedMolecule,
    atoms: &[usize],
    bonds: &[usize],
    processed: &mut [bool],
) {
    for &atom_id in atoms {
        molecule.atoms[atom_id].is_resonant = true;
        processed[atom_id] = true;
    }

    let mut sys_atoms = atoms.to_vec();
    let mut sys_bonds = bonds.to_vec();
    sys_atoms.sort_unstable();
    sys_bonds.sort_unstable();

    molecule.resonance_systems.push(ResonanceSystem {
        atom_ids: sys_atoms,
        bond_ids: sys_bonds,
    });
}

/// Finds the bond ID connecting two atoms. Panics if not found (internal consistency check).
fn find_bond_id(molecule: &AnnotatedMolecule, u: usize, v: usize) -> usize {
    molecule
        .bonds
        .iter()
        .find(|b| {
            (b.atom_ids.0 == u && b.atom_ids.1 == v) || (b.atom_ids.0 == v && b.atom_ids.1 == u)
        })
        .map(|b| b.id)
        .expect("Bond must exist between adjacent atoms")
}

/// Detects Carboxylate groups: C(=O)O-
fn detect_carboxylate_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let mut double_o = None;
        let mut single_o = None;

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                match order {
                    GraphBondOrder::Double => double_o = Some(neighbor_id),
                    GraphBondOrder::Single if molecule.atoms[neighbor_id].degree == 1 => {
                        single_o = Some(neighbor_id);
                    }
                    _ => {}
                }
            }
        }

        if let (Some(o1), Some(o2)) = (double_o, single_o) {
            let b1 = find_bond_id(molecule, c_idx, o1);
            let b2 = find_bond_id(molecule, c_idx, o2);
            register_core_system(molecule, &[c_idx, o1, o2], &[b1, b2], processed);
        }
    }
}

/// Detects Nitro groups: N(=O)O-
fn detect_nitro_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx] || molecule.atoms[n_idx].element != Element::N {
            continue;
        }

        let mut oxygen_neighbors = Vec::new();
        for &(neighbor_id, _) in &molecule.adjacency[n_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                oxygen_neighbors.push(neighbor_id);
            }
        }

        if oxygen_neighbors.len() == 2 {
            let o1 = oxygen_neighbors[0];
            let o2 = oxygen_neighbors[1];
            let b1 = find_bond_id(molecule, n_idx, o1);
            let b2 = find_bond_id(molecule, n_idx, o2);
            register_core_system(molecule, &[n_idx, o1, o2], &[b1, b2], processed);
        }
    }
}

/// Detects Guanidinium groups: C(N)(N)N+
fn detect_guanidinium_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let n_neighbors: Vec<_> = molecule.adjacency[c_idx]
            .iter()
            .filter(|(id, _)| molecule.atoms[*id].element == Element::N)
            .map(|(id, _)| *id)
            .collect();

        if n_neighbors.len() == 3 {
            let bonds: Vec<_> = n_neighbors
                .iter()
                .map(|&n_id| find_bond_id(molecule, c_idx, n_id))
                .collect();
            let mut atoms = vec![c_idx];
            atoms.extend(n_neighbors);
            register_core_system(molecule, &atoms, &bonds, processed);
        }
    }
}

/// Detects Amide groups: O=C-N
fn detect_amide_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let mut double_o = None;
        let mut single_n = None;

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            match (molecule.atoms[neighbor_id].element, order) {
                (Element::O, GraphBondOrder::Double) => double_o = Some(neighbor_id),
                (Element::N, GraphBondOrder::Single) => single_n = Some(neighbor_id),
                _ => {}
            }
        }

        if let (Some(o), Some(n)) = (double_o, single_n) {
            let b_co = find_bond_id(molecule, c_idx, o);
            let b_cn = find_bond_id(molecule, c_idx, n);
            register_core_system(molecule, &[c_idx, o, n], &[b_co, b_cn], processed);
        }
    }
}
