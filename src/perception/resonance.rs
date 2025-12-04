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
    detect_thiourea_groups(molecule, &mut processed);
    detect_amide_groups(molecule, &mut processed);
    detect_phosphate_groups(molecule, &mut processed);
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

/// Helper to record a detected resonance system for later topology emission.
fn push_resonance_system(molecule: &mut AnnotatedMolecule, atoms: &[usize], bonds: &[usize]) {
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
            for &atom_id in &[c_idx, o1, o2] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &[c_idx, o1, o2], &[b1, b2]);
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
            for &atom_id in &[n_idx, o1, o2] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &[n_idx, o1, o2], &[b1, b2]);
        }
    }
}

/// Detects Guanidinium groups: C(N)(N)N+
fn detect_guanidinium_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx]
            || molecule.atoms[c_idx].element != Element::C
            || molecule.atoms[c_idx].is_in_ring
        {
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
            for &atom_id in &atoms {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &atoms, &bonds);
        }
    }
}

/// Detects Thiourea cores: C(=S)(N)(N)
fn detect_thiourea_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let mut sulfur_neighbor = None;
        let mut nitrogen_neighbors = Vec::new();

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            match (molecule.atoms[neighbor_id].element, order) {
                (Element::S, GraphBondOrder::Double) => sulfur_neighbor = Some(neighbor_id),
                (Element::N, GraphBondOrder::Single) => nitrogen_neighbors.push(neighbor_id),
                _ => {}
            }
        }

        if sulfur_neighbor.is_some() && nitrogen_neighbors.len() == 2 {
            let n1 = nitrogen_neighbors[0];
            let n2 = nitrogen_neighbors[1];
            let b1 = find_bond_id(molecule, c_idx, n1);
            let b2 = find_bond_id(molecule, c_idx, n2);

            for &atom_id in &[c_idx, n1, n2] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }

            push_resonance_system(molecule, &[c_idx, n1, n2], &[b1, b2]);
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
            if molecule.atoms[c_idx].is_resonant || molecule.atoms[n].is_resonant {
                continue;
            }
            let b_co = find_bond_id(molecule, c_idx, o);
            let b_cn = find_bond_id(molecule, c_idx, n);
            for &atom_id in &[c_idx, o, n] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &[c_idx, o, n], &[b_co, b_cn]);
        }
    }
}

/// Detects Phosphate groups: P(=O)(O-)
fn detect_phosphate_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for p_idx in 0..molecule.atoms.len() {
        if molecule.atoms[p_idx].element != Element::P {
            continue;
        }

        let mut terminal_oxygens = Vec::new();

        for &(neighbor_id, _) in &molecule.adjacency[p_idx] {
            if molecule.atoms[neighbor_id].element == Element::O
                && molecule.atoms[neighbor_id].degree == 1
            {
                terminal_oxygens.push(neighbor_id);
            }
        }

        if terminal_oxygens.len() >= 2 {
            let o1 = terminal_oxygens[0];
            let o2 = terminal_oxygens[1];

            let b1 = find_bond_id(molecule, p_idx, o1);
            let b2 = find_bond_id(molecule, p_idx, o2);

            molecule.atoms[o1].is_resonant = true;
            molecule.atoms[o2].is_resonant = true;
            processed[o1] = true;
            processed[o2] = true;
            processed[p_idx] = true;

            push_resonance_system(molecule, &[p_idx, o1, o2], &[b1, b2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::properties::Element;

    fn build_molecule(
        elements: &[Element],
        bonds: &[(usize, usize, GraphBondOrder)],
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for &element in elements {
            graph.add_atom(element);
        }
        for &(u, v, order) in bonds {
            graph.add_bond(u, v, order).expect("bond endpoints exist");
        }
        AnnotatedMolecule::new(&graph).expect("graph must be chemically valid")
    }

    #[test]
    fn carboxylate_is_detected_as_core_system() {
        let mut molecule = build_molecule(
            &[Element::C, Element::O, Element::O],
            &[
                (0, 1, GraphBondOrder::Double),
                (0, 2, GraphBondOrder::Single),
            ],
        );
        molecule.atoms[0].degree = 2;
        molecule.atoms[1].degree = 1;
        molecule.atoms[2].degree = 1;

        perceive(&mut molecule).unwrap();

        assert!(
            molecule.atoms[0].is_resonant,
            "Carboxylate C should be resonant"
        );
        assert!(
            molecule.atoms[1].is_resonant,
            "Carboxylate O= should be resonant"
        );
        assert!(
            molecule.atoms[2].is_resonant,
            "Carboxylate O- should be resonant"
        );

        assert_eq!(molecule.resonance_systems.len(), 1);
        let system = &molecule.resonance_systems[0];
        assert_eq!(system.atom_ids, vec![0, 1, 2]);
        assert_eq!(system.bond_ids, vec![0, 1]);
    }

    #[test]
    fn nitro_is_detected_regardless_of_kekule_form() {
        let molecule1 = {
            let mut mol = build_molecule(
                &[Element::N, Element::O, Element::O],
                &[
                    (0, 1, GraphBondOrder::Double),
                    (0, 2, GraphBondOrder::Single),
                ],
            );
            perceive(&mut mol).unwrap();
            mol
        };

        let molecule2 = {
            let mut mol = build_molecule(
                &[Element::N, Element::O, Element::O],
                &[
                    (0, 1, GraphBondOrder::Single),
                    (0, 2, GraphBondOrder::Double),
                ],
            );
            perceive(&mut mol).unwrap();
            mol
        };

        for (i, molecule) in [molecule1, molecule2].iter().enumerate() {
            assert!(
                molecule.atoms[0].is_resonant,
                "Nitro N should be resonant (form {})",
                i + 1
            );
            assert!(
                molecule.atoms[1].is_resonant,
                "Nitro O should be resonant (form {})",
                i + 1
            );
            assert!(
                molecule.atoms[2].is_resonant,
                "Nitro O should be resonant (form {})",
                i + 1
            );

            assert_eq!(
                molecule.resonance_systems.len(),
                1,
                "Should detect 1 nitro system (form {})",
                i + 1
            );
            let system = &molecule.resonance_systems[0];
            assert_eq!(system.atom_ids, vec![0, 1, 2]);
            assert_eq!(system.bond_ids, vec![0, 1]);
        }
    }

    #[test]
    fn amide_is_detected_as_core_system() {
        let molecule = {
            let mut mol = build_molecule(
                &[Element::C, Element::O, Element::N],
                &[
                    (0, 1, GraphBondOrder::Double),
                    (0, 2, GraphBondOrder::Single),
                ],
            );
            perceive(&mut mol).unwrap();
            mol
        };

        assert!(molecule.atoms[0].is_resonant, "Amide C should be resonant");
        assert!(molecule.atoms[1].is_resonant, "Amide O should be resonant");
        assert!(molecule.atoms[2].is_resonant, "Amide N should be resonant");

        assert_eq!(molecule.resonance_systems.len(), 1);
        let system = &molecule.resonance_systems[0];
        assert_eq!(system.atom_ids, vec![0, 1, 2]);
        assert_eq!(system.bond_ids, vec![0, 1]);
    }

    #[test]
    fn guanidinium_is_detected_as_core_system() {
        let molecule = {
            let mut mol = build_molecule(
                &[Element::C, Element::N, Element::N, Element::N],
                &[
                    (0, 1, GraphBondOrder::Double),
                    (0, 2, GraphBondOrder::Single),
                    (0, 3, GraphBondOrder::Single),
                ],
            );
            perceive(&mut mol).unwrap();
            mol
        };

        assert!((0..=3).all(|i| molecule.atoms[i].is_resonant));
        assert_eq!(molecule.resonance_systems.len(), 1);
        let system = &molecule.resonance_systems[0];
        assert_eq!(system.atom_ids, vec![0, 1, 2, 3]);
        assert_eq!(system.bond_ids, vec![0, 1, 2]);
    }

    #[test]
    fn phenol_oxygen_is_marked_resonant_via_propagation() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C, Element::O, Element::H],
            &[
                (0, 1, GraphBondOrder::Single),
                (1, 2, GraphBondOrder::Single),
                (2, 3, GraphBondOrder::Single),
            ],
        );
        molecule.atoms[1].is_resonant = true;
        molecule.atoms[2].lone_pairs = 2;

        perceive(&mut molecule).unwrap();

        assert!(
            molecule.atoms[2].is_resonant,
            "Phenol oxygen should be marked resonant"
        );
        assert!(
            molecule.resonance_systems.is_empty(),
            "Peripheral resonance should not create a new core system"
        );
    }

    #[test]
    fn non_resonant_structures_are_ignored() {
        let molecule = {
            let mut mol = build_molecule(
                &[Element::C, Element::C, Element::O],
                &[
                    (0, 1, GraphBondOrder::Single),
                    (1, 2, GraphBondOrder::Single),
                ],
            );
            perceive(&mut mol).unwrap();
            mol
        };

        assert!(!molecule.atoms[0].is_resonant);
        assert!(!molecule.atoms[1].is_resonant);
        assert!(!molecule.atoms[2].is_resonant);
        assert!(molecule.resonance_systems.is_empty());
    }
}
