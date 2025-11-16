use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::{BondOrder, Element};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let conjugated_systems =
        pauling::find_resonance_systems(molecule).map_err(PerceptionError::PaulingError)?;

    for system in conjugated_systems {
        for atom_id in system.atoms {
            if let Some(atom) = molecule.atoms.get_mut(atom_id) {
                atom.is_in_conjugated_system = true;
            } else {
                return Err(PerceptionError::Other(format!(
                    "pauling library returned an invalid atom ID ({}) that is out of bounds",
                    atom_id
                )));
            }
        }
    }

    apply_local_resonance_patterns(molecule);

    Ok(())
}

fn apply_local_resonance_patterns(molecule: &mut AnnotatedMolecule) {
    mark_aromatic_atoms_conjugated(molecule);
    mark_amide_and_thioamide_systems(molecule);
    mark_sulfonamide_systems(molecule);
    suppress_halogen_oxyanion_conjugation(molecule);
    demote_sigma_bound_sulfurs(molecule);
}

fn mark_aromatic_atoms_conjugated(molecule: &mut AnnotatedMolecule) {
    for atom in &mut molecule.atoms {
        if atom.is_aromatic {
            atom.is_in_conjugated_system = true;
        }
    }
}

fn mark_amide_and_thioamide_systems(molecule: &mut AnnotatedMolecule) {
    for pivot_idx in 0..molecule.atoms.len() {
        if molecule.atoms[pivot_idx].element != Element::C {
            continue;
        }

        let pi_partners: Vec<_> = molecule.adjacency[pivot_idx]
            .iter()
            .filter(|&&(_, order)| order == BondOrder::Double)
            .filter_map(|&(neighbor_id, _)| {
                let neighbor = &molecule.atoms[neighbor_id];
                matches!(neighbor.element, Element::O | Element::S).then_some(neighbor_id)
            })
            .collect();

        if pi_partners.is_empty() {
            continue;
        }

        let hetero_donors: Vec<_> = molecule.adjacency[pivot_idx]
            .iter()
            .filter(|&&(_, order)| order == BondOrder::Single)
            .filter_map(|&(neighbor_id, _)| {
                let neighbor = &molecule.atoms[neighbor_id];
                if neighbor_id != pivot_idx
                    && matches!(neighbor.element, Element::N | Element::O | Element::S)
                    && neighbor.lone_pairs > 0
                {
                    Some(neighbor_id)
                } else {
                    None
                }
            })
            .collect();

        if hetero_donors.is_empty() {
            continue;
        }

        mark_atoms_conjugated(molecule, [pivot_idx]);

        for pi_partner in pi_partners {
            mark_atoms_conjugated(molecule, [pi_partner]);
            for &donor in &hetero_donors {
                mark_atoms_conjugated(molecule, [donor]);
            }
        }
    }
}

fn mark_sulfonamide_systems(molecule: &mut AnnotatedMolecule) {
    for s_idx in 0..molecule.atoms.len() {
        if molecule.atoms[s_idx].element != Element::S {
            continue;
        }

        let double_bonded_oxygens: Vec<_> = molecule.adjacency[s_idx]
            .iter()
            .filter(|&&(_, order)| order == BondOrder::Double)
            .filter_map(|&(neighbor_id, _)| {
                (molecule.atoms[neighbor_id].element == Element::O).then_some(neighbor_id)
            })
            .collect();

        if double_bonded_oxygens.len() < 2 {
            continue;
        }

        let sulfonamide_neighbors: Vec<_> = molecule.adjacency[s_idx]
            .iter()
            .filter(|&&(_, order)| order == BondOrder::Single)
            .filter_map(|&(neighbor_id, _)| {
                let neighbor = &molecule.atoms[neighbor_id];
                (neighbor.element == Element::N && neighbor.lone_pairs > 0).then_some(neighbor_id)
            })
            .collect();

        for neighbor_id in sulfonamide_neighbors {
            mark_atoms_conjugated(molecule, [s_idx, neighbor_id]);
        }
    }
}

fn suppress_halogen_oxyanion_conjugation(molecule: &mut AnnotatedMolecule) {
    for center_idx in 0..molecule.atoms.len() {
        if !matches!(
            molecule.atoms[center_idx].element,
            Element::Cl | Element::Br | Element::I
        ) {
            continue;
        }

        let oxygen_neighbors: Vec<_> = molecule.adjacency[center_idx]
            .iter()
            .filter_map(|&(neighbor_id, _)| {
                (molecule.atoms[neighbor_id].element == Element::O).then_some(neighbor_id)
            })
            .collect();

        if oxygen_neighbors.len() >= 3 {
            for oxygen_idx in oxygen_neighbors {
                if let Some(atom) = molecule.atoms.get_mut(oxygen_idx) {
                    atom.is_in_conjugated_system = false;
                }
            }
        }
    }
}

fn demote_sigma_bound_sulfurs(molecule: &mut AnnotatedMolecule) {
    for s_idx in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[s_idx];
        if atom.element != Element::S || !atom.is_in_conjugated_system {
            continue;
        }

        let has_pi_bond = molecule.adjacency[s_idx]
            .iter()
            .any(|&(_, order)| order != BondOrder::Single);

        if !has_pi_bond {
            if let Some(atom_mut) = molecule.atoms.get_mut(s_idx) {
                atom_mut.is_in_conjugated_system = false;
            }
        }
    }
}

fn mark_atoms_conjugated<const N: usize>(molecule: &mut AnnotatedMolecule, atom_ids: [usize; N]) {
    for atom_id in atom_ids {
        if let Some(atom) = molecule.atoms.get_mut(atom_id) {
            atom.is_in_conjugated_system = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::properties::{BondOrder, Element};
    use std::collections::BTreeSet;

    fn build_molecule(
        elements: &[Element],
        bonds: &[(usize, usize, BondOrder)],
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

    fn hydrocarbon(
        backbone_bonds: &[(usize, usize, BondOrder)],
        hydrogen_counts: &[u8],
    ) -> AnnotatedMolecule {
        let heavy_atoms = hydrogen_counts.len();
        let mut elements = vec![Element::C; heavy_atoms];
        let mut bonds = backbone_bonds.to_vec();
        let mut next_index = heavy_atoms;
        for (atom_idx, &hydrogens) in hydrogen_counts.iter().enumerate() {
            for _ in 0..hydrogens {
                elements.push(Element::H);
                bonds.push((atom_idx, next_index, BondOrder::Single));
                next_index += 1;
            }
        }
        build_molecule(&elements, &bonds)
    }

    fn run_resonance(mut molecule: AnnotatedMolecule) -> AnnotatedMolecule {
        perceive(&mut molecule).expect("resonance perception should succeed");
        molecule
    }

    fn assert_conjugated_atoms(molecule: &AnnotatedMolecule, expected: &[usize]) {
        let observed: BTreeSet<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| atom.is_in_conjugated_system.then_some(idx))
            .collect();
        let anticipated: BTreeSet<_> = expected.iter().copied().collect();
        assert_eq!(
            observed, anticipated,
            "unexpected conjugated atom assignment"
        );
    }

    fn butadiene() -> AnnotatedMolecule {
        hydrocarbon(
            &[
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Double),
            ],
            &[2, 1, 1, 2],
        )
    }

    fn benzene_ring() -> AnnotatedMolecule {
        hydrocarbon(
            &[
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Double),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Double),
                (5, 0, BondOrder::Single),
            ],
            &[1, 1, 1, 1, 1, 1],
        )
    }

    fn allyl_anion() -> AnnotatedMolecule {
        let mut molecule = hydrocarbon(
            &[(0, 1, BondOrder::Double), (1, 2, BondOrder::Single)],
            &[2, 1, 2],
        );
        molecule.atoms[2].formal_charge = -1;
        molecule
    }

    fn dual_diene_with_sp3_break() -> AnnotatedMolecule {
        hydrocarbon(
            &[
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Double),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Single),
                (5, 6, BondOrder::Double),
                (6, 7, BondOrder::Single),
                (7, 8, BondOrder::Double),
            ],
            &[2, 1, 1, 2, 2, 1, 1, 2, 2],
        )
    }

    fn hexane() -> AnnotatedMolecule {
        hydrocarbon(
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Single),
            ],
            &[3, 2, 2, 2, 2, 3],
        )
    }

    #[test]
    fn linear_diene_marks_expected_chain_atoms() {
        let molecule = run_resonance(butadiene());
        assert_conjugated_atoms(&molecule, &[0, 1, 2, 3]);
    }

    #[test]
    fn benzene_ring_forms_single_resonance_system() {
        let molecule = run_resonance(benzene_ring());
        assert_conjugated_atoms(&molecule, &[0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn allyl_anion_includes_anionic_carbon_in_conjugation() {
        let molecule = run_resonance(allyl_anion());
        assert_conjugated_atoms(&molecule, &[0, 1, 2]);
    }

    #[test]
    fn saturated_breaks_split_disconnected_conjugated_systems() {
        let molecule = run_resonance(dual_diene_with_sp3_break());
        assert_conjugated_atoms(&molecule, &[0, 1, 2, 3, 5, 6, 7, 8]);
    }

    #[test]
    fn saturated_hexane_has_no_conjugation() {
        let molecule = run_resonance(hexane());
        assert_conjugated_atoms(&molecule, &[]);
    }
}
