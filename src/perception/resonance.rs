use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;

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

    Ok(())
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
