use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::{BondOrder, Element};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let mut processed = vec![false; molecule.atoms.len()];

    assign_nitrone_groups(molecule, &mut processed)?;
    assign_nitro_groups(molecule, &mut processed)?;
    assign_sulfur_oxides(molecule, &mut processed)?;
    assign_phosphorus_oxides(molecule, &mut processed)?;
    assign_carboxylate_anions(molecule, &mut processed)?;
    assign_ammonium_and_iminium(molecule, &mut processed)?;
    assign_onium_ions(molecule, &mut processed)?;
    assign_phosphonium_ions(molecule, &mut processed)?;
    assign_enolate_phenate_anions(molecule, &mut processed)?;

    assign_general(molecule, &processed)?;

    Ok(())
}

fn assign_nitrone_groups(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx]
            || molecule.atoms[n_idx].element != Element::N
            || molecule.atoms[n_idx].degree != 3
        {
            continue;
        }

        let mut double_bond_c_idx = None;
        let mut single_bond_o_idx = None;
        let mut single_bond_c_idx = None;

        for &(neighbor_id, order) in &molecule.adjacency[n_idx] {
            match (molecule.atoms[neighbor_id].element, order) {
                (Element::C, BondOrder::Double) => double_bond_c_idx = Some(neighbor_id),
                (Element::O, BondOrder::Single) => single_bond_o_idx = Some(neighbor_id),
                (Element::C, BondOrder::Single) => single_bond_c_idx = Some(neighbor_id),
                _ => {}
            }
        }

        match (double_bond_c_idx, single_bond_o_idx, single_bond_c_idx) {
            (Some(c1), Some(o), Some(c2)) if !processed[c1] && !processed[o] && !processed[c2] => {
                molecule.atoms[n_idx].formal_charge = 1;
                molecule.atoms[n_idx].lone_pairs = 0;
                molecule.atoms[o].formal_charge = -1;
                molecule.atoms[o].lone_pairs = 3;

                processed[n_idx] = true;
                processed[o] = true;
            }
            _ => {}
        }
    }
    Ok(())
}

fn assign_nitro_groups(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx]
            || molecule.atoms[n_idx].element != Element::N
            || molecule.atoms[n_idx].degree != 3
        {
            continue;
        }

        let mut double_bond_o_idx = None;
        let mut single_bond_o_idx = None;
        let mut other_neighbor_count = 0;

        for &(neighbor_id, order) in &molecule.adjacency[n_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                match order {
                    BondOrder::Double if double_bond_o_idx.is_none() => {
                        double_bond_o_idx = Some(neighbor_id)
                    }
                    BondOrder::Single if single_bond_o_idx.is_none() => {
                        single_bond_o_idx = Some(neighbor_id)
                    }
                    _ => continue,
                }
            } else {
                other_neighbor_count += 1;
            }
        }

        match (double_bond_o_idx, single_bond_o_idx) {
            (Some(o1), Some(o2))
                if other_neighbor_count == 1 && !processed[o1] && !processed[o2] =>
            {
                molecule.atoms[n_idx].formal_charge = 1;
                molecule.atoms[n_idx].lone_pairs = 0;

                molecule.atoms[o1].formal_charge = 0;
                molecule.atoms[o1].lone_pairs = 2;

                molecule.atoms[o2].formal_charge = -1;
                molecule.atoms[o2].lone_pairs = 3;

                processed[n_idx] = true;
                processed[o1] = true;
                processed[o2] = true;
            }
            _ => {}
        }
    }
    Ok(())
}

fn assign_sulfur_oxides(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for s_idx in 0..molecule.atoms.len() {
        if processed[s_idx] || molecule.atoms[s_idx].element != Element::S {
            continue;
        }

        let double_bonded_oxygens: Vec<usize> = molecule.adjacency[s_idx]
            .iter()
            .filter(|&&(id, order)| {
                molecule.atoms[id].element == Element::O && order == BondOrder::Double
            })
            .map(|&(id, _)| id)
            .collect();

        match (molecule.atoms[s_idx].degree, double_bonded_oxygens.len()) {
            (3, 1) => {
                let o_idx = double_bonded_oxygens[0];
                if !processed[o_idx] {
                    molecule.atoms[s_idx].formal_charge = 1;
                    molecule.atoms[s_idx].lone_pairs = 1;
                    molecule.atoms[o_idx].formal_charge = -1;
                    molecule.atoms[o_idx].lone_pairs = 3;
                    processed[s_idx] = true;
                    processed[o_idx] = true;
                }
            }
            (4, 2) => {
                let o1_idx = double_bonded_oxygens[0];
                let o2_idx = double_bonded_oxygens[1];
                if !processed[o1_idx] && !processed[o2_idx] {
                    molecule.atoms[s_idx].formal_charge = 2;
                    molecule.atoms[s_idx].lone_pairs = 0;
                    molecule.atoms[o1_idx].formal_charge = -1;
                    molecule.atoms[o1_idx].lone_pairs = 3;
                    molecule.atoms[o2_idx].formal_charge = -1;
                    molecule.atoms[o2_idx].lone_pairs = 3;
                    processed[s_idx] = true;
                    processed[o1_idx] = true;
                    processed[o2_idx] = true;
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn assign_phosphorus_oxides(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for p_idx in 0..molecule.atoms.len() {
        if processed[p_idx]
            || molecule.atoms[p_idx].element != Element::P
            || molecule.atoms[p_idx].degree != 4
        {
            continue;
        }

        let double_bonded_oxygens: Vec<usize> = molecule.adjacency[p_idx]
            .iter()
            .filter(|&&(id, order)| {
                molecule.atoms[id].element == Element::O && order == BondOrder::Double
            })
            .map(|&(id, _)| id)
            .collect();

        if double_bonded_oxygens.len() == 1 {
            let o_idx = double_bonded_oxygens[0];
            if !processed[o_idx] {
                molecule.atoms[p_idx].formal_charge = 1;
                molecule.atoms[p_idx].lone_pairs = 0;
                molecule.atoms[o_idx].formal_charge = -1;
                molecule.atoms[o_idx].lone_pairs = 3;
                processed[p_idx] = true;
                processed[o_idx] = true;
            }
        }
    }
    Ok(())
}

fn assign_carboxylate_anions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx]
            || molecule.atoms[c_idx].element != Element::C
            || molecule.atoms[c_idx].degree != 3
        {
            continue;
        }

        let mut double_bond_o_idx = None;
        let mut single_bond_o_idx = None;

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                match order {
                    BondOrder::Double if double_bond_o_idx.is_none() => {
                        double_bond_o_idx = Some(neighbor_id)
                    }
                    BondOrder::Single if single_bond_o_idx.is_none() => {
                        single_bond_o_idx = Some(neighbor_id)
                    }
                    _ => continue,
                }
            }
        }

        match (double_bond_o_idx, single_bond_o_idx) {
            (Some(o1), Some(o2)) if !processed[o1] && !processed[o2] => {
                molecule.atoms[c_idx].formal_charge = 0;
                molecule.atoms[c_idx].lone_pairs = 0;

                molecule.atoms[o1].formal_charge = 0;
                molecule.atoms[o1].lone_pairs = 2;

                molecule.atoms[o2].formal_charge = -1;
                molecule.atoms[o2].lone_pairs = 3;

                processed[c_idx] = true;
                processed[o1] = true;
                processed[o2] = true;
            }
            _ => {}
        }
    }
    Ok(())
}

fn assign_ammonium_and_iminium(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (n_idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag || molecule.atoms[n_idx].element != Element::N {
            continue;
        }

        let degree = molecule.atoms[n_idx].degree;
        let has_double_bond = molecule.adjacency[n_idx]
            .iter()
            .any(|&(_, order)| order == BondOrder::Double);

        if degree == 4 || (degree == 3 && has_double_bond) {
            let atom = &mut molecule.atoms[n_idx];
            atom.formal_charge = 1;
            atom.lone_pairs = 0;
            *processed_flag = true;
        }
    }
    Ok(())
}

fn assign_onium_ions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag {
            continue;
        }

        let element = molecule.atoms[idx].element;
        let degree = molecule.atoms[idx].degree;
        if (element == Element::O || element == Element::S) && degree == 3 {
            let atom_mut = &mut molecule.atoms[idx];
            atom_mut.formal_charge = 1;
            atom_mut.lone_pairs = 1;
            *processed_flag = true;
        }
    }
    Ok(())
}

fn assign_phosphonium_ions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (p_idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag
            || molecule.atoms[p_idx].element != Element::P
            || molecule.atoms[p_idx].degree != 4
        {
            continue;
        }

        let has_double_bond_o = molecule.adjacency[p_idx].iter().any(|&(id, order)| {
            molecule.atoms[id].element == Element::O && order == BondOrder::Double
        });

        if !has_double_bond_o {
            let atom = &mut molecule.atoms[p_idx];
            atom.formal_charge = 1;
            atom.lone_pairs = 0;
            *processed_flag = true;
        }
    }
    Ok(())
}

fn assign_enolate_phenate_anions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (o_idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag
            || molecule.atoms[o_idx].element != Element::O
            || molecule.atoms[o_idx].degree != 1
        {
            continue;
        }

        let (neighbor_id, order) = molecule.adjacency[o_idx][0];

        if order != BondOrder::Single {
            continue;
        }

        let neighbor = &molecule.atoms[neighbor_id];
        if neighbor.element == Element::C {
            let neighbor_is_sp2 = molecule.adjacency[neighbor_id]
                .iter()
                .any(|&(_, order)| order == BondOrder::Double);

            if neighbor_is_sp2 {
                let atom = &mut molecule.atoms[o_idx];
                atom.formal_charge = -1;
                atom.lone_pairs = 3;
                *processed_flag = true;
            }
        }
    }
    Ok(())
}

fn assign_general(
    molecule: &mut AnnotatedMolecule,
    processed: &[bool],
) -> Result<(), PerceptionError> {
    for (i, processed_flag) in processed.iter().enumerate() {
        if *processed_flag {
            continue;
        }
        let element = molecule.atoms[i].element;

        let valence = element.valence_electrons().ok_or_else(|| {
            PerceptionError::Other(format!(
                "valence electrons not defined for element {:?}",
                element
            ))
        })?;

        let bonding_electrons: u8 = molecule.adjacency[i]
            .iter()
            .map(|&(_, order)| bond_order_to_valence(order))
            .sum();

        let mut lone_pairs = 0;

        let is_second_period = matches!(
            element,
            Element::B | Element::C | Element::N | Element::O | Element::F
        );

        if element == Element::H {
            let bonded_electrons = bonding_electrons.saturating_mul(2);
            if bonded_electrons <= 2 {
                lone_pairs = (2 - bonded_electrons) / 2;
            }
        } else if is_second_period {
            let bonded_electrons = bonding_electrons.saturating_mul(2);
            if bonded_electrons <= 8 {
                lone_pairs = (8 - bonded_electrons) / 2;
            }
        } else if valence >= bonding_electrons {
            lone_pairs = (valence - bonding_electrons) / 2;
        }

        let formal_charge = valence as i8 - bonding_electrons as i8 - (lone_pairs * 2) as i8;

        let atom_mut = &mut molecule.atoms[i];
        atom_mut.lone_pairs = lone_pairs;
        atom_mut.formal_charge = formal_charge;
    }
    Ok(())
}

fn bond_order_to_valence(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => panic!("Aromatic bonds should have been kekulized by now."),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::properties::Element;

    fn build_molecule(
        elements: &[Element],
        bonds: &[(usize, usize, BondOrder)],
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for &element in elements {
            graph.add_atom(element);
        }
        for &(u, v, order) in bonds {
            graph.add_bond(u, v, order).expect("valid bond");
        }
        AnnotatedMolecule::new(&graph).expect("graph construction should succeed")
    }

    fn run_perception(
        elements: &[Element],
        bonds: &[(usize, usize, BondOrder)],
    ) -> AnnotatedMolecule {
        let mut molecule = build_molecule(elements, bonds);
        perceive(&mut molecule).expect("perception should succeed");
        molecule
    }

    fn assert_atom_state(molecule: &AnnotatedMolecule, idx: usize, charge: i8, lone_pairs: u8) {
        let atom = &molecule.atoms[idx];
        assert_eq!(
            atom.formal_charge, charge,
            "unexpected charge on atom {} ({:?})",
            idx, atom.element
        );
        assert_eq!(
            atom.lone_pairs, lone_pairs,
            "unexpected lone pair count on atom {} ({:?})",
            idx, atom.element
        );
    }

    #[test]
    fn nitrone_groups_receive_expected_charges() {
        let elements = vec![
            Element::C,
            Element::N,
            Element::O,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (0, 5, BondOrder::Single),
            (3, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
            (3, 8, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 1, 1, 0);
        assert_atom_state(&molecule, 2, -1, 3);
    }

    #[test]
    fn nitro_groups_assign_expected_formal_charges() {
        let elements = vec![
            Element::C,
            Element::N,
            Element::O,
            Element::O,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (1, 2, BondOrder::Double),
            (1, 3, BondOrder::Single),
            (1, 0, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (0, 5, BondOrder::Single),
            (0, 6, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 1, 1, 0);
        assert_atom_state(&molecule, 2, 0, 2);
        assert_atom_state(&molecule, 3, -1, 3);
    }

    #[test]
    fn sulfoxide_pattern_sets_expected_charges() {
        let elements = vec![
            Element::S,
            Element::O,
            Element::C,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (2, 4, BondOrder::Single),
            (2, 5, BondOrder::Single),
            (2, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
            (3, 8, BondOrder::Single),
            (3, 9, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 0, 1, 1);
        assert_atom_state(&molecule, 1, -1, 3);
    }

    #[test]
    fn sulfone_pattern_assigns_double_anionic_oxygens() {
        let elements = vec![
            Element::S,
            Element::O,
            Element::O,
            Element::C,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Double),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (3, 5, BondOrder::Single),
            (3, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
            (4, 8, BondOrder::Single),
            (4, 9, BondOrder::Single),
            (4, 10, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 0, 2, 0);
        assert_atom_state(&molecule, 1, -1, 3);
        assert_atom_state(&molecule, 2, -1, 3);
    }

    #[test]
    fn phosphorus_oxide_assigns_positive_phosphorus() {
        let elements = vec![Element::P, Element::O, Element::H, Element::H, Element::H];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 0, 1, 0);
        assert_atom_state(&molecule, 1, -1, 3);
    }

    #[test]
    fn carboxylate_anion_marks_single_bonded_oxygen() {
        let elements = vec![
            Element::C,
            Element::O,
            Element::O,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (3, 5, BondOrder::Single),
            (3, 6, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 1, 0, 2);
        assert_atom_state(&molecule, 2, -1, 3);
        assert_atom_state(&molecule, 0, 0, 0);
    }

    #[test]
    fn ammonium_and_iminium_patterns_assign_positive_nitrogen() {
        let ammonium = run_perception(
            &[Element::N, Element::H, Element::H, Element::H, Element::H],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
            ],
        );
        assert_atom_state(&ammonium, 0, 1, 0);

        let elements = vec![
            Element::C,
            Element::N,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (0, 5, BondOrder::Single),
            (2, 6, BondOrder::Single),
            (2, 7, BondOrder::Single),
            (2, 8, BondOrder::Single),
        ];
        let iminium = run_perception(&elements, &bonds);
        assert_atom_state(&iminium, 1, 1, 0);
    }

    #[test]
    fn onium_and_phosphonium_patterns_assign_positive_charges() {
        let oxonium = run_perception(
            &[Element::O, Element::H, Element::H, Element::H],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        assert_atom_state(&oxonium, 0, 1, 1);

        let phosphonium = run_perception(
            &[Element::P, Element::H, Element::H, Element::H, Element::H],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
            ],
        );
        assert_atom_state(&phosphonium, 0, 1, 0);
    }

    #[test]
    fn enolate_detection_marks_anionic_oxygen() {
        let elements = vec![
            Element::O,
            Element::C,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (1, 5, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (2, 4, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 0, -1, 3);
    }

    #[test]
    fn general_rules_handle_small_neutral_molecules() {
        let water = run_perception(
            &[Element::O, Element::H, Element::H],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        assert_atom_state(&water, 0, 0, 2);
        assert_atom_state(&water, 1, 0, 0);
        assert_atom_state(&water, 2, 0, 0);

        let methane = run_perception(
            &[Element::C, Element::H, Element::H, Element::H, Element::H],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
            ],
        );
        assert_atom_state(&methane, 0, 0, 0);

        let ammonia = run_perception(
            &[Element::N, Element::H, Element::H, Element::H],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        assert_atom_state(&ammonia, 0, 0, 1);
    }

    #[test]
    fn general_rules_handle_carbonyl_and_carbon_dioxide() {
        let formaldehyde = run_perception(
            &[Element::C, Element::O, Element::H, Element::H],
            &[
                (0, 1, BondOrder::Double),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        assert_atom_state(&formaldehyde, 0, 0, 0);
        assert_atom_state(&formaldehyde, 1, 0, 2);

        let carbon_dioxide = run_perception(
            &[Element::O, Element::C, Element::O],
            &[(0, 1, BondOrder::Double), (1, 2, BondOrder::Double)],
        );
        assert_atom_state(&carbon_dioxide, 1, 0, 0);
        assert_atom_state(&carbon_dioxide, 0, 0, 2);
        assert_atom_state(&carbon_dioxide, 2, 0, 2);
    }

    #[test]
    fn general_rules_handle_acetamide_fragment() {
        let elements = vec![
            Element::C,
            Element::O,
            Element::C,
            Element::N,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (2, 4, BondOrder::Single),
            (2, 5, BondOrder::Single),
            (2, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
            (3, 8, BondOrder::Single),
        ];

        let molecule = run_perception(&elements, &bonds);
        assert_atom_state(&molecule, 0, 0, 0);
        assert_atom_state(&molecule, 1, 0, 2);
        assert_atom_state(&molecule, 3, 0, 1);
        assert_atom_state(&molecule, 2, 0, 0);
    }
}
