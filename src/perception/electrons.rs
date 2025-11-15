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

        if let (Some(c1), Some(o), Some(c2)) =
            (double_bond_c_idx, single_bond_o_idx, single_bond_c_idx)
        {
            if !processed[c1] && !processed[o] && !processed[c2] {
                molecule.atoms[n_idx].formal_charge = 1;
                molecule.atoms[n_idx].lone_pairs = 0;
                molecule.atoms[o].formal_charge = -1;
                molecule.atoms[o].lone_pairs = 3;

                processed[n_idx] = true;
                processed[o] = true;
            }
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

        if let (Some(o1), Some(o2)) = (double_bond_o_idx, single_bond_o_idx) {
            if other_neighbor_count == 1 && !processed[o1] && !processed[o2] {
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

        if let (Some(o1), Some(o2)) = (double_bond_o_idx, single_bond_o_idx) {
            if !processed[o1] && !processed[o2] {
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
        }
    }
    Ok(())
}

fn assign_ammonium_and_iminium(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx] || molecule.atoms[n_idx].element != Element::N {
            continue;
        }

        let atom = &molecule.atoms[n_idx];
        let has_double_bond = molecule.adjacency[n_idx]
            .iter()
            .any(|&(_, order)| order == BondOrder::Double);

        if atom.degree == 4 {
            molecule.atoms[n_idx].formal_charge = 1;
            molecule.atoms[n_idx].lone_pairs = 0;
            processed[n_idx] = true;
        } else if atom.degree == 3 && has_double_bond {
            molecule.atoms[n_idx].formal_charge = 1;
            molecule.atoms[n_idx].lone_pairs = 0;
            processed[n_idx] = true;
        }
    }
    Ok(())
}

fn assign_onium_ions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for idx in 0..molecule.atoms.len() {
        if processed[idx] {
            continue;
        }

        let atom = &molecule.atoms[idx];
        if (atom.element == Element::O || atom.element == Element::S) && atom.degree == 3 {
            molecule.atoms[idx].formal_charge = 1;
            molecule.atoms[idx].lone_pairs = 1;
            processed[idx] = true;
        }
    }
    Ok(())
}

fn assign_phosphonium_ions(
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

        let has_double_bond_o = molecule.adjacency[p_idx].iter().any(|&(id, order)| {
            molecule.atoms[id].element == Element::O && order == BondOrder::Double
        });

        if !has_double_bond_o {
            molecule.atoms[p_idx].formal_charge = 1;
            molecule.atoms[p_idx].lone_pairs = 0;
            processed[p_idx] = true;
        }
    }
    Ok(())
}

fn assign_general(
    molecule: &mut AnnotatedMolecule,
    processed: &[bool],
) -> Result<(), PerceptionError> {
    for i in 0..molecule.atoms.len() {
        if processed[i] {
            continue;
        }
        let atom = &molecule.atoms[i];

        let valence = atom.element.valence_electrons().ok_or_else(|| {
            PerceptionError::Other(format!(
                "valence electrons not defined for element {:?}",
                atom.element
            ))
        })?;

        let bonding_electrons: u8 = molecule.adjacency[i]
            .iter()
            .map(|&(_, order)| bond_order_to_valence(order))
            .sum();

        let mut lone_pairs = 0;

        let is_second_period = matches!(
            atom.element,
            Element::B | Element::C | Element::N | Element::O | Element::F
        );

        if atom.element == Element::H {
            let bonded_electrons = bonding_electrons.saturating_mul(2);
            if bonded_electrons <= 2 {
                lone_pairs = (2 - bonded_electrons) / 2;
            }
        } else if is_second_period {
            let bonded_electrons = bonding_electrons.saturating_mul(2);
            if bonded_electrons <= 8 {
                lone_pairs = (8 - bonded_electrons) / 2;
            }
        } else {
            if valence >= bonding_electrons {
                lone_pairs = (valence - bonding_electrons) / 2;
            }
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
