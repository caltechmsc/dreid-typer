use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::{BondOrder, Element};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let mut processed = vec![false; molecule.atoms.len()];

    assign_general(molecule, &processed)?;

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
