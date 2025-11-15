use super::model::{AnnotatedAtom, AnnotatedMolecule};
use crate::core::error::PerceptionError;
use crate::core::properties::{Element, Hybridization};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    for i in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[i];

        let steric_number = atom.degree + atom.lone_pairs;

        let hybridization = determine_hybridization(atom, steric_number)?;

        let atom_mut = &mut molecule.atoms[i];
        atom_mut.hybridization = hybridization;

        atom_mut.steric_number = match hybridization {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            _ => steric_number,
        };
    }
    Ok(())
}

fn determine_hybridization(
    atom: &AnnotatedAtom,
    steric_number: u8,
) -> Result<Hybridization, PerceptionError> {
    if is_non_hybridized_element(atom.element) {
        return Ok(Hybridization::None);
    }

    if atom.is_in_conjugated_system && !atom.is_anti_aromatic {
        if steric_number <= 3 || (steric_number == 4 && atom.lone_pairs > 0) {
            return Ok(Hybridization::Resonant);
        }
    }

    if atom.is_aromatic {
        return Ok(Hybridization::SP2);
    }

    match steric_number {
        4 => Ok(Hybridization::SP3),
        3 => Ok(Hybridization::SP2),
        2 => Ok(Hybridization::SP),
        0 | 1 => Ok(Hybridization::None),
        _ => Err(PerceptionError::HybridizationInference { atom_id: atom.id }),
    }
}

fn is_non_hybridized_element(element: Element) -> bool {
    matches!(
        element,
        Element::Li
            | Element::Na
            | Element::K
            | Element::Rb
            | Element::Cs
            | Element::Fr
            | Element::Be
            | Element::Mg
            | Element::Ca
            | Element::Sr
            | Element::Ba
            | Element::Ra
            | Element::F
            | Element::Cl
            | Element::Br
            | Element::I
            | Element::At
            | Element::He
            | Element::Ne
            | Element::Ar
            | Element::Kr
            | Element::Xe
            | Element::Rn
            | Element::Sc
            | Element::Ti
            | Element::V
            | Element::Cr
            | Element::Mn
            | Element::Fe
            | Element::Co
            | Element::Ni
            | Element::Cu
            | Element::Zn
            | Element::Y
            | Element::Zr
            | Element::Nb
            | Element::Mo
            | Element::Tc
            | Element::Ru
            | Element::Rh
            | Element::Pd
            | Element::Ag
            | Element::Cd
            | Element::Hf
            | Element::Ta
            | Element::W
            | Element::Re
            | Element::Os
            | Element::Ir
            | Element::Pt
            | Element::Au
            | Element::Hg
    )
}
