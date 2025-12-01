//! Assigns VSEPR-consistent hybridization states after conjugation and aromaticity analysis.
//!
//! This perception stage translates the structural annotations produced by earlier passes
//! (degree, lone pairs, conjugation, aromatic flags) into concrete `Hybridization` labels and
//! the corresponding steric numbers required by later typing decisions.

use super::model::{AnnotatedAtom, AnnotatedMolecule};
use crate::core::error::PerceptionError;
use crate::core::properties::{Element, Hybridization};

/// Updates every atom with its final hybridization and steric number assignments.
///
/// The procedure evaluates the VSEPR steric number for each annotated atom, reconciles it
/// with conjugation and aromatic flags, and records the resulting `Hybridization` along with
/// the adjusted steric number in-place.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule containing degrees, lone pairs, and resonance flags.
///
/// # Errors
///
/// Returns [`PerceptionError::HybridizationInference`] when an atom presents an unsupported
/// steric environment.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    for atom in &mut molecule.atoms {
        atom.hybridization = initial_hybridization(atom)?;
    }

    loop {
        let mut changes = 0;
        for i in 0..molecule.atoms.len() {
            if molecule.atoms[i].hybridization == Hybridization::SP3
                && molecule.atoms[i].lone_pairs > 0
            {
                let is_adjacent_to_pi_system =
                    molecule.adjacency[i].iter().any(|&(neighbor_id, _)| {
                        matches!(
                            molecule.atoms[neighbor_id].hybridization,
                            Hybridization::SP2 | Hybridization::SP | Hybridization::Resonant
                        )
                    });

                if is_adjacent_to_pi_system {
                    molecule.atoms[i].hybridization = Hybridization::Resonant;
                    molecule.atoms[i].is_resonant = true;
                    changes += 1;
                }
            }
        }
        if changes == 0 {
            break;
        }
    }

    for atom in &mut molecule.atoms {
        atom.steric_number = match atom.hybridization {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            Hybridization::None | Hybridization::Unknown => atom.degree + atom.lone_pairs,
        };
    }

    Ok(())
}

/// Determines the initial hybridization for a given atom, respecting resonance flags
/// before applying pure VSEPR steric-number logic.
fn initial_hybridization(atom: &AnnotatedAtom) -> Result<Hybridization, PerceptionError> {
    if is_non_hybridized_element(atom.element) {
        return Ok(Hybridization::None);
    }

    if atom.is_resonant && !atom.is_anti_aromatic {
        return Ok(Hybridization::Resonant);
    }

    let steric_number = atom.degree + atom.lone_pairs;
    match steric_number {
        4 => Ok(Hybridization::SP3),
        3 => Ok(Hybridization::SP2),
        2 => Ok(Hybridization::SP),
        0 | 1 => Ok(Hybridization::None),
        _ => Err(PerceptionError::HybridizationInference { atom_id: atom.id }),
    }
}

/// Detects elements that should stay in the `None` hybridization state regardless of geometry.
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
