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
/// * `molecule` - Annotated molecule containing degrees, lone pairs, and conjugation flags.
///
/// # Errors
///
/// Returns [`PerceptionError::HybridizationInference`] when an atom presents an unsupported
/// steric environment (e.g., steric numbers above four that cannot be explained by the
/// heuristics).
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

/// Chooses the best-fitting hybridization label for an atom and steric number pair.
///
/// The logic first filters out elements that never hybridize, then prioritizes resonance and
/// aromatic overrides before falling back to plain steric-number based VSEPR rules.
///
/// # Arguments
///
/// * `atom` - Annotated atom whose flags determine special-case handling.
/// * `steric_number` - Sum of σ-bonding domains and lone pairs computed upstream.
///
/// # Errors
///
/// Returns [`PerceptionError::HybridizationInference`] when no supported hybridization matches
/// the provided steric number, signalling inconsistent upstream annotations.
fn determine_hybridization(
    atom: &AnnotatedAtom,
    steric_number: u8,
) -> Result<Hybridization, PerceptionError> {
    if is_non_hybridized_element(atom.element) {
        return Ok(Hybridization::None);
    }

    if atom.is_resonant
        && !atom.is_anti_aromatic
        && (steric_number <= 3 || (steric_number == 4 && atom.lone_pairs > 0))
    {
        return Ok(Hybridization::Resonant);
    }

    match steric_number {
        4 => Ok(Hybridization::SP3),
        3 => Ok(Hybridization::SP2),
        2 => Ok(Hybridization::SP),
        0 | 1 => Ok(Hybridization::None),
        _ => Err(PerceptionError::HybridizationInference { atom_id: atom.id }),
    }
}

/// Detects elements that should stay in the `None` hybridization state regardless of geometry.
///
/// # Arguments
///
/// * `element` - Element under evaluation.
///
/// # Returns
///
/// `true` when the element belongs to the alkali, alkaline earth, halogen, noble gas, or
/// transition-metal sets that this perceiver treats as non-hybridizing.
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
