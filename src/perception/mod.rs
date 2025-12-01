//! Coordinates the sequential perception pipeline that annotates molecules prior to typing.
//!
//! This module wires the specialized perception stages—ring detection, Kekulé expansion,
//! electron bookkeeping, aromaticity, resonance, and hybridization—into a single pass that
//! populates an [`AnnotatedMolecule`] for downstream typing.

mod aromaticity;
mod electrons;
mod hybridization;
mod kekulize;
mod model;
mod resonance;
mod rings;

pub use model::{AnnotatedAtom, AnnotatedMolecule, ResonanceSystem};

use crate::core::error::{PerceptionError, TyperError};
use crate::core::graph::MolecularGraph;

type PerceptionStepFn = fn(&mut AnnotatedMolecule) -> Result<(), PerceptionError>;
type PerceptionStep = (&'static str, PerceptionStepFn);

/// Runs the full perception pipeline and returns an annotated molecule.
///
/// The function constructs an [`AnnotatedMolecule`] from the input graph, executes the fixed set
/// of perception stages in order, and reports any failure with the offending step name folded
/// into the [`TyperError`].
///
/// # Arguments
///
/// * `graph` - Validated molecular graph containing atoms and bonds.
///
/// # Returns
///
/// Fully annotated molecule that records ring membership, electron bookkeeping, aromaticity,
/// resonance, and hybridization properties for every atom.
///
/// # Errors
///
/// Returns [`TyperError::InvalidInput`] when the graph contains invalid bonding, or
/// [`TyperError::PerceptionFailed`] when any perception stage emits a [`PerceptionError`].
pub fn perceive(graph: &MolecularGraph) -> Result<AnnotatedMolecule, TyperError> {
    let mut molecule = AnnotatedMolecule::new(graph).map_err(TyperError::InvalidInput)?;

    let pipeline: [PerceptionStep; 6] = [
        ("Rings", rings::perceive),
        ("Kekulization", kekulize::perceive),
        ("Electrons", electrons::perceive),
        ("Aromaticity", aromaticity::perceive),
        ("Resonance", resonance::perceive),
        ("Hybridization", hybridization::perceive),
    ];

    for (name, step_fn) in pipeline {
        step_fn(&mut molecule).map_err(|source| TyperError::PerceptionFailed {
            step: name.to_string(),
            source,
        })?;
    }

    Ok(molecule)
}
