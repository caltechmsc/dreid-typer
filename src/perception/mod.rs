mod aromaticity;
mod electrons;
mod hybridization;
mod kekulize;
mod model;
mod resonance;
mod rings;

pub use model::AnnotatedMolecule;

use crate::core::error::{PerceptionError, TyperError};
use crate::core::graph::MolecularGraph;

pub fn perceive(graph: &MolecularGraph) -> Result<AnnotatedMolecule, TyperError> {
    let mut molecule = AnnotatedMolecule::new(graph).map_err(TyperError::InvalidInput)?;

    let pipeline: [(
        &str,
        fn(&mut AnnotatedMolecule) -> Result<(), PerceptionError>,
    ); 6] = [
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
