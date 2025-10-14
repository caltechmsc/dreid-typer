//! Orchestrates the molecular perception pipeline for the DREIDING force field.
//!
//! This module coordinates the three-phase perception process: electron counting, functional group
//! template application, and generic property detection. It serves as the main entry point for
//! converting a `MolecularGraph` into a richly annotated `ProcessingGraph` ready for atom typing.

use super::graph::{ProcessingGraph, RingInfo};
use super::perception::{
    apply_ring_annotations, perceive_electron_counts, perceive_generic_properties, perceive_rings,
};
use super::templates;
use crate::core::error::TyperError;
use crate::core::graph::MolecularGraph;

/// Result of the molecular perception pipeline.
///
/// Contains the fully annotated processing graph and ring information extracted during perception.
#[derive(Debug)]
pub struct PerceptionResult {
    /// The processing graph with all chemical annotations applied.
    pub processing_graph: ProcessingGraph,
    /// Information about all detected rings in the molecule.
    pub ring_info: RingInfo,
}

/// Runs the complete molecular perception pipeline.
///
/// Executes the three-phase perception process: electron counting, functional group templates,
/// and generic property detection to produce a fully annotated processing graph.
///
/// # Arguments
///
/// * `molecular_graph` - The input molecular graph to perceive.
///
/// # Returns
///
/// A `PerceptionResult` containing the annotated processing graph and ring information.
///
/// # Errors
///
/// Returns `TyperError::InvalidInputGraph` if graph conversion fails, or `TyperError::AnnotationFailed`
/// if any perception step encounters an error.
pub fn perceive(molecular_graph: &MolecularGraph) -> Result<PerceptionResult, TyperError> {
    // Phase 1: Electron count perception
    let mut processing_graph = perceive_electron_counts(molecular_graph)?;

    // Phase 2: Functional group template perception (expert system)
    if let Err(err) = templates::apply_functional_group_templates(&mut processing_graph) {
        return Err(TyperError::AnnotationFailed(err));
    }

    // Phase 3: Ring system and generic aromaticity perception (general system)
    let ring_info = perceive_rings(&processing_graph);
    apply_ring_annotations(&mut processing_graph, &ring_info);
    if let Err(err) = perceive_generic_properties(&mut processing_graph, &ring_info) {
        return Err(TyperError::AnnotationFailed(err));
    }

    Ok(PerceptionResult {
        processing_graph,
        ring_info,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element, Hybridization};

    #[test]
    fn perceiving_methane_yields_expected_annotations() {
        let mut mg = MolecularGraph::new();
        let c = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);

        mg.add_bond(c, h1, BondOrder::Single).unwrap();
        mg.add_bond(c, h2, BondOrder::Single).unwrap();
        mg.add_bond(c, h3, BondOrder::Single).unwrap();
        mg.add_bond(c, h4, BondOrder::Single).unwrap();

        let result = perceive(&mg).unwrap();

        assert_eq!(result.processing_graph.atoms.len(), 5);
        assert!(result.ring_info.0.is_empty());

        let carbon = &result.processing_graph.atoms[c];
        assert_eq!(carbon.steric_number, 4);
        assert_eq!(carbon.hybridization, Hybridization::SP3);
        assert!(!carbon.is_in_ring);
        assert!(!carbon.is_aromatic);
    }

    #[test]
    fn perceiving_benzene_marks_aromatic_ring() {
        let mut mg = MolecularGraph::new();
        let mut atoms = vec![];
        for _ in 0..6 {
            atoms.push(mg.add_atom(Element::C, 0));
        }
        for i in 0..6 {
            mg.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }

        let result = perceive(&mg).unwrap();

        assert_eq!(result.ring_info.0.len(), 1);
        assert!(
            result
                .processing_graph
                .atoms
                .iter()
                .all(|atom| atom.is_aromatic)
        );
        assert!(
            result
                .processing_graph
                .atoms
                .iter()
                .all(|atom| atom.hybridization == Hybridization::Resonant)
        );
    }
}
