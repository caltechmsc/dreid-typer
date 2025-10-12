use super::graph::{ProcessingGraph, RingInfo};
use super::perception::{
    apply_ring_annotations, perceive_electron_counts, perceive_generic_properties, perceive_rings,
};
use super::templates;
use crate::core::error::TyperError;
use crate::core::graph::MolecularGraph;

#[derive(Debug)]
pub struct PerceptionResult {
    pub processing_graph: ProcessingGraph,
    pub ring_info: RingInfo,
}

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
