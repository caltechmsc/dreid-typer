use super::graph::{ProcessingGraph, RingInfo};
use super::{aromaticity, hybridization, rings};
use crate::core::Hybridization;
use crate::core::error::{AnnotationError, TyperError};
use crate::core::graph::MolecularGraph;

pub(crate) fn process_graph(
    molecular_graph: &MolecularGraph,
) -> Result<ProcessingGraph, TyperError> {
    // --- Phase 1 & 2: Graph Construction and Base Annotation ---
    let mut graph = ProcessingGraph::new(molecular_graph).map_err(TyperError::InvalidInputGraph)?;

    // --- Phase 3.1: Ring System Perception ---
    let ring_info = rings::perceive_rings(&graph);
    apply_ring_info(&mut graph, &ring_info);

    // --- Phase 3.3: Aromaticity Perception ---
    aromaticity::perceive_aromaticity(&mut graph, &ring_info)?;

    // --- Phase 4: Final Hybridization Inference ---
    hybridization::infer_hybridization_for_all(&mut graph).map_err(TyperError::AnnotationFailed)?;

    // --- Final Validation ---
    if let Some(failed_atom) = graph
        .atoms
        .iter()
        .find(|a| a.hybridization == Hybridization::Unknown)
    {
        return Err(TyperError::AnnotationFailed(
            AnnotationError::HybridizationInference {
                atom_id: failed_atom.id,
            },
        ));
    }

    Ok(graph)
}

fn apply_ring_info(graph: &mut ProcessingGraph, ring_info: &RingInfo) {
    let mut atom_ring_sizes: Vec<Vec<u8>> = vec![vec![]; graph.atoms.len()];
    for ring in &ring_info.0 {
        let ring_len = ring.len() as u8;
        for &atom_id in ring {
            atom_ring_sizes[atom_id].push(ring_len);
        }
    }

    for (atom_id, atom) in graph.atoms.iter_mut().enumerate() {
        if !atom_ring_sizes[atom_id].is_empty() {
            atom.is_in_ring = true;
            atom.smallest_ring_size = atom_ring_sizes[atom_id].iter().min().copied();
        }
    }
}
