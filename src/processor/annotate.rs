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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{BondOrder, Element};

    fn build_and_process(
        num_atoms: usize,
        elements: &[(usize, Element)],
        bonds: &[(usize, usize, BondOrder)],
    ) -> ProcessingGraph {
        let mut mg = MolecularGraph::new();
        for _ in 0..num_atoms {
            mg.add_atom(Element::C);
        }
        for (i, el) in elements {
            mg.atoms[*i].element = *el;
        }
        for (u, v, order) in bonds {
            mg.add_bond(*u, *v, *order).unwrap();
        }

        process_graph(&mg).expect("Annotation pipeline failed during test setup")
    }

    #[test]
    fn test_ethane_annotation() {
        let pg = build_and_process(
            8,
            &[],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
                (1, 6, BondOrder::Single),
                (1, 7, BondOrder::Single),
            ],
        );

        let c1 = &pg.atoms[0];
        assert_eq!(c1.degree, 4);
        assert!(!c1.is_in_ring);
        assert!(!c1.is_aromatic);
        assert_eq!(c1.hybridization, Hybridization::SP3);
    }

    #[test]
    fn test_cyclohexane_annotation() {
        let pg = build_and_process(
            6,
            &[],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Single),
                (5, 0, BondOrder::Single),
            ],
        );

        for atom in &pg.atoms {
            assert_eq!(atom.degree, 2);
            assert!(atom.is_in_ring);
            assert_eq!(atom.smallest_ring_size, Some(6));
            assert!(!atom.is_aromatic);
            assert_eq!(atom.hybridization, Hybridization::SP3);
        }
    }

    #[test]
    fn test_benzene_annotation() {
        let pg = build_and_process(
            6,
            &[],
            &[
                (0, 1, BondOrder::Aromatic),
                (1, 2, BondOrder::Aromatic),
                (2, 3, BondOrder::Aromatic),
                (3, 4, BondOrder::Aromatic),
                (4, 5, BondOrder::Aromatic),
                (5, 0, BondOrder::Aromatic),
            ],
        );

        for atom in &pg.atoms {
            assert_eq!(atom.degree, 2);
            assert!(atom.is_in_ring);
            assert_eq!(atom.smallest_ring_size, Some(6));
            assert!(atom.is_aromatic);
            assert_eq!(atom.hybridization, Hybridization::Resonant);
        }
    }

    #[test]
    fn test_pyrrole_annotation() {
        let pg = build_and_process(
            5,
            &[(0, Element::N)],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 0, BondOrder::Single),
            ],
        );

        let nitrogen = &pg.atoms[0];
        assert_eq!(nitrogen.element, Element::N);
        assert_eq!(nitrogen.degree, 2);
        assert!(nitrogen.is_in_ring);
        assert_eq!(nitrogen.smallest_ring_size, Some(5));
        assert!(nitrogen.is_aromatic);
        assert_eq!(nitrogen.hybridization, Hybridization::Resonant);

        for i in 1..=4 {
            let carbon = &pg.atoms[i];
            assert_eq!(carbon.element, Element::C);
            assert!(carbon.is_in_ring);
            assert_eq!(carbon.smallest_ring_size, Some(5));
            assert!(carbon.is_aromatic);
            assert_eq!(carbon.hybridization, Hybridization::Resonant);
        }
    }

    #[test]
    fn test_boron_trifluoride_annotation() {
        let pg = build_and_process(
            4,
            &[
                (0, Element::B),
                (1, Element::F),
                (2, Element::F),
                (3, Element::F),
            ],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );

        let boron = &pg.atoms[0];
        assert_eq!(boron.element, Element::B);
        assert_eq!(boron.degree, 3);
        assert!(!boron.is_in_ring);
        assert!(!boron.is_aromatic);
        assert_eq!(boron.hybridization, Hybridization::SP2);

        for i in 1..=3 {
            let fluorine = &pg.atoms[i];
            assert_eq!(fluorine.element, Element::F);
            assert_eq!(fluorine.degree, 1);
            assert_eq!(fluorine.hybridization, Hybridization::None);
        }
    }
}
