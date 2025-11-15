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

type PerceptionStepFn = fn(&mut AnnotatedMolecule) -> Result<(), PerceptionError>;
type PerceptionStep = (&'static str, PerceptionStepFn);

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::properties::{BondOrder, Element, Hybridization};

    fn benzene_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let carbons: Vec<_> = (0..6).map(|_| graph.add_atom(Element::C)).collect();
        let hydrogens: Vec<_> = (0..6).map(|_| graph.add_atom(Element::H)).collect();

        let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for &(u, v) in &ring_edges {
            graph
                .add_bond(carbons[u], carbons[v], BondOrder::Aromatic)
                .expect("valid aromatic bond in benzene ring");
        }
        for i in 0..6 {
            graph
                .add_bond(carbons[i], hydrogens[i], BondOrder::Single)
                .expect("valid C-H bond");
        }

        graph
    }

    fn aromatic_bond_outside_ring_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let c1 = graph.add_atom(Element::C);
        let c2 = graph.add_atom(Element::C);
        let h1 = graph.add_atom(Element::H);
        let h2 = graph.add_atom(Element::H);

        graph
            .add_bond(c1, c2, BondOrder::Aromatic)
            .expect("valid aromatic bond");
        graph
            .add_bond(c1, h1, BondOrder::Single)
            .expect("valid C-H bond");
        graph
            .add_bond(c2, h2, BondOrder::Single)
            .expect("valid C-H bond");

        graph
    }

    #[test]
    fn perception_pipeline_assigns_benzene_properties() {
        let graph = benzene_graph();
        let molecule = perceive(&graph).expect("perception pipeline should succeed");

        assert_eq!(molecule.rings.len(), 1, "benzene must yield a single ring");
        for (idx, atom) in molecule.atoms.iter().enumerate() {
            match atom.element {
                Element::C => {
                    assert!(atom.is_in_ring, "carbon {idx} must be in the ring");
                    assert!(atom.is_aromatic, "carbon {idx} must be aromatic");
                    assert!(
                        atom.is_in_conjugated_system,
                        "carbon {idx} must be in conjugated system"
                    );
                    assert_eq!(
                        atom.hybridization,
                        Hybridization::Resonant,
                        "carbon {idx} should end Resonant"
                    );
                    assert_eq!(atom.steric_number, 3);
                }
                Element::H => {
                    assert_eq!(atom.hybridization, Hybridization::None);
                    assert_eq!(atom.steric_number, 1);
                }
                other => panic!("unexpected element in benzene fixture: {other:?}"),
            }
        }

        assert!(
            molecule
                .bonds
                .iter()
                .all(|bond| bond.order != BondOrder::Aromatic),
            "all aromatic bonds should be Kekulé-expanded"
        );
    }

    #[test]
    fn pipeline_reports_step_name_when_kekulization_fails() {
        let graph = aromatic_bond_outside_ring_graph();
        let err = perceive(&graph).expect_err("pipeline should fail before completion");

        match err {
            TyperError::PerceptionFailed { step, source } => {
                assert_eq!(step, "Kekulization");
                match source {
                    PerceptionError::KekulizationFailed { message } => {
                        assert!(
                            message.contains("not in a ring"),
                            "error message should describe the missing ring context"
                        );
                    }
                    other => panic!("unexpected perception error emitted: {other:?}"),
                }
            }
            other => panic!("expected PerceptionFailed error, got {other:?}"),
        }
    }
}
