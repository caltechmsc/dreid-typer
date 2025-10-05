use crate::core::Hybridization;
use crate::core::error::TyperError;
use crate::core::graph::{
    Angle, Atom, Bond, ImproperDihedral, MolecularGraph, MolecularTopology, ProperDihedral,
};
use crate::processor::ProcessingGraph;
use std::collections::HashSet;

pub(crate) fn build_topology(
    initial_graph: &MolecularGraph,
    processing_graph: &ProcessingGraph,
    atom_types: &[String],
) -> Result<MolecularTopology, TyperError> {
    let atoms = build_atoms(processing_graph, atom_types);
    let bonds = build_bonds(initial_graph);
    let angles = build_angles(processing_graph);
    let proper_dihedrals = build_proper_dihedrals(processing_graph);
    let improper_dihedrals = build_improper_dihedrals(processing_graph);

    Ok(MolecularTopology {
        atoms,
        bonds: bonds.into_iter().collect(),
        angles: angles.into_iter().collect(),
        proper_dihedrals: proper_dihedrals.into_iter().collect(),
        improper_dihedrals: improper_dihedrals.into_iter().collect(),
    })
}

fn build_atoms(processing_graph: &ProcessingGraph, atom_types: &[String]) -> Vec<Atom> {
    processing_graph
        .atoms
        .iter()
        .map(|atom_view| Atom {
            id: atom_view.id,
            element: atom_view.element,
            atom_type: atom_types[atom_view.id].clone(),
            hybridization: atom_view.hybridization,
        })
        .collect()
}

fn build_bonds(initial_graph: &MolecularGraph) -> HashSet<Bond> {
    initial_graph
        .bonds
        .iter()
        .map(|edge| {
            let (u, v) = edge.atom_ids;
            Bond::new(u, v, edge.order)
        })
        .collect()
}

fn build_angles(graph: &ProcessingGraph) -> HashSet<Angle> {
    let mut angles = HashSet::new();
    for j in 0..graph.atoms.len() {
        let neighbors = &graph.adjacency[j];
        if neighbors.len() < 2 {
            continue;
        }

        for i_idx in 0..neighbors.len() {
            for k_idx in (i_idx + 1)..neighbors.len() {
                let i = neighbors[i_idx].0;
                let k = neighbors[k_idx].0;
                angles.insert(Angle::new(i, j, k));
            }
        }
    }
    angles
}

fn build_proper_dihedrals(graph: &ProcessingGraph) -> HashSet<ProperDihedral> {
    let mut dihedrals = HashSet::new();
    for j in 0..graph.atoms.len() {
        for &(k, _) in &graph.adjacency[j] {
            if j >= k {
                continue;
            }

            for &(i, _) in &graph.adjacency[j] {
                if i == k {
                    continue;
                }

                for &(l, _) in &graph.adjacency[k] {
                    if l == j {
                        continue;
                    }

                    dihedrals.insert(ProperDihedral::new(i, j, k, l));
                }
            }
        }
    }
    dihedrals
}

fn build_improper_dihedrals(graph: &ProcessingGraph) -> HashSet<ImproperDihedral> {
    let mut dihedrals = HashSet::new();
    for i in 0..graph.atoms.len() {
        let atom = &graph.atoms[i];

        if atom.degree == 3 {
            if matches!(
                atom.hybridization,
                Hybridization::SP2 | Hybridization::Resonant
            ) {
                let neighbors = &graph.adjacency[i];
                let j = neighbors[0].0;
                let k = neighbors[1].0;
                let l = neighbors[2].0;

                dihedrals.insert(ImproperDihedral::new(j, k, i, l));
            }
        }
    }
    dihedrals
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::AtomNode;
    use crate::core::{BondOrder, Element};
    use crate::processor::{AtomView, ProcessingGraph};

    struct TestMolecule {
        initial_graph: MolecularGraph,
        processing_graph: ProcessingGraph,
        atom_types: Vec<String>,
    }

    impl TestMolecule {
        fn new(atoms: &[(Element, Hybridization, &str)]) -> Self {
            let mut initial_graph = MolecularGraph::new();
            let mut atom_views = Vec::new();
            let mut atom_types = Vec::new();

            for (i, (element, hybridization, atom_type)) in atoms.iter().enumerate() {
                initial_graph.atoms.push(AtomNode {
                    id: i,
                    element: *element,
                });
                atom_views.push(AtomView {
                    id: i,
                    element: *element,
                    degree: 0,
                    hybridization: *hybridization,
                    is_in_ring: false,
                    smallest_ring_size: None,
                    is_aromatic: *hybridization == Hybridization::Resonant,
                });
                atom_types.push(atom_type.to_string());
            }

            Self {
                initial_graph,
                processing_graph: ProcessingGraph {
                    atoms: atom_views,
                    adjacency: vec![],
                },
                atom_types,
            }
        }

        fn with_bond(mut self, u: usize, v: usize, order: BondOrder) -> Self {
            self.initial_graph.add_bond(u, v, order).unwrap();
            self
        }

        fn build(mut self) -> Self {
            let num_atoms = self.initial_graph.atoms.len();
            let mut adjacency = vec![vec![]; num_atoms];
            for bond in &self.initial_graph.bonds {
                let (u, v) = bond.atom_ids;
                adjacency[u].push((v, bond.order));
                adjacency[v].push((u, bond.order));
            }

            for i in 0..num_atoms {
                self.processing_graph.atoms[i].degree = adjacency[i].len() as u8;
            }
            self.processing_graph.adjacency = adjacency;

            self
        }
    }

    #[test]
    fn test_build_methane_topology() {
        let molecule = TestMolecule::new(&[
            (Element::C, Hybridization::SP3, "C_3"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
        ])
        .with_bond(0, 1, BondOrder::Single)
        .with_bond(0, 2, BondOrder::Single)
        .with_bond(0, 3, BondOrder::Single)
        .with_bond(0, 4, BondOrder::Single)
        .build();

        let topology = build_topology(
            &molecule.initial_graph,
            &molecule.processing_graph,
            &molecule.atom_types,
        )
        .unwrap();

        assert_eq!(topology.atoms.len(), 5);
        assert_eq!(topology.bonds.len(), 4);
        assert_eq!(topology.angles.len(), 6);
        assert_eq!(topology.proper_dihedrals.len(), 0);
        assert_eq!(topology.improper_dihedrals.len(), 0);
    }

    #[test]
    fn test_build_ethane_topology() {
        let molecule = TestMolecule::new(&[
            (Element::C, Hybridization::SP3, "C_3"),
            (Element::C, Hybridization::SP3, "C_3"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
        ])
        .with_bond(0, 1, BondOrder::Single)
        .with_bond(0, 2, BondOrder::Single)
        .with_bond(0, 3, BondOrder::Single)
        .with_bond(0, 4, BondOrder::Single)
        .with_bond(1, 5, BondOrder::Single)
        .with_bond(1, 6, BondOrder::Single)
        .with_bond(1, 7, BondOrder::Single)
        .build();

        let topology = build_topology(
            &molecule.initial_graph,
            &molecule.processing_graph,
            &molecule.atom_types,
        )
        .unwrap();

        assert_eq!(topology.atoms.len(), 8);
        assert_eq!(topology.bonds.len(), 7);
        assert_eq!(topology.angles.len(), 12);
        assert_eq!(topology.proper_dihedrals.len(), 9);
        assert_eq!(topology.improper_dihedrals.len(), 0);
    }

    #[test]
    fn test_build_ethylene_topology() {
        let molecule = TestMolecule::new(&[
            (Element::C, Hybridization::SP2, "C_2"),
            (Element::C, Hybridization::SP2, "C_2"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
        ])
        .with_bond(0, 1, BondOrder::Double)
        .with_bond(0, 2, BondOrder::Single)
        .with_bond(0, 3, BondOrder::Single)
        .with_bond(1, 4, BondOrder::Single)
        .with_bond(1, 5, BondOrder::Single)
        .build();

        let topology = build_topology(
            &molecule.initial_graph,
            &molecule.processing_graph,
            &molecule.atom_types,
        )
        .unwrap();

        assert_eq!(topology.bonds.len(), 5);
        assert_eq!(topology.angles.len(), 6);
        assert_eq!(topology.proper_dihedrals.len(), 4);
        assert_eq!(topology.improper_dihedrals.len(), 2);

        assert!(
            topology
                .improper_dihedrals
                .contains(&ImproperDihedral::new(1, 2, 0, 3))
        );
        assert!(
            topology
                .improper_dihedrals
                .contains(&ImproperDihedral::new(0, 4, 1, 5))
        );
    }

    #[test]
    fn test_build_formaldehyde_topology() {
        let molecule = TestMolecule::new(&[
            (Element::C, Hybridization::SP2, "C_2"),
            (Element::O, Hybridization::SP2, "O_2"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
        ])
        .with_bond(0, 1, BondOrder::Double)
        .with_bond(0, 2, BondOrder::Single)
        .with_bond(0, 3, BondOrder::Single)
        .build();

        let topology = build_topology(
            &molecule.initial_graph,
            &molecule.processing_graph,
            &molecule.atom_types,
        )
        .unwrap();

        assert_eq!(topology.bonds.len(), 3);
        assert_eq!(topology.angles.len(), 3);
        assert_eq!(topology.proper_dihedrals.len(), 0);
        assert_eq!(topology.improper_dihedrals.len(), 1);
        assert!(
            topology
                .improper_dihedrals
                .contains(&ImproperDihedral::new(1, 2, 0, 3))
        );
    }

    #[test]
    fn test_build_benzene_topology() {
        let mut molecule_builder = TestMolecule::new(&[
            (Element::C, Hybridization::Resonant, "C_R"),
            (Element::C, Hybridization::Resonant, "C_R"),
            (Element::C, Hybridization::Resonant, "C_R"),
            (Element::C, Hybridization::Resonant, "C_R"),
            (Element::C, Hybridization::Resonant, "C_R"),
            (Element::C, Hybridization::Resonant, "C_R"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
        ]);

        for i in 0..6 {
            molecule_builder = molecule_builder.with_bond(i, (i + 1) % 6, BondOrder::Aromatic);
        }
        for i in 0..6 {
            molecule_builder = molecule_builder.with_bond(i, i + 6, BondOrder::Single);
        }

        let molecule = molecule_builder.build();

        let topology = build_topology(
            &molecule.initial_graph,
            &molecule.processing_graph,
            &molecule.atom_types,
        )
        .unwrap();

        let num_atoms = 12;
        let num_bonds = 12;
        let num_angles = 18;
        let num_proper_dihedrals = 24;
        let num_improper_dihedrals = 6;

        assert_eq!(topology.atoms.len(), num_atoms);
        assert_eq!(topology.bonds.len(), num_bonds);
        assert_eq!(topology.angles.len(), num_angles);
        assert!(
            topology.proper_dihedrals.len() > 0,
            "Expected some proper dihedrals"
        );
        assert_eq!(topology.improper_dihedrals.len(), num_improper_dihedrals);

        assert!(
            topology
                .improper_dihedrals
                .contains(&ImproperDihedral::new(1, 5, 0, 6))
        );
    }

    #[test]
    fn test_build_diatomic_no_angles_or_dihedrals() {
        let molecule = TestMolecule::new(&[
            (Element::H, Hybridization::None, "H_"),
            (Element::H, Hybridization::None, "H_"),
        ])
        .with_bond(0, 1, BondOrder::Single)
        .build();

        let topology = build_topology(
            &molecule.initial_graph,
            &molecule.processing_graph,
            &molecule.atom_types,
        )
        .unwrap();

        assert_eq!(topology.bonds.len(), 1);
        assert_eq!(topology.angles.len(), 0);
        assert_eq!(topology.proper_dihedrals.len(), 0);
        assert_eq!(topology.improper_dihedrals.len(), 0);
    }
}
