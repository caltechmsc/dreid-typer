use super::model::{AnnotatedAtom, AnnotatedMolecule, Ring};
use crate::core::error::PerceptionError;
use crate::core::properties::BondOrder;
use std::collections::{HashMap, HashSet};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    if molecule.rings.is_empty() {
        return Ok(());
    }

    let ring_systems_indices = find_ring_systems(&molecule.rings);

    for system_indices in ring_systems_indices {
        let system_atoms: HashSet<usize> = system_indices
            .iter()
            .flat_map(|&i| molecule.rings[i].iter())
            .copied()
            .collect();

        let model = AromaticityModel::new(molecule, &system_atoms);

        if model.is_aromatic() {
            for &atom_id in &system_atoms {
                molecule.atoms[atom_id].is_aromatic = true;
            }
        } else if model.is_anti_aromatic() {
            for &atom_id in &system_atoms {
                molecule.atoms[atom_id].is_anti_aromatic = true;
            }
        }
    }

    Ok(())
}

struct AromaticityModel<'a> {
    molecule: &'a AnnotatedMolecule,
    atoms: HashSet<usize>,
    pi_electrons: Option<u32>,
    is_potentially_planar: bool,
}

impl<'a> AromaticityModel<'a> {
    fn new(molecule: &'a AnnotatedMolecule, system_atoms: &HashSet<usize>) -> Self {
        let mut model = Self {
            molecule,
            atoms: system_atoms.clone(),
            pi_electrons: None,
            is_potentially_planar: false,
        };
        model.evaluate();
        model
    }

    fn is_aromatic(&self) -> bool {
        if !self.is_potentially_planar {
            return false;
        }
        matches!(self.pi_electrons, Some(pi) if pi > 0 && (pi - 2) % 4 == 0)
    }

    fn is_anti_aromatic(&self) -> bool {
        if !self.is_potentially_planar {
            return false;
        }
        matches!(self.pi_electrons, Some(pi) if pi > 0 && pi % 4 == 0)
    }

    fn evaluate(&mut self) {
        if !self
            .atoms
            .iter()
            .all(|&id| is_potentially_planar(&self.molecule.atoms[id]))
        {
            self.is_potentially_planar = false;
            return;
        }
        self.is_potentially_planar = true;

        let mut pi_count = 0;
        for &atom_id in &self.atoms {
            if let Some(contribution) = self.count_pi_contribution(atom_id) {
                pi_count += contribution;
            } else {
                self.is_potentially_planar = false;
                self.pi_electrons = None;
                return;
            }
        }
        self.pi_electrons = Some(pi_count);
    }

    fn count_pi_contribution(&self, atom_id: usize) -> Option<u32> {
        let atom = &self.molecule.atoms[atom_id];

        let has_endocyclic_double_bond = self.molecule.adjacency[atom_id]
            .iter()
            .any(|&(n_id, order)| order == BondOrder::Double && self.atoms.contains(&n_id));

        let has_exocyclic_double_bond = self.molecule.adjacency[atom_id]
            .iter()
            .any(|&(n_id, order)| order == BondOrder::Double && !self.atoms.contains(&n_id));

        if has_endocyclic_double_bond {
            return Some(1);
        }

        if !has_exocyclic_double_bond && atom.lone_pairs > 0 {
            return Some(2);
        }

        if atom.formal_charge == -1 {
            return Some(2);
        }
        if atom.formal_charge == 1 {
            return Some(0);
        }

        if has_exocyclic_double_bond {
            return Some(1);
        }

        None
    }
}

fn is_potentially_planar(atom: &AnnotatedAtom) -> bool {
    let steric_number = atom.degree + atom.lone_pairs;
    match steric_number {
        0..=3 => true,
        4 => atom.lone_pairs > 0,
        _ => false,
    }
}

fn find_ring_systems(rings: &[Ring]) -> Vec<Vec<usize>> {
    if rings.is_empty() {
        return vec![];
    }

    let ring_adj = build_ring_adjacency(rings);
    let mut systems = Vec::new();
    let mut visited = vec![false; rings.len()];

    for i in 0..rings.len() {
        if !visited[i] {
            let mut current_system_indices = Vec::new();
            let mut stack = vec![i];
            visited[i] = true;

            while let Some(ring_idx) = stack.pop() {
                current_system_indices.push(ring_idx);
                for &neighbor_idx in &ring_adj[ring_idx] {
                    if !visited[neighbor_idx] {
                        visited[neighbor_idx] = true;
                        stack.push(neighbor_idx);
                    }
                }
            }
            systems.push(current_system_indices);
        }
    }
    systems
}

fn build_ring_adjacency(rings: &[Ring]) -> Vec<Vec<usize>> {
    let mut atom_to_rings: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ring_idx, ring) in rings.iter().enumerate() {
        for &atom_id in ring {
            atom_to_rings.entry(atom_id).or_default().push(ring_idx);
        }
    }

    let mut adj = vec![vec![]; rings.len()];
    for ring_indices in atom_to_rings.values() {
        if ring_indices.len() > 1 {
            for i in 0..ring_indices.len() {
                for j in (i + 1)..ring_indices.len() {
                    let r1 = ring_indices[i];
                    let r2 = ring_indices[j];
                    adj[r1].push(r2);
                    adj[r2].push(r1);
                }
            }
        }
    }

    for neighbors in adj.iter_mut() {
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    adj
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::properties::Element;

    #[derive(Clone, Copy)]
    struct AtomSpec {
        element: Element,
        formal_charge: i8,
        lone_pairs: u8,
        degree_override: Option<u8>,
    }

    impl AtomSpec {
        fn new(element: Element) -> Self {
            Self {
                element,
                formal_charge: 0,
                lone_pairs: 0,
                degree_override: None,
            }
        }

        fn with_charge(mut self, charge: i8) -> Self {
            self.formal_charge = charge;
            self
        }

        fn with_lone_pairs(mut self, lone_pairs: u8) -> Self {
            self.lone_pairs = lone_pairs;
            self
        }

        fn with_degree(mut self, degree: u8) -> Self {
            self.degree_override = Some(degree);
            self
        }
    }

    fn build_test_molecule(
        atom_specs: &[AtomSpec],
        bonds: &[(usize, usize, BondOrder)],
        rings: &[&[usize]],
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for spec in atom_specs {
            graph.add_atom(spec.element);
        }
        for &(a, b, order) in bonds {
            graph.add_bond(a, b, order).expect("valid bond definition");
        }

        let mut molecule = AnnotatedMolecule::new(&graph).expect("graph must be valid");
        molecule.rings = rings.iter().map(|ring| ring.to_vec()).collect();
        annotate_ring_flags(&mut molecule);
        apply_atom_specs(&mut molecule, atom_specs);
        molecule
    }

    fn annotate_ring_flags(molecule: &mut AnnotatedMolecule) {
        for ring in &molecule.rings {
            for &atom_id in ring {
                let atom = &mut molecule.atoms[atom_id];
                atom.is_in_ring = true;
            }
        }
    }

    fn apply_atom_specs(molecule: &mut AnnotatedMolecule, specs: &[AtomSpec]) {
        for (i, spec) in specs.iter().enumerate() {
            let atom = &mut molecule.atoms[i];
            atom.formal_charge = spec.formal_charge;
            atom.lone_pairs = spec.lone_pairs;
            if let Some(override_degree) = spec.degree_override {
                atom.degree = override_degree;
            }
            atom.steric_number = atom.degree + atom.lone_pairs;
        }
    }

    fn perceive_aromaticity(mut molecule: AnnotatedMolecule) -> AnnotatedMolecule {
        perceive(&mut molecule).expect("aromaticity perception should succeed");
        molecule
    }

    fn assert_flag_sets(
        molecule: &AnnotatedMolecule,
        expected_aromatic: &[usize],
        expected_anti: &[usize],
    ) {
        use std::collections::HashSet;
        let aromatic: HashSet<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| atom.is_aromatic.then_some(idx))
            .collect();
        let anti: HashSet<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| atom.is_anti_aromatic.then_some(idx))
            .collect();

        assert_eq!(
            aromatic,
            expected_aromatic.iter().copied().collect(),
            "unexpected aromatic atom assignment"
        );
        assert_eq!(
            anti,
            expected_anti.iter().copied().collect(),
            "unexpected anti-aromatic atom assignment"
        );

        for (idx, atom) in molecule.atoms.iter().enumerate() {
            assert!(
                !(atom.is_aromatic && atom.is_anti_aromatic),
                "atom {idx} cannot be aromatic and anti-aromatic simultaneously"
            );
        }
    }

    fn h() -> AtomSpec {
        AtomSpec::new(Element::H)
    }
    fn c() -> AtomSpec {
        AtomSpec::new(Element::C)
    }
    fn n_pyridine() -> AtomSpec {
        AtomSpec::new(Element::N).with_lone_pairs(1)
    }
    fn n_pyrrole() -> AtomSpec {
        AtomSpec::new(Element::N).with_lone_pairs(1)
    }
    fn o_furan() -> AtomSpec {
        AtomSpec::new(Element::O).with_lone_pairs(2)
    }
    fn o_carbonyl() -> AtomSpec {
        AtomSpec::new(Element::O).with_lone_pairs(2)
    }
    fn b_anion() -> AtomSpec {
        AtomSpec::new(Element::B).with_charge(-1)
    }

    fn benzene() -> AnnotatedMolecule {
        let atoms = vec![c(), c(), c(), c(), c(), c(), h(), h(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 0, BondOrder::Single),
            (0, 6, BondOrder::Single),
            (1, 7, BondOrder::Single),
            (2, 8, BondOrder::Single),
            (3, 9, BondOrder::Single),
            (4, 10, BondOrder::Single),
            (5, 11, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn pyrrole() -> AnnotatedMolecule {
        let atoms = vec![n_pyrrole(), c(), c(), c(), c(), h(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 0, BondOrder::Single),
            (0, 5, BondOrder::Single),
            (1, 6, BondOrder::Single),
            (2, 7, BondOrder::Single),
            (3, 8, BondOrder::Single),
            (4, 9, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4]])
    }

    fn borabenzene_anion() -> AnnotatedMolecule {
        let atoms = vec![b_anion(), c(), c(), c(), c(), c(), h(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 0, BondOrder::Single),
            (1, 6, BondOrder::Single),
            (2, 7, BondOrder::Single),
            (3, 8, BondOrder::Single),
            (4, 9, BondOrder::Single),
            (5, 10, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn cyclobutadiene() -> AnnotatedMolecule {
        let atoms = vec![c(), c(), c(), c(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 0, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (1, 5, BondOrder::Single),
            (2, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3]])
    }

    fn cyclooctatetraene_nonplanar() -> AnnotatedMolecule {
        let atoms: Vec<_> = (0..8)
            .map(|idx| {
                if idx % 2 == 0 {
                    c().with_degree(4)
                } else {
                    c()
                }
            })
            .chain((8..16).map(|_| h()))
            .collect();
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 6, BondOrder::Single),
            (6, 7, BondOrder::Double),
            (7, 0, BondOrder::Single),
            (0, 8, BondOrder::Single),
            (1, 9, BondOrder::Single),
            (2, 10, BondOrder::Single),
            (3, 11, BondOrder::Single),
            (4, 12, BondOrder::Single),
            (5, 13, BondOrder::Single),
            (6, 14, BondOrder::Single),
            (7, 15, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5, 6, 7]])
    }

    fn cyclohexane() -> AnnotatedMolecule {
        let atoms = (0..18)
            .map(|i| if i < 6 { c() } else { h() })
            .collect::<Vec<_>>();
        let bonds = vec![
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (5, 0, BondOrder::Single),
            (0, 6, BondOrder::Single),
            (0, 7, BondOrder::Single),
            (1, 8, BondOrder::Single),
            (1, 9, BondOrder::Single),
            (2, 10, BondOrder::Single),
            (2, 11, BondOrder::Single),
            (3, 12, BondOrder::Single),
            (3, 13, BondOrder::Single),
            (4, 14, BondOrder::Single),
            (4, 15, BondOrder::Single),
            (5, 16, BondOrder::Single),
            (5, 17, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn naphthalene() -> AnnotatedMolecule {
        let atoms = (0..18)
            .map(|i| if i < 10 { c() } else { h() })
            .collect::<Vec<_>>();
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 9, BondOrder::Single),
            (9, 8, BondOrder::Double),
            (8, 7, BondOrder::Single),
            (7, 6, BondOrder::Double),
            (6, 5, BondOrder::Single),
            (5, 0, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (0, 10, BondOrder::Single),
            (1, 11, BondOrder::Single),
            (2, 12, BondOrder::Single),
            (3, 13, BondOrder::Single),
            (6, 14, BondOrder::Single),
            (7, 15, BondOrder::Single),
            (8, 16, BondOrder::Single),
            (9, 17, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5], &[4, 5, 6, 7, 8, 9]])
    }

    fn anthracene() -> AnnotatedMolecule {
        let atoms = (0..24)
            .map(|i| if i < 14 { c() } else { h() })
            .collect::<Vec<_>>();
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 13, BondOrder::Single),
            (13, 12, BondOrder::Double),
            (12, 11, BondOrder::Single),
            (11, 10, BondOrder::Double),
            (10, 5, BondOrder::Single),
            (5, 0, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (10, 9, BondOrder::Single),
            (9, 8, BondOrder::Double),
            (8, 7, BondOrder::Single),
            (7, 6, BondOrder::Double),
            (6, 11, BondOrder::Single),
            (0, 14, BondOrder::Single),
            (1, 15, BondOrder::Single),
            (2, 16, BondOrder::Single),
            (3, 17, BondOrder::Single),
            (6, 18, BondOrder::Single),
            (7, 19, BondOrder::Single),
            (8, 20, BondOrder::Single),
            (9, 21, BondOrder::Single),
            (12, 22, BondOrder::Single),
            (13, 23, BondOrder::Single),
        ];
        build_test_molecule(
            &atoms,
            &bonds,
            &[
                &[0, 1, 2, 3, 4, 5],
                &[4, 5, 10, 11, 6, 7, 8, 9, 13, 12],
                &[6, 7, 8, 9, 10, 11],
            ],
        )
    }

    fn purine() -> AnnotatedMolecule {
        let atoms = vec![
            n_pyridine(),
            c(),
            n_pyridine(),
            c(),
            c(),
            n_pyrrole(),
            c(),
            n_pyridine(),
            c(),
            h(),
            h(),
            h(),
            h(),
        ];
        let bonds = vec![
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (5, 0, BondOrder::Single),
            (3, 8, BondOrder::Single),
            (8, 7, BondOrder::Double),
            (7, 6, BondOrder::Single),
            (6, 4, BondOrder::Double),
            (1, 9, BondOrder::Single),
            (5, 10, BondOrder::Single),
            (6, 11, BondOrder::Single),
            (8, 12, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5], &[3, 4, 6, 7, 8]])
    }

    fn alpha_pyrone() -> AnnotatedMolecule {
        let atoms = vec![
            c(),
            o_furan(),
            c(),
            c(),
            c(),
            c(),
            o_carbonyl(),
            h(),
            h(),
            h(),
            h(),
        ];
        let bonds = vec![
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 0, BondOrder::Single),
            (0, 6, BondOrder::Double),
            (2, 7, BondOrder::Single),
            (3, 8, BondOrder::Single),
            (4, 9, BondOrder::Single),
            (5, 10, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn pyrazole() -> AnnotatedMolecule {
        let atoms = vec![c(), n_pyrrole(), n_pyridine(), c(), c(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 0, BondOrder::Double),
            (0, 5, BondOrder::Single),
            (1, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
            (4, 8, BondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4]])
    }

    #[test]
    fn benzene_ring_is_aromatic() {
        let molecule = perceive_aromaticity(benzene());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5], &[]);
    }

    #[test]
    fn pyrrole_lone_pair_contributes_to_aromaticity() {
        let molecule = perceive_aromaticity(pyrrole());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4], &[]);
    }

    #[test]
    fn borabenzene_anion_is_aromatic() {
        let molecule = perceive_aromaticity(borabenzene_anion());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5], &[]);
    }

    #[test]
    fn cyclobutadiene_detected_as_antiaromatic() {
        let molecule = perceive_aromaticity(cyclobutadiene());
        assert_flag_sets(&molecule, &[], &[0, 1, 2, 3]);
    }

    #[test]
    fn cyclooctatetraene_rejected_due_to_non_planarity() {
        let molecule = perceive_aromaticity(cyclooctatetraene_nonplanar());
        assert_flag_sets(&molecule, &[], &[]);
    }

    #[test]
    fn cyclohexane_is_non_aromatic() {
        let molecule = perceive_aromaticity(cyclohexane());
        assert_flag_sets(&molecule, &[], &[]);
    }

    #[test]
    fn naphthalene_fused_rings_are_aromatic() {
        let molecule = perceive_aromaticity(naphthalene());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], &[]);
    }

    #[test]
    fn anthracene_three_ring_system_is_aromatic() {
        let molecule = perceive_aromaticity(anthracene());
        assert_flag_sets(&molecule, &(0..14).collect::<Vec<_>>(), &[]);
    }

    #[test]
    fn purine_dual_ring_system_is_aromatic() {
        let molecule = perceive_aromaticity(purine());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5, 6, 7, 8], &[]);
    }

    #[test]
    fn alpha_pyrone_is_non_aromatic() {
        let molecule = perceive_aromaticity(alpha_pyrone());
        assert_flag_sets(&molecule, &[], &[]);
    }

    #[test]
    fn pyrazole_is_aromatic() {
        let molecule = perceive_aromaticity(pyrazole());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4], &[]);
    }
}
