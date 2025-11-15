use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::BondOrder;
use std::collections::{HashMap, VecDeque};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let num_atoms = molecule.atoms.len();
    if num_atoms == 0 {
        return Ok(());
    }
    let num_bonds = molecule.bonds.len();
    let num_components = count_components(num_atoms, &molecule.adjacency);
    let cyclomatic_number = num_bonds as isize - num_atoms as isize + num_components as isize;

    if cyclomatic_number <= 0 {
        return Ok(());
    }

    let bond_id_to_index: HashMap<usize, usize> = molecule
        .bonds
        .iter()
        .map(|b| b.id)
        .enumerate()
        .map(|(i, id)| (id, i))
        .collect();

    let candidates = enumerate_cycle_candidates(molecule);

    let sssr =
        select_minimal_cycle_basis(candidates, cyclomatic_number as usize, &bond_id_to_index);

    annotate_molecule_with_rings(molecule, sssr);

    Ok(())
}

struct RingCandidate {
    atom_ids: Vec<usize>,
    bond_ids: Vec<usize>,
    len: usize,
}

fn enumerate_cycle_candidates(molecule: &AnnotatedMolecule) -> Vec<RingCandidate> {
    let mut candidates = Vec::new();

    for bond_to_remove in &molecule.bonds {
        if let Some(path) = shortest_path_bfs(
            molecule,
            bond_to_remove.atom_ids.0,
            bond_to_remove.atom_ids.1,
            Some(bond_to_remove.id),
        ) {
            let mut atom_ids = path.atom_ids;
            let mut bond_ids = path.bond_ids;
            atom_ids.push(bond_to_remove.atom_ids.1);
            bond_ids.push(bond_to_remove.id);

            candidates.push(RingCandidate {
                atom_ids,
                bond_ids,
                len: path.len + 1,
            });
        }
    }
    candidates
}

fn select_minimal_cycle_basis(
    mut candidates: Vec<RingCandidate>,
    cyclomatic_number: usize,
    bond_id_to_index: &HashMap<usize, usize>,
) -> Vec<RingCandidate> {
    candidates.sort_by_key(|c| c.len);

    let mut selected_rings = Vec::new();
    let mut basis: Vec<(BitVec, usize)> = Vec::new();

    for ring in candidates {
        let mut bitvec = BitVec::from_bond_ids(&ring.bond_ids, bond_id_to_index);

        for (basis_vec, pivot) in &basis {
            if bitvec.test(*pivot) {
                bitvec.xor(basis_vec);
            }
        }

        if let Some(pivot) = bitvec.leading_one() {
            basis.push((bitvec, pivot));
            basis.sort_by_key(|&(_, p)| p);
            selected_rings.push(ring);

            if selected_rings.len() == cyclomatic_number {
                break;
            }
        }
    }
    selected_rings
}

fn annotate_molecule_with_rings(molecule: &mut AnnotatedMolecule, rings: Vec<RingCandidate>) {
    for ring in rings {
        let ring_size = ring.atom_ids.len() as u8;
        for atom_id in ring.atom_ids {
            let props = &mut molecule.atoms[atom_id];
            props.is_in_ring = true;
            match props.smallest_ring_size {
                Some(current_size) if ring_size < current_size => {
                    props.smallest_ring_size = Some(ring_size);
                }
                None => {
                    props.smallest_ring_size = Some(ring_size);
                }
                _ => {}
            }
        }
    }
}

struct PathData {
    atom_ids: Vec<usize>,
    bond_ids: Vec<usize>,
    len: usize,
}

fn shortest_path_bfs(
    molecule: &AnnotatedMolecule,
    start_id: usize,
    end_id: usize,
    excluded_bond_id: Option<usize>,
) -> Option<PathData> {
    let mut queue = VecDeque::new();
    let mut visited = vec![false; molecule.atoms.len()];
    let mut parent: Vec<Option<(usize, usize)>> = vec![None; molecule.atoms.len()];

    visited[start_id] = true;
    queue.push_back(start_id);

    'outer: while let Some(current_id) = queue.pop_front() {
        for &(neighbor_id, _bond_order) in &molecule.adjacency[current_id] {
            let bond = molecule
                .bonds
                .iter()
                .find(|b| {
                    (b.atom_ids.0 == current_id && b.atom_ids.1 == neighbor_id)
                        || (b.atom_ids.0 == neighbor_id && b.atom_ids.1 == current_id)
                })
                .unwrap();

            if Some(bond.id) == excluded_bond_id {
                continue;
            }
            if !visited[neighbor_id] {
                visited[neighbor_id] = true;
                parent[neighbor_id] = Some((current_id, bond.id));
                queue.push_back(neighbor_id);

                if neighbor_id == end_id {
                    break 'outer;
                }
            }
        }
    }

    if !visited[end_id] {
        return None;
    }

    let mut atom_ids = Vec::new();
    let mut bond_ids = Vec::new();
    let mut len = 0;
    let mut cursor = end_id;

    while let Some((prev_id, bond_id)) = parent[cursor] {
        atom_ids.push(prev_id);
        bond_ids.push(bond_id);
        len += 1;
        cursor = prev_id;
        if cursor == start_id {
            break;
        }
    }
    atom_ids.reverse();
    bond_ids.reverse();

    Some(PathData {
        atom_ids,
        bond_ids,
        len,
    })
}

fn count_components(num_atoms: usize, adjacency: &[Vec<(usize, BondOrder)>]) -> usize {
    let mut visited = vec![false; num_atoms];
    let mut components = 0;
    for i in 0..num_atoms {
        if !visited[i] {
            components += 1;
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(current) = stack.pop() {
                for &(neighbor_id, _) in &adjacency[current] {
                    if !visited[neighbor_id] {
                        visited[neighbor_id] = true;
                        stack.push(neighbor_id);
                    }
                }
            }
        }
    }
    components
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct BitVec {
    data: Vec<u64>,
    size: usize,
}

impl BitVec {
    fn new(size: usize) -> Self {
        let words = size.div_ceil(64);
        Self {
            data: vec![0; words],
            size,
        }
    }

    fn from_bond_ids(bond_ids: &[usize], bond_id_to_index: &HashMap<usize, usize>) -> Self {
        let mut bitvec = Self::new(bond_id_to_index.len());
        for bond_id in bond_ids {
            if let Some(&index) = bond_id_to_index.get(bond_id) {
                let word_idx = index / 64;
                let bit_idx = index % 64;
                if word_idx < bitvec.data.len() {
                    bitvec.data[word_idx] |= 1u64 << bit_idx;
                }
            }
        }
        bitvec
    }

    fn xor(&mut self, other: &Self) {
        debug_assert_eq!(self.data.len(), other.data.len());
        for (a, b) in self.data.iter_mut().zip(&other.data) {
            *a ^= *b;
        }
    }

    fn test(&self, index: usize) -> bool {
        if index >= self.size {
            return false;
        }
        let word_idx = index / 64;
        let bit_idx = index % 64;
        (self.data[word_idx] & (1u64 << bit_idx)) != 0
    }

    fn leading_one(&self) -> Option<usize> {
        for (word_idx_rev, &word) in self.data.iter().enumerate().rev() {
            if word != 0 {
                let bit_pos_in_word = 63 - word.leading_zeros() as usize;
                return Some(word_idx_rev * 64 + bit_pos_in_word);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::properties::{BondOrder, Element};

    fn chain_graph(len: usize) -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        for _ in 0..len {
            graph.add_atom(Element::C);
        }
        for i in 0..len.saturating_sub(1) {
            graph
                .add_bond(i, i + 1, BondOrder::Single)
                .expect("valid bond");
        }
        graph
    }

    fn cycle_graph(len: usize) -> MolecularGraph {
        assert!(len >= 3, "cycles require at least three atoms");
        let mut graph = MolecularGraph::new();
        for _ in 0..len {
            graph.add_atom(Element::C);
        }
        for i in 0..len {
            let next = (i + 1) % len;
            graph
                .add_bond(i, next, BondOrder::Single)
                .expect("valid cycle bond");
        }
        graph
    }

    #[test]
    fn perceive_skips_acyclic_molecules() {
        let chain = chain_graph(4);
        let mut molecule = AnnotatedMolecule::new(&chain).expect("graph is valid");

        perceive(&mut molecule).expect("perception should succeed");

        assert!(
            molecule
                .atoms
                .iter()
                .all(|atom| !atom.is_in_ring && atom.smallest_ring_size.is_none())
        );
    }

    #[test]
    fn perceive_marks_ring_atoms_with_smallest_ring_size() {
        let square = cycle_graph(4);
        let mut molecule = AnnotatedMolecule::new(&square).expect("graph is valid");

        perceive(&mut molecule).expect("perception should succeed");

        for atom in &molecule.atoms {
            assert!(atom.is_in_ring, "atom {} should be in ring", atom.id);
            assert_eq!(atom.smallest_ring_size, Some(4));
        }
    }

    #[test]
    fn shortest_path_bfs_finds_alternative_route_when_edge_removed() {
        let triangle = cycle_graph(3);
        let molecule = AnnotatedMolecule::new(&triangle).expect("graph is valid");

        let removed_bond_id = molecule
            .bonds
            .iter()
            .find(|bond| bond.atom_ids == (0, 1) || bond.atom_ids == (1, 0))
            .map(|bond| bond.id)
            .expect("triangle should contain 0-1 bond");

        let path = shortest_path_bfs(&molecule, 0, 1, Some(removed_bond_id))
            .expect("path exists through third atom");

        assert_eq!(path.len, 2);
        assert_eq!(path.atom_ids, vec![0, 2]);
        assert_eq!(path.bond_ids.len(), 2);
    }

    #[test]
    fn count_components_detects_disconnected_fragments() {
        let adjacency = vec![
            vec![(1, BondOrder::Single)],
            vec![(0, BondOrder::Single)],
            vec![(3, BondOrder::Single)],
            vec![(2, BondOrder::Single)],
        ];

        assert_eq!(count_components(4, &adjacency), 2);
    }

    #[test]
    fn bitvec_supports_xor_and_leading_one() {
        let mut bond_map = HashMap::new();
        bond_map.insert(10usize, 0usize);
        bond_map.insert(20usize, 1usize);
        bond_map.insert(30usize, 2usize);

        let mut a = BitVec::from_bond_ids(&[10, 30], &bond_map);
        let b = BitVec::from_bond_ids(&[20, 30], &bond_map);

        assert!(a.test(0));
        assert!(a.test(2));
        assert!(b.test(1));
        assert!(b.test(2));

        a.xor(&b);
        assert!(a.test(0));
        assert!(a.test(1));
        assert!(!a.test(2));

        assert_eq!(a.leading_one(), Some(1));
    }
}
