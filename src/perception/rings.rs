//! Detects rings and records small-set cycle representatives for subsequent perception stages.
//!
//! This module builds a minimal cycle basis from the molecular graph so aromaticity, resonance,
//! and hybridization passes can quickly determine ring membership and sizes.

use super::model::{AnnotatedMolecule, Ring};
use crate::core::error::PerceptionError;
use crate::core::properties::GraphBondOrder;
use std::collections::{HashMap, VecDeque};

/// Computes ring information for the supplied annotated molecule.
///
/// Runs connected-component counting, enumerates simple cycle candidates, chooses a minimal cycle
/// basis, and marks atoms with ring membership metadata.
///
/// # Arguments
///
/// * `molecule` - Mutable annotated molecule that will receive ring annotations.
///
/// # Returns
///
/// `Ok(())` whether rings exist or not; the step keeps the `Result` signature to match the broader
/// perception pipeline but currently never emits [`PerceptionError`].
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

    let sssr_candidates =
        select_minimal_cycle_basis(candidates, cyclomatic_number as usize, &bond_id_to_index);

    let final_rings: Vec<Ring> = sssr_candidates
        .into_iter()
        .map(|c| {
            let mut atom_ids = c.atom_ids;
            atom_ids.sort_unstable();
            atom_ids
        })
        .collect();
    molecule.rings = final_rings;

    annotate_atoms_with_ring_info(molecule);

    Ok(())
}

/// Cycle descriptor storing both atom and bond identifiers.
struct RingCandidate {
    /// Ordered atom identifiers along the candidate cycle.
    atom_ids: Vec<usize>,
    /// Ordered bond identifiers along the candidate cycle.
    bond_ids: Vec<usize>,
    /// Cycle length measured in edges.
    len: usize,
}

/// Enumerates simple cycles by removing each bond and searching for alternate paths.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose adjacency and bonds will be analyzed.
///
/// # Returns
///
/// Collection of candidate rings containing atom and bond identifiers.
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

/// Selects up to `cyclomatic_number` cycles forming a minimal basis using Gaussian elimination.
///
/// # Arguments
///
/// * `candidates` - Candidate cycles sorted by length.
/// * `cyclomatic_number` - Target number of independent cycles to keep.
/// * `bond_id_to_index` - Mapping from bond IDs to dense indices for bit-vector math.
///
/// # Returns
///
/// Minimal cycle basis expressed as a subset of the input candidates.
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

/// Marks atoms as ring members and records their smallest ring size.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule updated in-place.
fn annotate_atoms_with_ring_info(molecule: &mut AnnotatedMolecule) {
    for ring in &molecule.rings {
        let ring_size = ring.len() as u8;
        for &atom_id in ring {
            if let Some(props) = molecule.atoms.get_mut(atom_id) {
                props.is_in_ring = true;

                let current_smallest = props.smallest_ring_size.get_or_insert(ring_size);
                if ring_size < *current_smallest {
                    *current_smallest = ring_size;
                }
            }
        }
    }
}

/// Stores the path discovered between two atoms when a bond is removed.
struct PathData {
    /// Atom identifiers along the path (excluding the destination, which is implied).
    atom_ids: Vec<usize>,
    /// Bond identifiers traversed along the path.
    bond_ids: Vec<usize>,
    /// Path length measured in edges.
    len: usize,
}

/// BFS-based shortest path search that optionally excludes one bond from consideration.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule providing adjacency and bond metadata.
/// * `start_id` - Starting atom identifier.
/// * `end_id` - Destination atom identifier.
/// * `excluded_bond_id` - Optional bond ID to ignore, simulating its removal.
///
/// # Returns
///
/// A [`PathData`] instance if a route exists or `None` when the nodes disconnect.
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

/// Counts the number of connected components in the molecular graph.
///
/// # Arguments
///
/// * `num_atoms` - Number of atoms in the graph.
/// * `adjacency` - Neighbor list for each atom.
///
/// # Returns
///
/// Count of disjoint components.
fn count_components(num_atoms: usize, adjacency: &[Vec<(usize, GraphBondOrder)>]) -> usize {
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

/// Sparse bit vector used for Gaussian elimination over GF(2).
#[derive(Clone, Debug, PartialEq, Eq)]
struct BitVec {
    /// Packed bit storage.
    data: Vec<u64>,
    /// Number of significant bits represented by the vector.
    size: usize,
}

impl BitVec {
    /// Creates a zero-initialized bit vector with the requested number of bits.
    ///
    /// # Arguments
    ///
    /// * `size` - Number of significant bits to represent.
    fn new(size: usize) -> Self {
        let words = size.div_ceil(64);
        Self {
            data: vec![0; words],
            size,
        }
    }

    /// Builds a bit vector with ones at positions mapped from the provided bond identifiers.
    ///
    /// # Arguments
    ///
    /// * `bond_ids` - Bonds participating in the cycle.
    /// * `bond_id_to_index` - Dense index lookup for bonds.
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

    /// XORs the bit vector with another vector of equal length.
    ///
    /// # Arguments
    ///
    /// * `other` - Right-hand operand.
    fn xor(&mut self, other: &Self) {
        debug_assert_eq!(self.data.len(), other.data.len());
        for (a, b) in self.data.iter_mut().zip(&other.data) {
            *a ^= *b;
        }
    }

    /// Tests whether the bit at `index` is set.
    ///
    /// # Arguments
    ///
    /// * `index` - Position to inspect.
    fn test(&self, index: usize) -> bool {
        if index >= self.size {
            return false;
        }
        let word_idx = index / 64;
        let bit_idx = index % 64;
        (self.data[word_idx] & (1u64 << bit_idx)) != 0
    }

    /// Returns the index of the most significant set bit, if any.
    ///
    /// # Returns
    ///
    /// `Some(index)` when a set bit exists or `None` when the vector is zero.
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
