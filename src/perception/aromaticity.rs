//! Evaluates fused ring systems to determine whether the atoms are aromatic, anti-aromatic, or neither.
//!
//! The module groups rings, builds localized models that count π-electrons under planarity
//! assumptions, and sets per-atom flags. Importantly, it also registers aromatic rings
//! as `ResonanceSystem`s so that their bonds are treated as resonant in the final topology.

use super::model::{AnnotatedAtom, AnnotatedMolecule, ResonanceSystem, Ring};
use crate::core::error::PerceptionError;
use crate::core::properties::GraphBondOrder;
use std::collections::{HashMap, HashSet};

/// Runs aromaticity perception over all ring systems present in the molecule.
///
/// The procedure clusters rings that share atoms, evaluates each cluster as a whole, falls back to
/// ring-by-ring evaluation when mixed behavior occurs, and annotates atoms as aromatic or
/// anti-aromatic accordingly. Confirmed aromatic systems are added to the molecule's
/// resonance systems list.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose atom flags should be updated.
///
/// # Returns
///
/// `Ok(())` once every ring system has been processed or when no rings exist.
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
            apply_aromaticity(molecule, &system_atoms);
        } else if model.is_anti_aromatic() {
            for &atom_id in &system_atoms {
                molecule.atoms[atom_id].is_anti_aromatic = true;
            }
        } else {
            evaluate_rings_individually(molecule, &system_indices);
        }
    }

    Ok(())
}

/// Marks all atoms and bonds in the system as aromatic and resonant.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to mutate.
/// * `system_atoms` - Atom IDs representing the aromatic system.
fn apply_aromaticity(molecule: &mut AnnotatedMolecule, system_atoms: &HashSet<usize>) {
    for &atom_id in system_atoms {
        let atom = &mut molecule.atoms[atom_id];
        atom.is_aromatic = true;
        atom.is_resonant = true;
    }

    let mut bond_ids = Vec::new();
    for bond in &molecule.bonds {
        if system_atoms.contains(&bond.atom_ids.0) && system_atoms.contains(&bond.atom_ids.1) {
            bond_ids.push(bond.id);
        }
    }

    let mut atom_ids: Vec<usize> = system_atoms.iter().copied().collect();
    atom_ids.sort_unstable();
    bond_ids.sort_unstable();

    molecule
        .resonance_systems
        .push(ResonanceSystem { atom_ids, bond_ids });
}

/// Evaluates each ring independently when a fused system lacks uniform behavior.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to mutate.
/// * `system_indices` - Indices of rings belonging to the fused system.
fn evaluate_rings_individually(molecule: &mut AnnotatedMolecule, system_indices: &[usize]) {
    for &ring_idx in system_indices {
        let ring_atoms: HashSet<_> = molecule.rings[ring_idx].iter().copied().collect();
        let ring_model = AromaticityModel::new(molecule, &ring_atoms);

        if ring_model.is_aromatic() {
            apply_aromaticity(molecule, &ring_atoms);
        } else if ring_model.is_anti_aromatic() {
            for &atom_id in &ring_atoms {
                molecule.atoms[atom_id].is_anti_aromatic = true;
            }
        }
    }
}

/// Local model capturing the atoms and π-electron count for a ring system.
struct AromaticityModel<'a> {
    /// Annotated molecule providing adjacency information.
    molecule: &'a AnnotatedMolecule,
    /// Atom IDs forming the current system under evaluation.
    atoms: HashSet<usize>,
    /// Computed π-electron count, if evaluation succeeded.
    pi_electrons: Option<u32>,
    /// Flag describing whether the atoms satisfy the planarity heuristic.
    is_potentially_planar: bool,
}

impl<'a> AromaticityModel<'a> {
    /// Constructs the model and immediately evaluates planarity and π-electrons.
    ///
    /// # Arguments
    ///
    /// * `molecule` - Annotated molecule backing the model.
    /// * `system_atoms` - Atom IDs representing a ring system.
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

    /// Returns `true` when the Huckel 4n+2 rule is satisfied.
    fn is_aromatic(&self) -> bool {
        if !self.is_potentially_planar {
            return false;
        }
        matches!(self.pi_electrons, Some(pi) if pi > 0 && (pi - 2) % 4 == 0)
    }

    /// Returns `true` when the Huckel 4n rule indicates anti-aromaticity.
    fn is_anti_aromatic(&self) -> bool {
        if !self.is_potentially_planar || self.has_cross_conjugation() {
            return false;
        }
        matches!(self.pi_electrons, Some(pi) if pi > 0 && pi % 4 == 0)
    }

    /// Computes planarity and π-electron counts, caching the results on the struct.
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

    /// Checks if an atom participates in a double bond within the ring system.
    fn atom_has_endocyclic_double(&self, atom_id: usize) -> bool {
        self.molecule.adjacency[atom_id]
            .iter()
            .any(|&(n_id, order)| order == GraphBondOrder::Double && self.atoms.contains(&n_id))
    }

    /// Checks if an atom carries a double bond outside the ring system.
    fn atom_has_exocyclic_double(&self, atom_id: usize) -> bool {
        self.molecule.adjacency[atom_id]
            .iter()
            .any(|&(n_id, order)| order == GraphBondOrder::Double && !self.atoms.contains(&n_id))
    }

    /// Detects cross-conjugation, which prevents anti-aromatic classification.
    fn has_cross_conjugation(&self) -> bool {
        self.atoms
            .iter()
            .any(|&atom_id| self.atom_has_exocyclic_double(atom_id))
    }

    /// Computes each atom's π contribution using bond, lone-pair, and resonance flags.
    fn count_pi_contribution(&self, atom_id: usize) -> Option<u32> {
        let atom = &self.molecule.atoms[atom_id];
        let has_endocyclic_double_bond = self.atom_has_endocyclic_double(atom_id);
        let has_exocyclic_double_bond = self.atom_has_exocyclic_double(atom_id);

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

        if atom.is_resonant && atom.is_in_ring {
            return Some(1);
        }

        None
    }
}

/// Heuristic planarity test derived from steric number rules.
fn is_potentially_planar(atom: &AnnotatedAtom) -> bool {
    let steric_number = atom.degree + atom.lone_pairs;
    match steric_number {
        0..=3 => true,
        4 => atom.lone_pairs > 0,
        _ => false,
    }
}

/// Groups rings into fused systems via shared atoms.
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

/// Builds an adjacency list between rings that share at least one atom.
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
