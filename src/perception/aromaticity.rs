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
