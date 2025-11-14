use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::{BondOrder, Element};
use std::collections::{HashMap, VecDeque};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let aromatic_bonds: Vec<usize> = molecule
        .bonds
        .iter()
        .filter(|b| b.order == BondOrder::Aromatic)
        .map(|b| b.id)
        .collect();

    if aromatic_bonds.is_empty() {
        return Ok(());
    }

    validate_aromatic_bonds_in_rings(molecule, &aromatic_bonds)?;

    let aromatic_systems = find_aromatic_systems(molecule, &aromatic_bonds);

    let mut new_bond_orders = HashMap::new();

    for system_bonds in aromatic_systems {
        let mut solver = KekuleSolver::new(molecule, &system_bonds);
        match solver.solve() {
            Some(solution) => {
                new_bond_orders.extend(solution);
            }
            None => {
                return Err(PerceptionError::KekulizationFailed {
                    message: "could not find a valid Kekulé structure for an aromatic system"
                        .to_string(),
                });
            }
        }
    }

    for (bond_id, new_order) in new_bond_orders {
        let bond = molecule.bonds.iter_mut().find(|b| b.id == bond_id).unwrap();
        let (u, v) = bond.atom_ids;
        bond.order = new_order;

        for (neighbor_id, order) in molecule.adjacency[u].iter_mut() {
            if *neighbor_id == v {
                *order = new_order;
                break;
            }
        }
        for (neighbor_id, order) in molecule.adjacency[v].iter_mut() {
            if *neighbor_id == u {
                *order = new_order;
                break;
            }
        }
    }

    Ok(())
}

struct KekuleSolver<'a> {
    molecule: &'a AnnotatedMolecule,
    bond_indices: Vec<usize>,
    assignments: Vec<Option<BondOrder>>,
}

impl<'a> KekuleSolver<'a> {
    fn new(molecule: &'a AnnotatedMolecule, system_bond_ids: &[usize]) -> Self {
        let bond_indices: Vec<usize> = system_bond_ids
            .iter()
            .map(|id| molecule.bonds.iter().position(|b| b.id == *id).unwrap())
            .collect();

        Self {
            molecule,
            bond_indices: bond_indices.clone(),
            assignments: vec![None; bond_indices.len()],
        }
    }

    fn solve(&mut self) -> Option<HashMap<usize, BondOrder>> {
        if self.backtrack(0) {
            let solution = self
                .bond_indices
                .iter()
                .enumerate()
                .map(|(i, &bond_idx)| {
                    let bond_id = self.molecule.bonds[bond_idx].id;
                    let order = self.assignments[i].unwrap();
                    (bond_id, order)
                })
                .collect();
            Some(solution)
        } else {
            None
        }
    }

    fn backtrack(&mut self, k: usize) -> bool {
        if k == self.assignments.len() {
            return true;
        }

        for &order_choice in &[BondOrder::Double, BondOrder::Single] {
            self.assignments[k] = Some(order_choice);

            if self.is_consistent(k) {
                if self.backtrack(k + 1) {
                    return true;
                }
            }
        }

        self.assignments[k] = None;
        false
    }

    fn is_consistent(&self, k: usize) -> bool {
        let bond_idx = self.bond_indices[k];
        let (u, v) = self.molecule.bonds[bond_idx].atom_ids;

        self.is_valence_ok(u) && self.is_valence_ok(v)
    }

    fn is_valence_ok(&self, atom_id: usize) -> bool {
        let max_valence = get_max_valence(self.molecule.atoms[atom_id].element);
        let mut current_valence = 0;

        for (neighbor_id, initial_order) in &self.molecule.adjacency[atom_id] {
            let bond = self
                .molecule
                .bonds
                .iter()
                .find(|b| {
                    (b.atom_ids.0 == atom_id && b.atom_ids.1 == *neighbor_id)
                        || (b.atom_ids.0 == *neighbor_id && b.atom_ids.1 == atom_id)
                })
                .unwrap();

            if *initial_order == BondOrder::Aromatic {
                if let Some(pos) = self
                    .bond_indices
                    .iter()
                    .position(|&idx| self.molecule.bonds[idx].id == bond.id)
                {
                    if let Some(assigned_order) = self.assignments[pos] {
                        current_valence += bond_order_to_valence(assigned_order);
                    }
                }
            } else {
                current_valence += bond_order_to_valence(*initial_order);
            }
        }

        current_valence <= max_valence
    }
}

fn find_aromatic_systems(
    molecule: &AnnotatedMolecule,
    aromatic_bonds: &[usize],
) -> Vec<Vec<usize>> {
    let mut systems = Vec::new();
    let mut visited_bonds = HashMap::new();

    for &start_bond_id in aromatic_bonds {
        if visited_bonds.contains_key(&start_bond_id) {
            continue;
        }

        let mut current_system = Vec::new();
        let mut queue = VecDeque::new();

        queue.push_back(start_bond_id);
        visited_bonds.insert(start_bond_id, true);

        while let Some(bond_id) = queue.pop_front() {
            current_system.push(bond_id);
            let bond = molecule.bonds.iter().find(|b| b.id == bond_id).unwrap();
            let (u, v) = bond.atom_ids;

            for atom_id in [u, v] {
                for (neighbor_id, order) in &molecule.adjacency[atom_id] {
                    if *order == BondOrder::Aromatic {
                        let neighbor_bond = molecule
                            .bonds
                            .iter()
                            .find(|b| {
                                (b.atom_ids.0 == atom_id && b.atom_ids.1 == *neighbor_id)
                                    || (b.atom_ids.0 == *neighbor_id && b.atom_ids.1 == atom_id)
                            })
                            .unwrap();

                        if !visited_bonds.contains_key(&neighbor_bond.id) {
                            visited_bonds.insert(neighbor_bond.id, true);
                            queue.push_back(neighbor_bond.id);
                        }
                    }
                }
            }
        }
        systems.push(current_system);
    }
    systems
}

fn validate_aromatic_bonds_in_rings(
    molecule: &AnnotatedMolecule,
    aromatic_bonds: &[usize],
) -> Result<(), PerceptionError> {
    for &bond_id in aromatic_bonds {
        let bond = molecule.bonds.iter().find(|b| b.id == bond_id).unwrap();
        let (u, v) = bond.atom_ids;
        if !molecule.atoms[u].is_in_ring || !molecule.atoms[v].is_in_ring {
            return Err(PerceptionError::KekulizationFailed {
                message: format!(
                    "aromatic bond (ID {}) found with at least one atom not in a ring",
                    bond_id
                ),
            });
        }
    }
    Ok(())
}

fn get_max_valence(element: Element) -> u8 {
    match element {
        Element::H | Element::F | Element::Cl | Element::Br | Element::I => 1,
        Element::O | Element::S => 2,
        Element::N | Element::P => 3,
        Element::C | Element::Si => 4,
        Element::B => 3,
        _ => 8,
    }
}

fn bond_order_to_valence(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 0,
    }
}
