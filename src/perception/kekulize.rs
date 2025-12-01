//! Resolves aromatic bonds into concrete single/double assignments via a Kekulé solver.
//!
//! The logic here isolates aromatic systems, validates that they sit inside rings, and runs a
//! backtracking search that respects valence limits and heteroatom allowances before updating the
//! annotated molecule in-place.

use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::{Element, GraphBondOrder};
use std::collections::{HashMap, VecDeque, hash_map::Entry};

/// Converts aromatic bonds inside the molecule to alternating single/double assignments.
///
/// The pass validates that every aromatic bond belongs to a ring, partitions the bonds into
/// connected systems, and then runs a Kekulé solver per system.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose bond orders and adjacency lists are mutated in place.
///
/// # Returns
///
/// `Ok(())` when every aromatic system receives a valid assignment.
///
/// # Errors
///
/// Returns [`PerceptionError::KekulizationFailed`] when an aromatic bond lies outside a ring or no
/// valid alternating assignment exists for a system.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let mut aromatic_bonds = Vec::new();
    for bond in &molecule.bonds {
        if bond.order == GraphBondOrder::Aromatic {
            aromatic_bonds.push(bond.id);
        }
    }

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

/// Backtracking assignment helper that finds valid bond orders for one aromatic system.
struct KekuleSolver<'a> {
    molecule: &'a AnnotatedMolecule,
    bond_indices: Vec<usize>,
    assignments: Vec<Option<GraphBondOrder>>,
}

impl<'a> KekuleSolver<'a> {
    /// Creates a solver scoped to the provided bond identifiers.
    ///
    /// # Arguments
    ///
    /// * `molecule` - Annotated molecule providing bond/atom metadata.
    /// * `system_bond_ids` - Aromatic bond IDs belonging to one connected system.
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

    /// Attempts to assign single/double orders to every bond in the system.
    ///
    /// # Returns
    ///
    /// Map of bond IDs to resolved orders, or `None` if no assignment satisfies the constraints.
    fn solve(&mut self) -> Option<HashMap<usize, GraphBondOrder>> {
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

    /// Recursively assigns orders while pruning inconsistent branches.
    ///
    /// # Arguments
    ///
    /// * `k` - Index of the bond currently being assigned.
    ///
    /// # Returns
    ///
    /// `true` when a full assignment is found downstream.
    fn backtrack(&mut self, k: usize) -> bool {
        if k == self.assignments.len() {
            return true;
        }

        for &order_choice in &[GraphBondOrder::Double, GraphBondOrder::Single] {
            self.assignments[k] = Some(order_choice);

            if self.is_consistent(k) && self.backtrack(k + 1) {
                return true;
            }
        }

        self.assignments[k] = None;
        false
    }

    /// Checks whether the partial assignment at index `k` respects valence constraints.
    fn is_consistent(&self, k: usize) -> bool {
        let bond_idx = self.bond_indices[k];
        let (u, v) = self.molecule.bonds[bond_idx].atom_ids;

        self.is_valence_ok(u) && self.is_valence_ok(v)
    }

    /// Verifies that the atom's accumulated valence does not exceed its maximum allowed value.
    ///
    /// Lone allowances let aromatic nitrogens/phosphorus carry at most one double bond while still
    /// being counted as having valence one toward aromatic contributions.
    fn is_valence_ok(&self, atom_id: usize) -> bool {
        let max_valence = get_max_valence(self.molecule.atoms[atom_id].element);
        let mut current_valence = 0;
        let mut aromatic_double_allowance = 0u8;
        let element = self.molecule.atoms[atom_id].element;

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

            if *initial_order == GraphBondOrder::Aromatic {
                if let Some(assigned_order) = self
                    .bond_indices
                    .iter()
                    .position(|&idx| self.molecule.bonds[idx].id == bond.id)
                    .and_then(|pos| self.assignments[pos])
                {
                    let contribution = if matches!(element, Element::N | Element::P)
                        && assigned_order == GraphBondOrder::Double
                    {
                        if aromatic_double_allowance >= 1 {
                            return false;
                        }
                        aromatic_double_allowance += 1;
                        1
                    } else {
                        bond_order_to_valence(assigned_order)
                    };
                    current_valence += contribution;
                }
            } else {
                current_valence += bond_order_to_valence(*initial_order);
            }
        }

        current_valence <= max_valence
    }
}

/// Groups aromatic bonds into connected systems for independent solving.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule providing adjacency and bond metadata.
/// * `aromatic_bonds` - All bond IDs currently marked as aromatic.
///
/// # Returns
///
/// Vector of bond-ID lists, each representing one connected aromatic system.
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
                    if *order == GraphBondOrder::Aromatic {
                        let neighbor_bond = molecule
                            .bonds
                            .iter()
                            .find(|b| {
                                (b.atom_ids.0 == atom_id && b.atom_ids.1 == *neighbor_id)
                                    || (b.atom_ids.0 == *neighbor_id && b.atom_ids.1 == atom_id)
                            })
                            .unwrap();

                        if let Entry::Vacant(entry) = visited_bonds.entry(neighbor_bond.id) {
                            entry.insert(true);
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

/// Ensures every aromatic bond resides entirely inside a perceived ring.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule containing ring flags.
/// * `aromatic_bonds` - Bond IDs that must be validated.
///
/// # Errors
///
/// Returns [`PerceptionError::KekulizationFailed`] when any bond touches a non-ring atom.
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

/// Returns the maximum valence allowed for the provided element.
///
/// # Arguments
///
/// * `element` - Element whose typical valence limit should be enforced.
///
/// # Returns
///
/// Maximum valence counted toward aromatic assignments.
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

/// Converts a bond order into its valence contribution.
///
/// # Arguments
///
/// * `order` - Bond order being counted toward valence.
///
/// # Returns
///
/// Integer contribution consistent with typical valence bookkeeping.
fn bond_order_to_valence(order: GraphBondOrder) -> u8 {
    match order {
        GraphBondOrder::Single => 1,
        GraphBondOrder::Double => 2,
        GraphBondOrder::Triple => 3,
        GraphBondOrder::Aromatic => 0,
    }
}
