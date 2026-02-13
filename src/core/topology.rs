//! Represents the canonical topology emitted after typing completes.
//!
//! This module defines the output structures (`MolecularTopology`, `Bond`, etc.)
//! that are ready for serialization or consumption by force field engines.
//! These types use `TopologyBondOrder`, which includes physical properties like
//! resonance, unlike the input graph.

use super::properties::{Element, Hybridization, TopologyBondOrder};

/// Canonical topology produced after the typer assigns atom types and torsions.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MolecularTopology {
    /// A list of all atoms with their final assigned properties.
    pub atoms: Vec<Atom>,
    /// A list of all bonds.
    pub bonds: Vec<Bond>,
    /// A list of all three-atom angles.
    pub angles: Vec<Angle>,
    /// A list of all four-atom torsions around rotatable bonds.
    pub torsions: Vec<Torsion>,
    /// A list of all four-atom inversions for planar centers.
    pub inversions: Vec<Inversion>,
}

/// Atom entry emitted in the final topology, combining identity and typing.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// The unique identifier of the atom.
    pub id: usize,
    /// The chemical element.
    pub element: Element,
    /// The final, assigned DREIDING atom type string.
    pub atom_type: String,
    /// The perceived hybridization state.
    pub hybridization: Hybridization,
}

/// Bond entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    /// The IDs of the two atoms, sorted to ensure a canonical representation.
    pub atom_ids: (usize, usize),
    /// The order of the bond.
    pub order: TopologyBondOrder,
}

impl Bond {
    /// Creates a new bond with atom IDs sorted to a canonical order.
    pub fn new(id1: usize, id2: usize, order: TopologyBondOrder) -> Self {
        let atom_ids = if id1 < id2 { (id1, id2) } else { (id2, id1) };
        Self { atom_ids, order }
    }
}

/// Angle entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Angle {
    /// The IDs of the three atoms (`end1`, `center`, `end2`), with end atoms sorted.
    pub atom_ids: (usize, usize, usize),
}

impl Angle {
    /// Creates a new angle with end atoms sorted to a canonical order.
    pub fn new(id1: usize, center_id: usize, id2: usize) -> Self {
        let atom_ids = if id1 < id2 {
            (id1, center_id, id2)
        } else {
            (id2, center_id, id1)
        };
        Self { atom_ids }
    }
}

/// Torsion entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Torsion {
    /// The IDs of the four atoms (`i`, `j`, `k`, `l`) where `j-k` is the rotatable bond,
    /// with `i-l` sorted.
    pub atom_ids: (usize, usize, usize, usize),
}

impl Torsion {
    /// Creates a new torsion with terminal atoms sorted to a canonical order.
    pub fn new(i: usize, j: usize, k: usize, l: usize) -> Self {
        let fwd = (i, j, k, l);
        let rev = (l, k, j, i);
        let atom_ids = if fwd <= rev { fwd } else { rev };
        Self { atom_ids }
    }
}

/// Inversion entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Inversion {
    /// The IDs of the four atoms (`center`, `axis`, `plane1`, `plane2`)
    /// where `center` is the inversion center, `axis` is the unique neighbor
    /// defining the axis, with `plane1` and `plane2` sorted.
    pub atom_ids: (usize, usize, usize, usize),
}

impl Inversion {
    /// Creates a new inversion with plane atoms sorted to a canonical order.
    pub fn new(center: usize, axis: usize, plane1: usize, plane2: usize) -> Self {
        let (p1, p2) = if plane1 < plane2 {
            (plane1, plane2)
        } else {
            (plane2, plane1)
        };
        Self {
            atom_ids: (center, axis, p1, p2),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::properties::TopologyBondOrder;

    #[test]
    fn bond_new_sorts_atom_ids() {
        let bond = Bond::new(4, 1, TopologyBondOrder::Triple);
        assert_eq!(bond.atom_ids, (1, 4));
        assert_eq!(bond.order, TopologyBondOrder::Triple);
    }

    #[test]
    fn angle_new_orders_terminal_atoms() {
        let angle = Angle::new(7, 3, 2);
        assert_eq!(angle.atom_ids, (2, 3, 7));
    }

    #[test]
    fn torsion_new_canonicalizes_orientation() {
        let forward = Torsion::new(1, 2, 3, 4);
        let reversed = Torsion::new(4, 3, 2, 1);

        assert_eq!(forward.atom_ids, reversed.atom_ids);
        assert_eq!(forward.atom_ids, (1, 2, 3, 4));
    }

    #[test]
    fn inversion_new_sorts_only_plane_atoms() {
        let inv = Inversion::new(5, 9, 4, 1);
        assert_eq!(inv.atom_ids, (5, 9, 1, 4));

        let inv2 = Inversion::new(5, 1, 9, 4);
        assert_eq!(inv2.atom_ids, (5, 1, 4, 9));
        assert_ne!(inv.atom_ids, inv2.atom_ids);
    }

    #[test]
    fn inversion_three_terms_per_center_are_distinct() {
        let inv1 = Inversion::new(0, 1, 2, 3);
        let inv2 = Inversion::new(0, 2, 1, 3);
        let inv3 = Inversion::new(0, 3, 1, 2);

        assert_eq!(inv1.atom_ids, (0, 1, 2, 3));
        assert_eq!(inv2.atom_ids, (0, 2, 1, 3));
        assert_eq!(inv3.atom_ids, (0, 3, 1, 2));

        assert_ne!(inv1, inv2);
        assert_ne!(inv2, inv3);
        assert_ne!(inv1, inv3);
    }

    #[test]
    fn molecular_topology_default_is_empty() {
        let topology = MolecularTopology::default();

        assert!(topology.atoms.is_empty());
        assert!(topology.bonds.is_empty());
        assert!(topology.angles.is_empty());
        assert!(topology.torsions.is_empty());
        assert!(topology.inversions.is_empty());
    }
}
