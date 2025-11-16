//! Represents the molecular connectivity primitives shared across perception and
//! typing, plus the canonical topology emitted after typing completes.
//!
//! The types exposed here allow callers to construct raw graphs, run the
//! perception pipeline, and inspect the resulting atoms, bonds, angles, and
//! torsions in a consistent, serializable format.

use super::error::GraphValidationError;
use super::properties::{BondOrder, Element, Hybridization};

/// Stores the identifier and element for a single atom within a
/// [`MolecularGraph`].
#[derive(Debug, Clone)]
pub struct AtomNode {
    /// Zero-based identifier assigned when the atom is inserted into the graph.
    pub id: usize,
    /// Chemical element represented by this node.
    pub element: Element,
}

/// Captures a bond between two atoms inside a [`MolecularGraph`].
#[derive(Debug, Clone)]
pub struct BondEdge {
    /// Zero-based identifier assigned sequentially as bonds are inserted.
    pub id: usize,
    /// Tuple of atom IDs that this bond connects.
    pub atom_ids: (usize, usize),
    /// Bond multiplicity recorded for the edge.
    pub order: BondOrder,
}

/// Mutable graph of atoms and bonds supplied to the perception pipeline.
#[derive(Debug, Clone, Default)]
pub struct MolecularGraph {
    /// Collection of all atoms currently present in the graph.
    pub atoms: Vec<AtomNode>,
    /// Collection of all bonds currently present in the graph.
    pub bonds: Vec<BondEdge>,
}

impl MolecularGraph {
    /// Creates an empty graph ready for atom and bond insertion.
    ///
    /// # Returns
    ///
    /// A `MolecularGraph` with zero atoms and bonds.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::MolecularGraph;
    /// let graph = MolecularGraph::new();
    /// assert!(graph.atoms.is_empty());
    /// assert!(graph.bonds.is_empty());
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a new atom with the provided [`Element`] and returns its ID.
    ///
    /// # Arguments
    ///
    /// * `element` - Chemical element to assign to the node.
    ///
    /// # Returns
    ///
    /// The zero-based identifier for the newly inserted atom.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{Element, MolecularGraph};
    /// let mut graph = MolecularGraph::new();
    /// let carbon_id = graph.add_atom(Element::C);
    /// assert_eq!(carbon_id, 0);
    /// ```
    pub fn add_atom(&mut self, element: Element) -> usize {
        let id = self.atoms.len();
        self.atoms.push(AtomNode { id, element });
        id
    }

    /// Adds a bond between two existing atoms.
    ///
    /// # Arguments
    ///
    /// * `atom1_id` - Identifier of the first atom.
    /// * `atom2_id` - Identifier of the second atom.
    /// * `order` - Bond multiplicity to record.
    ///
    /// # Returns
    ///
    /// The zero-based identifier assigned to the new bond.
    ///
    /// # Errors
    ///
    /// Returns [`GraphValidationError::MissingAtom`] if either atom ID has not
    /// been inserted, or [`GraphValidationError::SelfBondingAtom`] if both IDs
    /// refer to the same atom.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{BondOrder, Element, MolecularGraph};
    /// let mut graph = MolecularGraph::new();
    /// let c = graph.add_atom(Element::C);
    /// let o = graph.add_atom(Element::O);
    /// let bond_id = graph.add_bond(c, o, BondOrder::Double).unwrap();
    /// assert_eq!(bond_id, 0);
    /// ```
    pub fn add_bond(
        &mut self,
        atom1_id: usize,
        atom2_id: usize,
        order: BondOrder,
    ) -> Result<usize, GraphValidationError> {
        if atom1_id >= self.atoms.len() {
            return Err(GraphValidationError::MissingAtom { atom_id: atom1_id });
        }
        if atom2_id >= self.atoms.len() {
            return Err(GraphValidationError::MissingAtom { atom_id: atom2_id });
        }
        if atom1_id == atom2_id {
            return Err(GraphValidationError::SelfBondingAtom { atom_id: atom1_id });
        }

        let id = self.bonds.len();
        self.bonds.push(BondEdge {
            id,
            atom_ids: (atom1_id, atom2_id),
            order,
        });
        Ok(id)
    }
}

/// Canonical topology produced after the typer assigns atom types and torsions.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MolecularTopology {
    /// A list of all atoms with their final assigned properties.
    pub atoms: Vec<Atom>,
    /// A list of all bonds.
    pub bonds: Vec<Bond>,
    /// A list of all three-atom angles.
    pub angles: Vec<Angle>,
    /// A list of all proper dihedral angles (torsions).
    pub propers: Vec<ProperDihedral>,
    /// A list of all improper dihedral angles (out-of-plane bends).
    pub impropers: Vec<ImproperDihedral>,
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    /// The IDs of the two atoms, sorted to ensure a canonical representation.
    pub atom_ids: (usize, usize),
    /// The order of the bond.
    pub order: BondOrder,
}

impl Bond {
    pub fn new(id1: usize, id2: usize, order: BondOrder) -> Self {
        let atom_ids = if id1 < id2 { (id1, id2) } else { (id2, id1) };
        Self { atom_ids, order }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Angle {
    /// The IDs of the three atoms (`end1`, `center`, `end2`), with end atoms sorted.
    pub atom_ids: (usize, usize, usize),
}

impl Angle {
    pub fn new(id1: usize, center_id: usize, id2: usize) -> Self {
        let atom_ids = if id1 < id2 {
            (id1, center_id, id2)
        } else {
            (id2, center_id, id1)
        };
        Self { atom_ids }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ProperDihedral {
    /// The IDs of the four atoms (`a-b-c-d`), sorted lexicographically.
    pub atom_ids: (usize, usize, usize, usize),
}

impl ProperDihedral {
    pub fn new(a: usize, b: usize, c: usize, d: usize) -> Self {
        let fwd = (a, b, c, d);
        let rev = (d, c, b, a);
        let atom_ids = if fwd <= rev { fwd } else { rev };
        Self { atom_ids }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ImproperDihedral {
    /// The IDs of the four atoms (`plane1`, `plane2`, `center`, `plane3`),
    /// with plane atoms sorted.
    pub atom_ids: (usize, usize, usize, usize),
}

impl ImproperDihedral {
    pub fn new(p1: usize, p2: usize, center: usize, p3: usize) -> Self {
        let mut plane_ids = [p1, p2, p3];
        plane_ids.sort_unstable();
        let atom_ids = (plane_ids[0], plane_ids[1], center, plane_ids[2]);
        Self { atom_ids }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::error::GraphValidationError;
    use crate::core::properties::{BondOrder, Element};

    fn graph_with_atoms(elements: &[Element]) -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        for &element in elements {
            graph.add_atom(element);
        }
        graph
    }

    #[test]
    fn molecular_graph_add_atom_assigns_sequential_ids() {
        let mut graph = MolecularGraph::new();

        let carbon_id = graph.add_atom(Element::C);
        let oxygen_id = graph.add_atom(Element::O);
        let nitrogen_id = graph.add_atom(Element::N);

        assert_eq!(carbon_id, 0);
        assert_eq!(oxygen_id, 1);
        assert_eq!(nitrogen_id, 2);
        assert_eq!(graph.atoms.len(), 3);
        assert_eq!(graph.atoms[0].element, Element::C);
        assert_eq!(graph.atoms[1].element, Element::O);
        assert_eq!(graph.atoms[2].element, Element::N);
    }

    #[test]
    fn molecular_graph_add_bond_registers_edge() {
        let mut graph = graph_with_atoms(&[Element::C, Element::O]);

        let bond_id = graph
            .add_bond(0, 1, BondOrder::Double)
            .expect("adding a valid bond should succeed");

        assert_eq!(bond_id, 0);
        assert_eq!(graph.bonds.len(), 1);
        assert_eq!(graph.bonds[0].atom_ids, (0, 1));
        assert_eq!(graph.bonds[0].order, BondOrder::Double);
    }

    #[test]
    fn molecular_graph_rejects_bond_with_missing_atom() {
        let mut graph = graph_with_atoms(&[Element::C]);

        let err = graph
            .add_bond(0, 1, BondOrder::Single)
            .expect_err("bonding to a missing atom should fail");

        match err {
            GraphValidationError::MissingAtom { atom_id } => assert_eq!(atom_id, 1),
            _ => panic!("unexpected error returned: {err:?}"),
        }
    }

    #[test]
    fn molecular_graph_rejects_self_bonding() {
        let mut graph = graph_with_atoms(&[Element::C]);

        let err = graph
            .add_bond(0, 0, BondOrder::Single)
            .expect_err("self bonds should be rejected");

        match err {
            GraphValidationError::SelfBondingAtom { atom_id } => assert_eq!(atom_id, 0),
            _ => panic!("unexpected error returned: {err:?}"),
        }
    }

    #[test]
    fn bond_new_sorts_atom_ids() {
        let bond = Bond::new(4, 1, BondOrder::Triple);
        assert_eq!(bond.atom_ids, (1, 4));
        assert_eq!(bond.order, BondOrder::Triple);
    }

    #[test]
    fn angle_new_orders_terminal_atoms() {
        let angle = Angle::new(7, 3, 2);
        assert_eq!(angle.atom_ids, (2, 3, 7));
    }

    #[test]
    fn proper_dihedral_new_canonicalizes_orientation() {
        let forward = ProperDihedral::new(1, 2, 3, 4);
        let reversed = ProperDihedral::new(4, 3, 2, 1);

        assert_eq!(forward.atom_ids, reversed.atom_ids);
        assert_eq!(forward.atom_ids, (1, 2, 3, 4));
    }

    #[test]
    fn improper_dihedral_new_sorts_plane_atoms() {
        let improper = ImproperDihedral::new(9, 1, 5, 4);
        assert_eq!(improper.atom_ids, (1, 4, 5, 9));
    }

    #[test]
    fn molecular_topology_default_is_empty() {
        let topology = MolecularTopology::default();

        assert!(topology.atoms.is_empty());
        assert!(topology.bonds.is_empty());
        assert!(topology.angles.is_empty());
        assert!(topology.propers.is_empty());
        assert!(topology.impropers.is_empty());
    }
}
