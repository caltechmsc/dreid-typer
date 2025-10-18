//! Core data structures for representing molecular graphs and topologies.
//!
//! This module defines the fundamental types used throughout the dreid-typer library,
//! including the input `MolecularGraph` for chemical connectivity and the output
//! `MolecularTopology` with assigned atom types and perceived structural elements.
//! These structures form the basis of the three-phase pipeline: perception, typing, and building.

use super::{BondOrder, Element, Hybridization};

/// Represents an atom in a molecular graph with its basic properties.
///
/// This struct is used as a node in the `MolecularGraph` to store essential
/// information about each atom, including its unique identifier and chemical element.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AtomNode {
    /// The unique identifier for this atom within the graph.
    pub id: usize,
    /// The chemical element of this atom.
    pub element: Element,
}

/// Represents a bond between two atoms in a molecular graph.
///
/// This struct is used as an edge in the `MolecularGraph` to define the connectivity
/// and bond order between atoms.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BondEdge {
    /// The unique identifier for this bond within the graph.
    pub id: usize,
    /// The identifiers of the two atoms connected by this bond, stored as a tuple.
    pub atom_ids: (usize, usize),
    /// The order of this bond (single, double, triple, etc.).
    pub order: BondOrder,
}

/// A simple representation of a molecule's connectivity as a graph.
///
/// This struct serves as the input to the dreid-typer pipeline, containing only
/// the basic chemical connectivity information: atoms with their elements and
/// formal charges, and bonds with their orders. It does not include any derived
/// properties like hybridization or atom types.
///
/// # Examples
///
/// Creating a simple methane molecule:
///
/// ```
/// use dreid_typer::{MolecularGraph, Element, BondOrder};
///
/// let mut graph = MolecularGraph::new();
/// let carbon = graph.add_atom(Element::C, 0);
/// let hydrogen1 = graph.add_atom(Element::H, 0);
/// let hydrogen2 = graph.add_atom(Element::H, 0);
/// let hydrogen3 = graph.add_atom(Element::H, 0);
/// let hydrogen4 = graph.add_atom(Element::H, 0);
///
/// graph.add_bond(carbon, hydrogen1, BondOrder::Single).unwrap();
/// graph.add_bond(carbon, hydrogen2, BondOrder::Single).unwrap();
/// graph.add_bond(carbon, hydrogen3, BondOrder::Single).unwrap();
/// graph.add_bond(carbon, hydrogen4, BondOrder::Single).unwrap();
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MolecularGraph {
    /// The list of atoms in the molecule.
    pub atoms: Vec<AtomNode>,
    /// The list of bonds connecting the atoms.
    pub bonds: Vec<BondEdge>,
}

impl MolecularGraph {
    /// Creates a new, empty molecular graph.
    ///
    /// # Returns
    ///
    /// A new `MolecularGraph` instance with no atoms or bonds.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::MolecularGraph;
    ///
    /// let graph = MolecularGraph::new();
    /// assert!(graph.atoms.is_empty());
    /// assert!(graph.bonds.is_empty());
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a new atom to the molecular graph.
    ///
    /// The atom is assigned a unique ID based on the current number of atoms,
    /// and its properties are stored in the graph.
    ///
    /// # Arguments
    ///
    /// * `element` - The chemical element of the atom.
    /// * `formal_charge` - The formal charge of the atom.
    ///
    /// # Returns
    ///
    /// The unique ID assigned to the newly added atom.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{MolecularGraph, Element};
    ///
    /// let mut graph = MolecularGraph::new();
    /// let atom_id = graph.add_atom(Element::C, 0);
    /// assert_eq!(atom_id, 0);
    /// assert_eq!(graph.atoms.len(), 1);
    /// ```
    pub fn add_atom(&mut self, element: Element, formal_charge: i8) -> usize {
        let id = self.atoms.len();
        self.atoms.push(AtomNode {
            id,
            element,
            formal_charge,
        });
        id
    }

    /// Adds a new bond between two atoms in the molecular graph.
    ///
    /// The bond is validated to ensure both atom IDs are valid and not the same.
    /// If successful, the bond is assigned a unique ID and added to the graph.
    ///
    /// # Arguments
    ///
    /// * `atom1_id` - The ID of the first atom in the bond.
    /// * `atom2_id` - The ID of the second atom in the bond.
    /// * `order` - The bond order (single, double, etc.).
    ///
    /// # Returns
    ///
    /// A `Result` containing the unique ID of the newly added bond on success,
    /// or an error message as a string slice on failure.
    ///
    /// # Errors
    ///
    /// Returns an error if either atom ID is out of bounds or if the atoms are the same.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{MolecularGraph, Element, BondOrder};
    ///
    /// let mut graph = MolecularGraph::new();
    /// let atom1 = graph.add_atom(Element::C, 0);
    /// let atom2 = graph.add_atom(Element::C, 0);
    /// let bond_id = graph.add_bond(atom1, atom2, BondOrder::Single).unwrap();
    /// assert_eq!(bond_id, 0);
    /// ```
    pub fn add_bond(
        &mut self,
        atom1_id: usize,
        atom2_id: usize,
        order: BondOrder,
    ) -> Result<usize, &'static str> {
        if atom1_id >= self.atoms.len() || atom2_id >= self.atoms.len() {
            return Err("Cannot add bond: atom ID is out of bounds");
        }
        if atom1_id == atom2_id {
            return Err("Cannot add bond: an atom cannot bond to itself");
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

/// Represents an atom in the final molecular topology with assigned properties.
///
/// This struct extends the basic atom information with the assigned DREIDING
/// atom type and hybridization state, as determined by the typing phase.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// The unique identifier for this atom.
    pub id: usize,
    /// The chemical element of this atom.
    pub element: Element,
    /// The assigned DREIDING force field atom type.
    pub atom_type: String,
    /// The hybridization state of this atom.
    pub hybridization: Hybridization,
}

/// Represents a bond in the molecular topology.
///
/// This struct defines a bond between two atoms with its order, used in the
/// final topology output.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    /// The identifiers of the two atoms connected by this bond, sorted in ascending order.
    pub atom_ids: (usize, usize),
    /// The order of this bond.
    pub order: BondOrder,
}

impl Bond {
    /// Creates a new bond between two atoms.
    ///
    /// The atom IDs are automatically sorted to ensure consistent representation.
    ///
    /// # Arguments
    ///
    /// * `id1` - The ID of the first atom.
    /// * `id2` - The ID of the second atom.
    /// * `order` - The bond order.
    ///
    /// # Returns
    ///
    /// A new `Bond` instance with sorted atom IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{Bond, BondOrder};
    ///
    /// let bond = Bond::new(2, 1, BondOrder::Single);
    /// assert_eq!(bond.atom_ids, (1, 2));
    /// ```
    pub fn new(id1: usize, id2: usize, order: BondOrder) -> Self {
        if id1 < id2 {
            Self {
                atom_ids: (id1, id2),
                order,
            }
        } else {
            Self {
                atom_ids: (id2, id1),
                order,
            }
        }
    }
}

/// Represents an angle formed by three atoms in the molecular topology.
///
/// An angle is defined by two bonds sharing a common central atom.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Angle {
    /// The identifiers of the three atoms forming the angle, with the central atom in the middle.
    pub atom_ids: (usize, usize, usize),
}

impl Angle {
    /// Creates a new angle between three atoms.
    ///
    /// The outer atom IDs are sorted to ensure consistent representation.
    ///
    /// # Arguments
    ///
    /// * `id1` - The ID of the first outer atom.
    /// * `center_id` - The ID of the central atom.
    /// * `id2` - The ID of the second outer atom.
    ///
    /// # Returns
    ///
    /// A new `Angle` instance with sorted outer atom IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::Angle;
    ///
    /// let angle = Angle::new(3, 1, 2);
    /// assert_eq!(angle.atom_ids, (2, 1, 3));
    /// ```
    pub fn new(id1: usize, center_id: usize, id2: usize) -> Self {
        if id1 < id2 {
            Self {
                atom_ids: (id1, center_id, id2),
            }
        } else {
            Self {
                atom_ids: (id2, center_id, id1),
            }
        }
    }
}

/// Represents a proper dihedral angle in the molecular topology.
///
/// A proper dihedral is defined by four atoms in a chain, measuring the
/// torsion angle around the central bond.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ProperDihedral {
    /// The identifiers of the four atoms defining the dihedral, in order.
    pub atom_ids: (usize, usize, usize, usize),
}

impl ProperDihedral {
    /// Creates a new proper dihedral from four atoms.
    ///
    /// The atom sequence is chosen to be lexicographically smallest between
    /// the forward and reverse representations.
    ///
    /// # Arguments
    ///
    /// * `id1` - The ID of the first atom.
    /// * `id2` - The ID of the second atom.
    /// * `id3` - The ID of the third atom.
    /// * `id4` - The ID of the fourth atom.
    ///
    /// # Returns
    ///
    /// A new `ProperDihedral` instance.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::ProperDihedral;
    ///
    /// let dihedral = ProperDihedral::new(4, 3, 2, 1);
    /// assert_eq!(dihedral.atom_ids, (1, 2, 3, 4));
    /// ```
    pub fn new(id1: usize, id2: usize, id3: usize, id4: usize) -> Self {
        let forward = (id1, id2, id3, id4);
        let reverse = (id4, id3, id2, id1);
        if forward <= reverse {
            Self { atom_ids: forward }
        } else {
            Self { atom_ids: reverse }
        }
    }
}

/// Represents an improper dihedral angle in the molecular topology.
///
/// An improper dihedral measures the out-of-plane angle for atoms like those
/// in carbonyl groups or aromatic rings.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ImproperDihedral {
    /// The identifiers of the four atoms defining the improper dihedral.
    pub atom_ids: (usize, usize, usize, usize),
}

impl ImproperDihedral {
    /// Creates a new improper dihedral from four atoms.
    ///
    /// The three plane atoms are sorted to ensure consistent representation.
    ///
    /// # Arguments
    ///
    /// * `plane_id1` - The ID of the first plane atom.
    /// * `plane_id2` - The ID of the second plane atom.
    /// * `center_id` - The ID of the central atom.
    /// * `plane_id3` - The ID of the third plane atom.
    ///
    /// # Returns
    ///
    /// A new `ImproperDihedral` instance with sorted plane atoms.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::ImproperDihedral;
    ///
    /// let dihedral = ImproperDihedral::new(4, 2, 1, 3);
    /// assert_eq!(dihedral.atom_ids, (2, 3, 1, 4));
    /// ```
    pub fn new(plane_id1: usize, plane_id2: usize, center_id: usize, plane_id3: usize) -> Self {
        let mut plane_ids = [plane_id1, plane_id2, plane_id3];
        plane_ids.sort_unstable();
        Self {
            atom_ids: (plane_ids[0], plane_ids[1], center_id, plane_ids[2]),
        }
    }
}

/// The complete molecular topology with assigned atom types and structural elements.
///
/// This struct represents the final output of the dreid-typer pipeline,
/// containing all atoms with their DREIDING types, plus the perceived bonds,
/// angles, and dihedrals needed for force field calculations.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MolecularTopology {
    /// The list of atoms with their assigned types and properties.
    pub atoms: Vec<Atom>,
    /// The list of bonds in the molecule.
    pub bonds: Vec<Bond>,
    /// The list of angles formed by bonded atom triplets.
    pub angles: Vec<Angle>,
    /// The list of proper dihedral angles.
    pub proper_dihedrals: Vec<ProperDihedral>,
    /// The list of improper dihedral angles.
    pub improper_dihedrals: Vec<ImproperDihedral>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{BondOrder, Element};

    #[test]
    fn molecular_graph_new_is_empty() {
        let graph = MolecularGraph::new();
        assert!(graph.atoms.is_empty());
        assert!(graph.bonds.is_empty());
    }

    #[test]
    fn molecular_graph_add_atom() {
        let mut graph = MolecularGraph::new();
        let atom_id_1 = graph.add_atom(Element::C, 0);
        assert_eq!(atom_id_1, 0);
        assert_eq!(graph.atoms.len(), 1);
        assert_eq!(graph.atoms[0].id, 0);
        assert_eq!(graph.atoms[0].element, Element::C);
        assert_eq!(graph.atoms[0].formal_charge, 0);

        let atom_id_2 = graph.add_atom(Element::H, 0);
        assert_eq!(atom_id_2, 1);
        assert_eq!(graph.atoms.len(), 2);
        assert_eq!(graph.atoms[1].id, 1);
        assert_eq!(graph.atoms[1].element, Element::H);
        assert_eq!(graph.atoms[1].formal_charge, 0);
    }

    #[test]
    fn molecular_graph_add_bond_succeeds() {
        let mut graph = MolecularGraph::new();
        graph.add_atom(Element::C, 0);
        graph.add_atom(Element::C, 0);
        let bond_id = graph.add_bond(0, 1, BondOrder::Single).unwrap();
        assert_eq!(bond_id, 0);
        assert_eq!(graph.bonds.len(), 1);
        assert_eq!(graph.bonds[0].atom_ids, (0, 1));
        assert_eq!(graph.bonds[0].order, BondOrder::Single);
    }

    #[test]
    fn molecular_graph_add_bond_with_out_of_bounds_atom_id() {
        let mut graph = MolecularGraph::new();
        graph.add_atom(Element::C, 0);
        let result = graph.add_bond(0, 1, BondOrder::Single);
        assert!(result.is_err());
    }

    #[test]
    fn molecular_graph_add_bond_to_itself() {
        let mut graph = MolecularGraph::new();
        graph.add_atom(Element::C, 0);
        let result = graph.add_bond(0, 0, BondOrder::Single);
        assert!(result.is_err());
    }

    #[test]
    fn bond_new_sorts_atom_ids() {
        let bond = Bond::new(2, 1, BondOrder::Single);
        assert_eq!(bond.atom_ids, (1, 2));
    }

    #[test]
    fn bond_new_with_pre_sorted_atom_ids() {
        let bond = Bond::new(1, 2, BondOrder::Single);
        assert_eq!(bond.atom_ids, (1, 2));
    }

    #[test]
    fn angle_new_sorts_outer_atom_ids() {
        let angle = Angle::new(3, 1, 2);
        assert_eq!(angle.atom_ids, (2, 1, 3));
    }

    #[test]
    fn angle_new_with_pre_sorted_outer_atom_ids() {
        let angle = Angle::new(2, 1, 3);
        assert_eq!(angle.atom_ids, (2, 1, 3));
    }

    #[test]
    fn proper_dihedral_new_chooses_lexicographically_smallest_representation() {
        let dihedral = ProperDihedral::new(4, 3, 2, 1);
        assert_eq!(dihedral.atom_ids, (1, 2, 3, 4));
    }

    #[test]
    fn proper_dihedral_new_with_lexicographically_smallest_representation() {
        let dihedral = ProperDihedral::new(1, 2, 3, 4);
        assert_eq!(dihedral.atom_ids, (1, 2, 3, 4));
    }

    #[test]
    fn improper_dihedral_new_sorts_plane_ids() {
        let dihedral = ImproperDihedral::new(4, 2, 1, 3);
        assert_eq!(dihedral.atom_ids, (2, 3, 1, 4));
    }

    #[test]
    fn improper_dihedral_new_with_sorted_plane_ids() {
        let dihedral = ImproperDihedral::new(1, 2, 3, 4);
        assert_eq!(dihedral.atom_ids, (1, 2, 3, 4));
    }
}
