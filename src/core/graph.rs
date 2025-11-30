//! Defines the molecular graph structure used for building molecules.
//!
//! This module provides the `MolecularGraph` container, which serves as the
//! primary input interface for the library. It captures atoms and bonds with
//! `GraphBondOrder` connectivity (Single, Double, Triple, Aromatic) before
//! perception begins.

use super::error::GraphValidationError;
use super::properties::{Element, GraphBondOrder};

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
    pub order: GraphBondOrder,
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
    /// use dreid_typer::{GraphBondOrder, Element, MolecularGraph};
    /// let mut graph = MolecularGraph::new();
    /// let c = graph.add_atom(Element::C);
    /// let o = graph.add_atom(Element::O);
    /// let bond_id = graph.add_bond(c, o, GraphBondOrder::Double).unwrap();
    /// assert_eq!(bond_id, 0);
    /// ```
    pub fn add_bond(
        &mut self,
        atom1_id: usize,
        atom2_id: usize,
        order: GraphBondOrder,
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
