//! Internal model shared across perception stages to annotate atoms, bonds, rings, and adjacency.
//!
//! The structures defined here wrap the raw `MolecularGraph` with mutable fields that each
//! perception pass enriches before the typing engine consumes them.

use crate::core::error::GraphValidationError;
use crate::core::graph::{BondEdge, MolecularGraph};
use crate::core::properties::{Element, GraphBondOrder, Hybridization};

/// Perception-friendly atom record that stores both graph identity and inferred properties.
#[derive(Debug, Clone)]
pub struct AnnotatedAtom {
    /// Zero-based identifier matching the source [`MolecularGraph`].
    pub id: usize,
    /// Chemical element of the atom.
    pub element: Element,

    /// Current formal charge assigned by electron perception.
    pub formal_charge: i8,
    /// Number of lone pairs tracked for hybridization and resonance logic.
    pub lone_pairs: u8,
    /// Graph degree computed during adjacency building.
    pub degree: u8,

    /// Whether the atom lies on any ring identified so far.
    pub is_in_ring: bool,
    /// Size of the smallest ring containing the atom, if any.
    pub smallest_ring_size: Option<u8>,

    /// Flag set once aromaticity perception confirms Huckel criteria for this atom.
    pub is_aromatic: bool,
    /// Flag set for anti-aromatic atoms that should avoid resonance promotion.
    pub is_anti_aromatic: bool,

    /// Marks atoms that participate in a resonance system (aromatic or functional group).
    pub is_resonant: bool,

    /// Steric number derived from lone pairs and neighbors for VSEPR calculations.
    pub steric_number: u8,
    /// Current hybridization assignment, defaulting to [`Hybridization::Unknown`].
    pub hybridization: Hybridization,
}

/// Convenience alias representing a ring as a list of atom identifiers.
pub type Ring = Vec<usize>;

/// Represents a subset of atoms and bonds that form a delocalized electron system.
///
/// This structure is populated by the `aromaticity` (for rings) and `resonance`
/// (for linear conjugated groups like carboxylates) passes. It is used by the
/// builder to identify which bonds should be emitted with `TopologyBondOrder::Resonant`.
#[derive(Debug, Clone)]
pub struct ResonanceSystem {
    /// IDs of atoms participating in this resonance system.
    pub atom_ids: Vec<usize>,
    /// IDs of bonds (`BondEdge::id`) that are part of the delocalized system.
    pub bond_ids: Vec<usize>,
}

/// Annotated molecular container combining atom metadata, bonds, adjacency, and ring sets.
#[derive(Debug, Clone)]
pub struct AnnotatedMolecule {
    /// All atoms with perception-specific annotations.
    pub atoms: Vec<AnnotatedAtom>,
    /// Copy of the graph bonds to provide stable IDs and connectivity.
    pub bonds: Vec<BondEdge>,
    /// Adjacency list capturing neighbor IDs and bond orders.
    pub adjacency: Vec<Vec<(usize, GraphBondOrder)>>,
    /// Collection of rings discovered during perception.
    pub rings: Vec<Ring>,
    /// Collection of all identified resonance systems.
    pub resonance_systems: Vec<ResonanceSystem>,
}

impl AnnotatedMolecule {
    /// Builds an annotated molecule from a validated [`MolecularGraph`].
    ///
    /// Initializes adjacency lists for every atom, clones the bond list, and seeds default
    /// annotations that later perception passes will populate.
    ///
    /// # Arguments
    ///
    /// * `graph` - Source molecular graph supplied by the caller.
    ///
    /// # Returns
    ///
    /// A new [`AnnotatedMolecule`] containing adjacency, cloned bonds, and zeroed annotations.
    ///
    /// # Errors
    ///
    /// Returns [`GraphValidationError::MissingAtom`] if any bond endpoint references an atom index
    /// outside the graph's atom list.
    pub fn new(graph: &MolecularGraph) -> Result<Self, GraphValidationError> {
        let mut adjacency = vec![vec![]; graph.atoms.len()];
        for bond in &graph.bonds {
            let (u, v) = bond.atom_ids;

            if u >= graph.atoms.len() {
                return Err(GraphValidationError::MissingAtom { atom_id: u });
            }
            if v >= graph.atoms.len() {
                return Err(GraphValidationError::MissingAtom { atom_id: v });
            }

            adjacency[u].push((v, bond.order));
            adjacency[v].push((u, bond.order));
        }

        let atoms = graph
            .atoms
            .iter()
            .map(|node| AnnotatedAtom {
                id: node.id,
                element: node.element,
                degree: adjacency[node.id].len() as u8,
                formal_charge: 0,
                lone_pairs: 0,
                is_in_ring: false,
                smallest_ring_size: None,
                is_aromatic: false,
                is_anti_aromatic: false,
                is_resonant: false,
                steric_number: 0,
                hybridization: Hybridization::Unknown,
            })
            .collect();

        Ok(Self {
            atoms,
            bonds: graph.bonds.clone(),
            adjacency,
            rings: Vec::new(),
            resonance_systems: Vec::new(),
        })
    }
}
