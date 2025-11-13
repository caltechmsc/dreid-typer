use super::error::GraphValidationError;
use super::properties::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone)]
pub struct AtomNode {
    pub id: usize,
    pub element: Element,
}

#[derive(Debug, Clone)]
pub struct BondEdge {
    pub id: usize,
    pub atom_ids: (usize, usize),
    pub order: BondOrder,
}

#[derive(Debug, Clone, Default)]
pub struct MolecularGraph {
    pub atoms: Vec<AtomNode>,
    pub bonds: Vec<BondEdge>,
}

impl MolecularGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_atom(&mut self, element: Element) -> usize {
        let id = self.atoms.len();
        self.atoms.push(AtomNode { id, element });
        id
    }

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
