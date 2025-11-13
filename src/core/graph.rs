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
