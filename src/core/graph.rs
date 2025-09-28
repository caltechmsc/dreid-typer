use super::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AtomNode {
    pub id: usize,
    pub element: Element,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BondEdge {
    pub id: usize,
    pub atom_ids: (usize, usize),
    pub order: BondOrder,
}

#[derive(Debug, Clone, PartialEq, Default)]
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

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub id: usize,
    pub element: Element,
    pub atom_type: String,
    pub hybridization: Hybridization,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    pub atom_ids: (usize, usize),
    pub order: BondOrder,
}

impl Bond {
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Angle {
    pub atom_ids: (usize, usize, usize),
}

impl Angle {
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ProperDihedral {
    pub atom_ids: (usize, usize, usize, usize),
}

impl ProperDihedral {
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ImproperDihedral {
    pub atom_ids: (usize, usize, usize, usize),
}

impl ImproperDihedral {
    pub fn new(plane_id1: usize, plane_id2: usize, center_id: usize, plane_id3: usize) -> Self {
        let mut plane_ids = [plane_id1, plane_id2, plane_id3];
        plane_ids.sort_unstable();
        Self {
            atom_ids: (plane_ids[0], plane_ids[1], center_id, plane_ids[2]),
        }
    }
}
