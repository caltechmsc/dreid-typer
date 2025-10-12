use super::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AtomNode {
    pub id: usize,
    pub element: Element,
    pub formal_charge: i8,
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

    pub fn add_atom(&mut self, element: Element, formal_charge: i8) -> usize {
        let id = self.atoms.len();
        self.atoms.push(AtomNode {
            id,
            element,
            formal_charge,
        });
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

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MolecularTopology {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub proper_dihedrals: Vec<ProperDihedral>,
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
        let atom_id_1 = graph.add_atom(Element::C);
        assert_eq!(atom_id_1, 0);
        assert_eq!(graph.atoms.len(), 1);
        assert_eq!(graph.atoms[0].id, 0);
        assert_eq!(graph.atoms[0].element, Element::C);

        let atom_id_2 = graph.add_atom(Element::H);
        assert_eq!(atom_id_2, 1);
        assert_eq!(graph.atoms.len(), 2);
        assert_eq!(graph.atoms[1].id, 1);
        assert_eq!(graph.atoms[1].element, Element::H);
    }

    #[test]
    fn molecular_graph_add_bond_succeeds() {
        let mut graph = MolecularGraph::new();
        graph.add_atom(Element::C);
        graph.add_atom(Element::C);
        let bond_id = graph.add_bond(0, 1, BondOrder::Single).unwrap();
        assert_eq!(bond_id, 0);
        assert_eq!(graph.bonds.len(), 1);
        assert_eq!(graph.bonds[0].atom_ids, (0, 1));
        assert_eq!(graph.bonds[0].order, BondOrder::Single);
    }

    #[test]
    fn molecular_graph_add_bond_with_out_of_bounds_atom_id() {
        let mut graph = MolecularGraph::new();
        graph.add_atom(Element::C);
        let result = graph.add_bond(0, 1, BondOrder::Single);
        assert!(result.is_err());
    }

    #[test]
    fn molecular_graph_add_bond_to_itself() {
        let mut graph = MolecularGraph::new();
        graph.add_atom(Element::C);
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
