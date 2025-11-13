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
