use crate::core::error::GraphValidationError;
use crate::core::graph::{BondEdge, MolecularGraph};
use crate::core::properties::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone)]
pub struct AnnotatedAtom {
    pub id: usize,
    pub element: Element,

    pub formal_charge: i8,
    pub lone_pairs: u8,
    pub degree: u8,

    pub is_in_ring: bool,
    pub smallest_ring_size: Option<u8>,

    pub is_aromatic: bool,
    pub is_anti_aromatic: bool,

    pub is_in_conjugated_system: bool,
    pub is_resonant: bool,
    pub steric_number: u8,
    pub hybridization: Hybridization,
}

pub type Ring = Vec<usize>;

#[derive(Debug, Clone)]
pub struct AnnotatedMolecule {
    pub atoms: Vec<AnnotatedAtom>,
    pub bonds: Vec<BondEdge>,
    pub adjacency: Vec<Vec<(usize, BondOrder)>>,
    pub rings: Vec<Ring>,
}

impl AnnotatedMolecule {
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
                is_in_conjugated_system: false,
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
        })
    }
}

use pauling::{
    AtomId as PaulingId, BondOrder as PaulingOrder, Element as PaulingElement, traits::AtomView,
    traits::BondView, traits::MoleculeGraph,
};

impl From<Element> for PaulingElement {
    fn from(element: Element) -> Self {
        match element {
            Element::H => PaulingElement::H,
            Element::He => PaulingElement::He,
            Element::Li => PaulingElement::Li,
            Element::Be => PaulingElement::Be,
            Element::B => PaulingElement::B,
            Element::C => PaulingElement::C,
            Element::N => PaulingElement::N,
            Element::O => PaulingElement::O,
            Element::F => PaulingElement::F,
            Element::Ne => PaulingElement::Ne,
            Element::Na => PaulingElement::Na,
            Element::Mg => PaulingElement::Mg,
            Element::Al => PaulingElement::Al,
            Element::Si => PaulingElement::Si,
            Element::P => PaulingElement::P,
            Element::S => PaulingElement::S,
            Element::Cl => PaulingElement::Cl,
            Element::Ar => PaulingElement::Ar,
            Element::K => PaulingElement::K,
            Element::Ca => PaulingElement::Ca,
            Element::Sc => PaulingElement::Sc,
            Element::Ti => PaulingElement::Ti,
            Element::V => PaulingElement::V,
            Element::Cr => PaulingElement::Cr,
            Element::Mn => PaulingElement::Mn,
            Element::Fe => PaulingElement::Fe,
            Element::Co => PaulingElement::Co,
            Element::Ni => PaulingElement::Ni,
            Element::Cu => PaulingElement::Cu,
            Element::Zn => PaulingElement::Zn,
            Element::Ga => PaulingElement::Ga,
            Element::Ge => PaulingElement::Ge,
            Element::As => PaulingElement::As,
            Element::Se => PaulingElement::Se,
            Element::Br => PaulingElement::Br,
            Element::Kr => PaulingElement::Kr,
            Element::Rb => PaulingElement::Rb,
            Element::Sr => PaulingElement::Sr,
            Element::Y => PaulingElement::Y,
            Element::Zr => PaulingElement::Zr,
            Element::Nb => PaulingElement::Nb,
            Element::Mo => PaulingElement::Mo,
            Element::Tc => PaulingElement::Tc,
            Element::Ru => PaulingElement::Ru,
            Element::Rh => PaulingElement::Rh,
            Element::Pd => PaulingElement::Pd,
            Element::Ag => PaulingElement::Ag,
            Element::Cd => PaulingElement::Cd,
            Element::In => PaulingElement::In,
            Element::Sn => PaulingElement::Sn,
            Element::Sb => PaulingElement::Sb,
            Element::Te => PaulingElement::Te,
            Element::I => PaulingElement::I,
            Element::Xe => PaulingElement::Xe,
            Element::Cs => PaulingElement::Cs,
            Element::Ba => PaulingElement::Ba,
            Element::La => PaulingElement::La,
            Element::Ce => PaulingElement::Ce,
            Element::Pr => PaulingElement::Pr,
            Element::Nd => PaulingElement::Nd,
            Element::Pm => PaulingElement::Pm,
            Element::Sm => PaulingElement::Sm,
            Element::Eu => PaulingElement::Eu,
            Element::Gd => PaulingElement::Gd,
            Element::Tb => PaulingElement::Tb,
            Element::Dy => PaulingElement::Dy,
            Element::Ho => PaulingElement::Ho,
            Element::Er => PaulingElement::Er,
            Element::Tm => PaulingElement::Tm,
            Element::Yb => PaulingElement::Yb,
            Element::Lu => PaulingElement::Lu,
            Element::Hf => PaulingElement::Hf,
            Element::Ta => PaulingElement::Ta,
            Element::W => PaulingElement::W,
            Element::Re => PaulingElement::Re,
            Element::Os => PaulingElement::Os,
            Element::Ir => PaulingElement::Ir,
            Element::Pt => PaulingElement::Pt,
            Element::Au => PaulingElement::Au,
            Element::Hg => PaulingElement::Hg,
            Element::Tl => PaulingElement::Tl,
            Element::Pb => PaulingElement::Pb,
            Element::Bi => PaulingElement::Bi,
            Element::Po => PaulingElement::Po,
            Element::At => PaulingElement::At,
            Element::Rn => PaulingElement::Rn,
            Element::Fr => PaulingElement::Fr,
            Element::Ra => PaulingElement::Ra,
            Element::Ac => PaulingElement::Ac,
            Element::Th => PaulingElement::Th,
            Element::Pa => PaulingElement::Pa,
            Element::U => PaulingElement::U,
            Element::Np => PaulingElement::Np,
            Element::Pu => PaulingElement::Pu,
            Element::Am => PaulingElement::Am,
            Element::Cm => PaulingElement::Cm,
            Element::Bk => PaulingElement::Bk,
            Element::Cf => PaulingElement::Cf,
            Element::Es => PaulingElement::Es,
            Element::Fm => PaulingElement::Fm,
            Element::Md => PaulingElement::Md,
            Element::No => PaulingElement::No,
            Element::Lr => PaulingElement::Lr,
            Element::Rf => PaulingElement::Rf,
            Element::Db => PaulingElement::Db,
            Element::Sg => PaulingElement::Sg,
            Element::Bh => PaulingElement::Bh,
            Element::Hs => PaulingElement::Hs,
            Element::Mt => PaulingElement::Mt,
            Element::Ds => PaulingElement::Ds,
            Element::Rg => PaulingElement::Rg,
            Element::Cn => PaulingElement::Cn,
            Element::Nh => PaulingElement::Nh,
            Element::Fl => PaulingElement::Fl,
            Element::Mc => PaulingElement::Mc,
            Element::Lv => PaulingElement::Lv,
            Element::Ts => PaulingElement::Ts,
            Element::Og => PaulingElement::Og,
        }
    }
}

impl From<BondOrder> for PaulingOrder {
    fn from(order: BondOrder) -> Self {
        match order {
            BondOrder::Single => PaulingOrder::Single,
            BondOrder::Double => PaulingOrder::Double,
            BondOrder::Triple => PaulingOrder::Triple,
            BondOrder::Aromatic => PaulingOrder::Aromatic,
        }
    }
}

impl AtomView for AnnotatedAtom {
    fn id(&self) -> PaulingId {
        self.id
    }

    fn element(&self) -> PaulingElement {
        self.element.into()
    }

    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
}

impl BondView for BondEdge {
    fn id(&self) -> PaulingId {
        self.id
    }

    fn order(&self) -> PaulingOrder {
        self.order.into()
    }

    fn start_atom_id(&self) -> PaulingId {
        self.atom_ids.0
    }

    fn end_atom_id(&self) -> PaulingId {
        self.atom_ids.1
    }
}

impl MoleculeGraph for AnnotatedMolecule {
    type Atom = AnnotatedAtom;
    type Bond = BondEdge;

    fn atoms(&self) -> impl Iterator<Item = &Self::Atom> {
        self.atoms.iter()
    }

    fn bonds(&self) -> impl Iterator<Item = &Self::Bond> {
        self.bonds.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::AtomNode;
    use pauling::traits::{AtomView, BondView, MoleculeGraph};

    fn water_like_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let oxygen = graph.add_atom(Element::O);
        let hydrogen1 = graph.add_atom(Element::H);
        let hydrogen2 = graph.add_atom(Element::H);

        graph
            .add_bond(oxygen, hydrogen1, BondOrder::Single)
            .expect("valid bond should be created");
        graph
            .add_bond(oxygen, hydrogen2, BondOrder::Single)
            .expect("valid bond should be created");

        graph
    }

    #[test]
    fn annotated_molecule_new_populates_adjacency_and_atoms() {
        let graph = water_like_graph();
        let molecule = AnnotatedMolecule::new(&graph).expect("graph should be valid");

        let adjacency_sizes: Vec<_> = molecule
            .adjacency
            .iter()
            .map(|neighbors| neighbors.len())
            .collect();
        assert_eq!(adjacency_sizes, vec![2, 1, 1]);

        let oxygen_neighbors = &molecule.adjacency[0];
        assert!(oxygen_neighbors.contains(&(1, BondOrder::Single)));
        assert!(oxygen_neighbors.contains(&(2, BondOrder::Single)));

        let oxygen = &molecule.atoms[0];
        assert_eq!(oxygen.element, Element::O);
        assert_eq!(oxygen.degree, 2);
        assert_eq!(oxygen.formal_charge, 0);
        assert_eq!(oxygen.lone_pairs, 0);
        assert!(!oxygen.is_in_ring);
        assert_eq!(oxygen.smallest_ring_size, None);
        assert!(!oxygen.is_aromatic);
        assert!(!oxygen.is_anti_aromatic);
        assert!(!oxygen.is_in_conjugated_system);
        assert!(!oxygen.is_resonant);
        assert_eq!(oxygen.steric_number, 0);
        assert_eq!(oxygen.hybridization, Hybridization::Unknown);
    }

    #[test]
    fn annotated_molecule_new_detects_invalid_bond_endpoints() {
        let graph = MolecularGraph {
            atoms: vec![AtomNode {
                id: 0,
                element: Element::C,
            }],
            bonds: vec![BondEdge {
                id: 0,
                atom_ids: (0, 2),
                order: BondOrder::Single,
            }],
        };

        let err = AnnotatedMolecule::new(&graph).expect_err("invalid bond must fail");

        match err {
            GraphValidationError::MissingAtom { atom_id } => assert_eq!(atom_id, 2),
            other => panic!("unexpected error returned: {other:?}"),
        }
    }

    #[test]
    fn annotated_atom_satisfies_atom_view_contract() {
        let atom = AnnotatedAtom {
            id: 7,
            element: Element::N,
            formal_charge: -1,
            lone_pairs: 2,
            degree: 3,
            is_in_ring: true,
            smallest_ring_size: Some(5),
            is_aromatic: true,
            is_anti_aromatic: false,
            is_in_conjugated_system: true,
            is_resonant: true,
            steric_number: 3,
            hybridization: Hybridization::SP2,
        };

        assert_eq!(atom.id(), 7);
        assert_eq!(atom.element(), pauling::Element::N);
        assert_eq!(atom.formal_charge(), -1);
    }

    #[test]
    fn bond_edge_satisfies_bond_view_contract() {
        let bond = BondEdge {
            id: 3,
            atom_ids: (9, 4),
            order: BondOrder::Aromatic,
        };

        assert_eq!(bond.id(), 3);
        assert_eq!(bond.start_atom_id(), 9);
        assert_eq!(bond.end_atom_id(), 4);
        assert_eq!(bond.order(), pauling::BondOrder::Aromatic);
    }

    #[test]
    fn annotated_molecule_molecule_graph_iterators_match_source_graph() {
        let graph = water_like_graph();
        let molecule = AnnotatedMolecule::new(&graph).unwrap();

        let atom_ids: Vec<_> = molecule.atoms().map(|atom| atom.id()).collect();
        assert_eq!(atom_ids, vec![0, 1, 2]);

        let bond_pairs: Vec<_> = molecule
            .bonds()
            .map(|bond| (bond.id(), bond.start_atom_id(), bond.end_atom_id()))
            .collect();
        assert_eq!(bond_pairs.len(), 2);
        assert!(bond_pairs.contains(&(0, 0, 1)));
        assert!(bond_pairs.contains(&(1, 0, 2)));
    }
}
