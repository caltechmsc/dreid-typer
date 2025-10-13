use dreid_typer::{BondOrder, Element};

#[derive(Debug)]
pub struct AtomBlueprint {
    pub label: &'static str,
    pub element: Element,
    pub charge: i8,
    pub expected_type: &'static str,
}

#[derive(Debug)]
pub struct BondBlueprint {
    pub atom1_label: &'static str,
    pub atom2_label: &'static str,
    pub order: BondOrder,
}

#[derive(Debug)]
pub struct MoleculeTestCase {
    pub name: &'static str,
    pub atoms: &'static [AtomBlueprint],
    pub bonds: &'static [BondBlueprint],
}
