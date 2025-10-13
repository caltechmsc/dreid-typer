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

pub const GLYCINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Glycine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            charge: 1,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            charge: 0,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            charge: 0,
            expected_type: "C_2",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            charge: -1,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            charge: 0,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            charge: 0,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            charge: 0,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            charge: 0,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA1",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HA2",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
    ],
    bonds: &[
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: BondOrder::Double,
        },
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA2",
            order: BondOrder::Single,
        },
    ],
};

pub const ALANINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Alanine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            charge: 1,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            charge: 0,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            charge: 0,
            expected_type: "C_2",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            charge: -1,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            charge: 0,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "CB",
            element: Element::C,
            charge: 0,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            charge: 0,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            charge: 0,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            charge: 0,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB3",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
    ],
    bonds: &[
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: BondOrder::Double,
        },
        BondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB3",
            order: BondOrder::Single,
        },
    ],
};
