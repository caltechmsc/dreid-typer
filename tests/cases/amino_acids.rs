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

pub const VALINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Valine Zwitterion",
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
            label: "CG1",
            element: Element::C,
            charge: 0,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG2",
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
            label: "HB",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG11",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG12",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG13",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG21",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG22",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG23",
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
            atom1_label: "CB",
            atom2_label: "CG1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
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
            atom2_label: "HB",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG11",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG12",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG13",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: BondOrder::Single,
        },
    ],
};

pub const LEUCINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Leucine Zwitterion",
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
            label: "CG",
            element: Element::C,
            charge: 0,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD1",
            element: Element::C,
            charge: 0,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD2",
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
            label: "HG",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD11",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD12",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD13",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD21",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD22",
            element: Element::H,
            charge: 0,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD23",
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
            atom1_label: "CB",
            atom2_label: "CG",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD2",
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
            atom1_label: "CG",
            atom2_label: "HG",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD11",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD12",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD13",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD21",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD22",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD23",
            order: BondOrder::Single,
        },
    ],
};
