use super::super::{AtomBlueprint, BondBlueprint, MoleculeTestCase};
use dreid_typer::{BondOrder, Element};

pub const GLYCINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Glycine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HA2",
            element: Element::H,
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
