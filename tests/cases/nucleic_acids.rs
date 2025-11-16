use super::super::{AtomBlueprint, BondBlueprint, MoleculeTestCase};
use dreid_typer::{BondOrder, Element};

pub const URACIL: MoleculeTestCase = MoleculeTestCase {
    name: "Uracil",
    atoms: &[
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N3",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "O4",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H5",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H6",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        BondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "N3",
            atom2_label: "C4",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C6",
            atom2_label: "N1",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C2",
            atom2_label: "O2",
            order: BondOrder::Double,
        },
        BondBlueprint {
            atom1_label: "C4",
            atom2_label: "O4",
            order: BondOrder::Double,
        },
        BondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N3",
            atom2_label: "H3",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C5",
            atom2_label: "H5",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C6",
            atom2_label: "H6",
            order: BondOrder::Single,
        },
    ],
};

pub const THYMINE: MoleculeTestCase = MoleculeTestCase {
    name: "Thymine",
    atoms: &[
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N3",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "O4",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H6",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H71",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H72",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H73",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        BondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "N3",
            atom2_label: "C4",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C6",
            atom2_label: "N1",
            order: BondOrder::Aromatic,
        },
        BondBlueprint {
            atom1_label: "C2",
            atom2_label: "O2",
            order: BondOrder::Double,
        },
        BondBlueprint {
            atom1_label: "C4",
            atom2_label: "O4",
            order: BondOrder::Double,
        },
        BondBlueprint {
            atom1_label: "C5",
            atom2_label: "C7",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "N3",
            atom2_label: "H3",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C6",
            atom2_label: "H6",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C7",
            atom2_label: "H71",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C7",
            atom2_label: "H72",
            order: BondOrder::Single,
        },
        BondBlueprint {
            atom1_label: "C7",
            atom2_label: "H73",
            order: BondOrder::Single,
        },
    ],
};
