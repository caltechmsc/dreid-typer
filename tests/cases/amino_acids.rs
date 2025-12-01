use super::super::{AtomBlueprint, InputBondBlueprint, MoleculeTestCase, OutputBondBlueprint};
use dreid_typer::{Element, GraphBondOrder, TopologyBondOrder};

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
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
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
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA2",
            order: TopologyBondOrder::Single,
        },
    ],
};
