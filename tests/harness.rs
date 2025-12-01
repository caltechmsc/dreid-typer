#[path = "cases/mod.rs"]
pub mod cases;

use dreid_typer::{
    Element, GraphBondOrder, MolecularGraph, MolecularTopology, TopologyBondOrder, assign_topology,
};
use std::collections::{HashMap, HashSet};

#[derive(Debug)]
pub struct AtomBlueprint {
    pub label: &'static str,
    pub element: Element,
    pub expected_type: &'static str,
}

#[derive(Debug)]
pub struct InputBondBlueprint {
    pub atom1_label: &'static str,
    pub atom2_label: &'static str,
    pub order: GraphBondOrder,
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct OutputBondBlueprint {
    pub atom1_label: &'static str,
    pub atom2_label: &'static str,
    pub order: TopologyBondOrder,
}

#[derive(Debug)]
pub struct MoleculeTestCase {
    pub name: &'static str,
    pub atoms: &'static [AtomBlueprint],
    pub bonds: &'static [InputBondBlueprint],
    pub expected_bonds: &'static [OutputBondBlueprint],
}

pub struct LabeledMolecule {
    graph: MolecularGraph,
    labels: HashMap<&'static str, usize>,
}

impl LabeledMolecule {
    pub fn graph(&self) -> &MolecularGraph {
        &self.graph
    }

    pub fn id(&self, label: &'static str) -> usize {
        *self
            .labels
            .get(label)
            .unwrap_or_else(|| panic!("Unknown atom label: {}", label))
    }
}

pub fn run_molecule_test_case(case: &MoleculeTestCase) {
    let molecule = build_from_blueprint(case);

    let topology = assign_topology(molecule.graph())
        .unwrap_or_else(|err| panic!("Topology assignment failed for '{}': {:?}", case.name, err));

    verify_atom_types(&topology, &molecule, case);
    verify_bond_orders(&topology, &molecule, case);
}

fn build_from_blueprint(case: &MoleculeTestCase) -> LabeledMolecule {
    let mut graph = MolecularGraph::new();
    let mut labels = HashMap::new();

    for atom_bp in case.atoms {
        let id = graph.add_atom(atom_bp.element);
        if labels.insert(atom_bp.label, id).is_some() {
            panic!(
                "Molecule '{}': Duplicate atom label '{}'",
                case.name, atom_bp.label
            );
        }
    }

    for bond_bp in case.bonds {
        let id1 = *labels
            .get(bond_bp.atom1_label)
            .unwrap_or_else(|| panic!("Label '{}' not found", bond_bp.atom1_label));
        let id2 = *labels
            .get(bond_bp.atom2_label)
            .unwrap_or_else(|| panic!("Label '{}' not found", bond_bp.atom2_label));
        graph.add_bond(id1, id2, bond_bp.order).unwrap();
    }

    LabeledMolecule { graph, labels }
}

fn verify_atom_types(
    topology: &MolecularTopology,
    molecule: &LabeledMolecule,
    case: &MoleculeTestCase,
) {
    let actual_types: HashMap<usize, &str> = topology
        .atoms
        .iter()
        .map(|a| (a.id, a.atom_type.as_str()))
        .collect();

    let mut all_heavy_atoms_tested = true;

    for atom_bp in case.atoms {
        let atom_id = molecule.id(atom_bp.label);
        let actual_type = actual_types.get(&atom_id).unwrap_or(&"UNTYPED");

        assert_eq!(
            *actual_type, atom_bp.expected_type,
            "\n --- Test Failure ---\nMolecule: '{}'\nAtom Label: '{}'\n  Expected Type: '{}'\n  Actual Type:   '{}'\n -------------------- \n",
            case.name, atom_bp.label, atom_bp.expected_type, actual_type
        );
    }

    for atom in &topology.atoms {
        if atom.element != Element::H {
            if !case.atoms.iter().any(|ab| molecule.id(ab.label) == atom.id) {
                all_heavy_atoms_tested = false;
                eprintln!(
                    "Warning: Heavy atom with ID {} was not checked in test case '{}'",
                    atom.id, case.name
                );
            }
        }
    }
    assert!(
        all_heavy_atoms_tested,
        "One or more heavy atoms were not defined in the test case blueprint for '{}'",
        case.name
    );
}

fn verify_bond_orders(
    topology: &MolecularTopology,
    molecule: &LabeledMolecule,
    case: &MoleculeTestCase,
) {
    let expected_bonds: HashSet<_> = case
        .expected_bonds
        .iter()
        .map(|bp| {
            let mut ids = [molecule.id(bp.atom1_label), molecule.id(bp.atom2_label)];
            ids.sort_unstable();
            ((ids[0], ids[1]), bp.order)
        })
        .collect();

    let actual_bonds: HashSet<_> = topology
        .bonds
        .iter()
        .map(|bond| (bond.atom_ids, bond.order))
        .collect();

    assert_eq!(
        actual_bonds, expected_bonds,
        "\n --- Test Failure ---\nMolecule: '{}'\nBond order mismatch.\n  Expected: {:?}\n  Actual:   {:?}\n -------------------- \n",
        case.name, expected_bonds, actual_bonds
    );
}
