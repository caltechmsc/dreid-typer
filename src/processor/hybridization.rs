use super::graph::{AtomView, ProcessingGraph};
use crate::core::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ProvisionalHybridization {
    SP,
    SP2,
    SP3,
}

pub(super) fn calculate_provisional_hybridization(
    atom: &AtomView,
    graph: &ProcessingGraph,
) -> ProvisionalHybridization {
    let pi_electrons_from_bonds = graph.adjacency[atom.id]
        .iter()
        .map(|(_, order)| order.pi_contribution())
        .sum::<u8>();

    let mut prov_hyb = match pi_electrons_from_bonds {
        2.. => ProvisionalHybridization::SP,
        1 => ProvisionalHybridization::SP2,
        0 => ProvisionalHybridization::SP3,
    };

    if prov_hyb == ProvisionalHybridization::SP3 {
        match atom.element {
            Element::N | Element::P | Element::As if atom.degree == 3 => {
                prov_hyb = ProvisionalHybridization::SP2;
            }
            Element::O | Element::S | Element::Se if atom.degree == 2 => {
                prov_hyb = ProvisionalHybridization::SP2;
            }
            Element::B if atom.degree == 3 => {
                prov_hyb = ProvisionalHybridization::SP2;
            }
            _ => {}
        }
    }

    prov_hyb
}

pub(crate) fn infer_hybridization_for_all(graph: &mut ProcessingGraph) {
    for i in 0..graph.atoms.len() {
        let final_hyb = infer_single_hybridization(&graph.atoms[i], graph);
        graph.atoms[i].hybridization = final_hyb;
    }

    debug_assert!(
        !graph
            .atoms
            .iter()
            .any(|a| a.hybridization == Hybridization::Unknown),
        "Hybridization inference failed: one or more atoms remain Unknown."
    );
}

fn infer_single_hybridization(atom: &AtomView, graph: &ProcessingGraph) -> Hybridization {
    if atom.is_aromatic {
        return Hybridization::Resonant;
    }

    if atom.degree == 0 {
        match atom.element {
            Element::Na | Element::Ca | Element::Fe | Element::Zn => return Hybridization::None,
            _ => return Hybridization::None,
        }
    }

    match atom.element {
        Element::H
        | Element::F
        | Element::Cl
        | Element::Br
        | Element::I
        | Element::He
        | Element::Ne
        | Element::Ar
        | Element::Kr
        | Element::Xe
        | Element::Na
        | Element::Ca
        | Element::Fe
        | Element::Zn => {
            return Hybridization::None;
        }
        _ => {}
    }

    let pi_electrons_from_bonds = graph.adjacency[atom.id]
        .iter()
        .map(|(_, order)| order.pi_contribution())
        .sum::<u8>();

    match pi_electrons_from_bonds {
        2.. => Hybridization::SP,
        1 => Hybridization::SP2,
        0 => match atom.element {
            Element::B if atom.degree == 3 => Hybridization::SP2,
            _ => Hybridization::SP3,
        },
    }
}

impl BondOrder {
    pub(super) fn pi_contribution(&self) -> u8 {
        match self {
            BondOrder::Single => 0,
            BondOrder::Double => 1,
            BondOrder::Triple => 2,
            BondOrder::Aromatic => 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;

    fn setup_test_graph(
        elements: Vec<Element>,
        bonds: Vec<(usize, usize, BondOrder)>,
    ) -> ProcessingGraph {
        let mut mg = MolecularGraph::new();
        for element in elements {
            mg.add_atom(element);
        }
        for (u, v, order) in bonds {
            mg.add_bond(u, v, order).unwrap();
        }
        ProcessingGraph::new(&mg).unwrap()
    }

    #[test]
    fn pi_contribution_from_single_bond_is_zero() {
        assert_eq!(BondOrder::Single.pi_contribution(), 0);
    }

    #[test]
    fn pi_contribution_from_double_bond_is_one() {
        assert_eq!(BondOrder::Double.pi_contribution(), 1);
    }

    #[test]
    fn pi_contribution_from_triple_bond_is_two() {
        assert_eq!(BondOrder::Triple.pi_contribution(), 2);
    }

    #[test]
    fn pi_contribution_from_aromatic_bond_is_one() {
        assert_eq!(BondOrder::Aromatic.pi_contribution(), 1);
    }

    #[test]
    fn infer_hybridization_for_methane_carbon_is_sp3() {
        let graph = setup_test_graph(
            vec![Element::C, Element::H, Element::H, Element::H, Element::H],
            vec![
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
            ],
        );
        let carbon = &graph.atoms[0];
        assert_eq!(
            infer_single_hybridization(carbon, &graph),
            Hybridization::SP3
        );
    }

    #[test]
    fn infer_hybridization_for_ethene_carbon_is_sp2() {
        let graph = setup_test_graph(
            vec![
                Element::C,
                Element::C,
                Element::H,
                Element::H,
                Element::H,
                Element::H,
            ],
            vec![
                (0, 1, BondOrder::Double),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
            ],
        );
        let carbon = &graph.atoms[0];
        assert_eq!(
            infer_single_hybridization(carbon, &graph),
            Hybridization::SP2
        );
    }

    #[test]
    fn infer_hybridization_for_ethyne_carbon_is_sp() {
        let graph = setup_test_graph(
            vec![Element::C, Element::C, Element::H, Element::H],
            vec![
                (0, 1, BondOrder::Triple),
                (0, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
            ],
        );
        let carbon = &graph.atoms[0];
        assert_eq!(
            infer_single_hybridization(carbon, &graph),
            Hybridization::SP
        );
    }

    #[test]
    fn infer_hybridization_for_boron_trifluoride_boron_is_sp2() {
        let graph = setup_test_graph(
            vec![Element::B, Element::F, Element::F, Element::F],
            vec![
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let boron = &graph.atoms[0];
        assert_eq!(
            infer_single_hybridization(boron, &graph),
            Hybridization::SP2
        );
    }

    #[test]
    fn infer_hybridization_for_aromatic_atom_is_resonant() {
        let mut graph = setup_test_graph(vec![Element::C], vec![]);
        graph.atoms[0].is_aromatic = true;
        let carbon = &graph.atoms[0];
        assert_eq!(
            infer_single_hybridization(carbon, &graph),
            Hybridization::Resonant
        );
    }

    #[test]
    fn infer_hybridization_for_halogen_is_none() {
        let graph = setup_test_graph(
            vec![Element::C, Element::F],
            vec![(0, 1, BondOrder::Single)],
        );
        let fluorine = &graph.atoms[1];
        assert_eq!(
            infer_single_hybridization(fluorine, &graph),
            Hybridization::None
        );
    }

    #[test]
    fn infer_hybridization_for_isolated_metal_is_none() {
        let graph = setup_test_graph(vec![Element::Na], vec![]);
        let sodium = &graph.atoms[0];
        assert_eq!(
            infer_single_hybridization(sodium, &graph),
            Hybridization::None
        );
    }

    #[test]
    fn calculate_provisional_hybridization_for_sp3_nitrogen_is_sp2() {
        let graph = setup_test_graph(
            vec![Element::N, Element::H, Element::H, Element::H],
            vec![
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let nitrogen = &graph.atoms[0];
        assert_eq!(
            calculate_provisional_hybridization(nitrogen, &graph),
            ProvisionalHybridization::SP2
        );
    }

    #[test]
    fn calculate_provisional_hybridization_for_sp3_oxygen_is_sp2() {
        let graph = setup_test_graph(
            vec![Element::O, Element::H, Element::H],
            vec![(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        let oxygen = &graph.atoms[0];
        assert_eq!(
            calculate_provisional_hybridization(oxygen, &graph),
            ProvisionalHybridization::SP2
        );
    }

    #[test]
    fn calculate_provisional_hybridization_for_sp3_boron_is_sp2() {
        let graph = setup_test_graph(
            vec![Element::B, Element::H, Element::H, Element::H],
            vec![
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let boron = &graph.atoms[0];
        assert_eq!(
            calculate_provisional_hybridization(boron, &graph),
            ProvisionalHybridization::SP2
        );
    }

    #[test]
    fn calculate_provisional_hybridization_for_double_bond_is_sp2() {
        let graph = setup_test_graph(
            vec![Element::C, Element::C],
            vec![(0, 1, BondOrder::Double)],
        );
        let carbon = &graph.atoms[0];
        assert_eq!(
            calculate_provisional_hybridization(carbon, &graph),
            ProvisionalHybridization::SP2
        );
    }

    #[test]
    fn calculate_provisional_hybridization_for_triple_bond_is_sp() {
        let graph = setup_test_graph(
            vec![Element::C, Element::C],
            vec![(0, 1, BondOrder::Triple)],
        );
        let carbon = &graph.atoms[0];
        assert_eq!(
            calculate_provisional_hybridization(carbon, &graph),
            ProvisionalHybridization::SP
        );
    }

    #[test]
    fn infer_hybridization_for_all_in_water() {
        let mut graph = setup_test_graph(
            vec![Element::O, Element::H, Element::H],
            vec![(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        infer_hybridization_for_all(&mut graph);
        assert_eq!(graph.atoms[0].hybridization, Hybridization::SP3);
        assert_eq!(graph.atoms[1].hybridization, Hybridization::None);
        assert_eq!(graph.atoms[2].hybridization, Hybridization::None);
    }
}
