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
