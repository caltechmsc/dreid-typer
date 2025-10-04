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
