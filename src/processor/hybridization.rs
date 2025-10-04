use super::graph::{AtomView, ProcessingGraph};
use crate::core::{BondOrder, Element, Hybridization};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ProvisionalHybridization {
    SP,
    SP2,
    SP3,
}
