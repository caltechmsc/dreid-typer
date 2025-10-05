use super::graph::{AtomView, ProcessingGraph};
use crate::core::Element;
use crate::core::error::AssignmentError;
use crate::rules::{Conditions, Rule};
use std::collections::HashMap;

pub(crate) struct TyperEngine<'a> {
    graph: &'a ProcessingGraph,
    rules: Vec<&'a Rule>,
    atom_states: Vec<Option<(String, i32)>>,
    rounds_completed: u32,
}
