mod graph;
mod perception;
mod pipeline;
mod templates;
mod typer;

pub(crate) use graph::{AtomView, ProcessingGraph};
pub(crate) use pipeline::perceive;
pub(crate) use typer::assign_types;
