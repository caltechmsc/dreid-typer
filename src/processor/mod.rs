//! Molecular perception and typing pipeline for the DREIDING force field.
//!
//! This module implements the core processing pipeline that transforms a basic molecular graph
//! into a fully typed molecular topology. It consists of perception algorithms for chemical
//! properties and a rule-based typing engine for atom type assignment.

mod graph;
mod perception;
mod pipeline;
mod templates;
mod typer;

pub(crate) use graph::{AtomView, ProcessingGraph};
pub(crate) use pipeline::perceive;
pub(crate) use typer::assign_types;
