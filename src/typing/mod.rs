//! Hosts the DREIDING typing pipeline, including rule parsing and rule application engines.
//!
//! This namespace exposes the rule schema (`rules`) and the iterative assignment engine
//! (`engine`) used by `assign_topology`.

/// Typing engine that evaluates rules over annotated molecules.
pub mod engine;
/// Rule definitions and parsing utilities.
pub mod rules;
