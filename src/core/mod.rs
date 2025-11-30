//! Foundational primitives shared across graph validation, perception, and typing.
//!
//! The `core` module houses the basic data types—errors, graph containers, and
//! chemical properties—that higher layers build upon when inferring topology.

/// Error types describing validation, perception, and typing failure modes.
pub mod error;
