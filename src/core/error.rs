//! Error types describing the failure modes of graph validation, chemical perception, and typing.
//!
//! These enums aggregate lower-level issues so that library consumers can bubble up a single
//! `TyperError` while still inspecting fine-grained context when needed.

use thiserror::Error;

/// Root error emitted by every fallible operation in the typing pipeline.
///
/// Each variant wraps a more specific error that pinpoints the subsystem that failed, allowing
/// callers to recover or log richer diagnostics without losing ergonomic `Result` signatures.
#[derive(Debug, Error)]
pub enum TyperError {
    /// Input validation of the `MolecularGraph` failed before perception could start.
    #[error("invalid input graph")]
    InvalidInput(#[from] GraphValidationError),

    /// Parsing of the DREIDING typing rules TOML payload did not succeed.
    #[error("failed to parse typing rules")]
    RuleParse(#[from] toml::de::Error),

    /// A specific chemical perception stage reported a failure.
    #[error("chemical perception failed during '{step}' step")]
    PerceptionFailed {
        /// Name of the perception step (e.g., "aromaticity" or "hybridization").
        step: String,
        /// Root perception error that triggered the failure.
        #[source]
        source: PerceptionError,
    },

    /// The typing engine exhausted its rounds before assigning all atom types.
    #[error("atom typing failed")]
    AssignmentFailed(#[from] AssignmentError),
}

/// Errors that describe structural or logical issues with the input `MolecularGraph`.
///
/// These failures are detected before any chemical reasoning is attempted so that malformed inputs
/// can be rejected early with precise diagnostics.
#[derive(Debug, Error)]
pub enum GraphValidationError {
    /// A bond references an atom identifier that is missing from the graph.
    #[error("bond references a non-existent atom with ID {atom_id}")]
    MissingAtom {
        /// Identifier of the atom that could not be found.
        atom_id: usize,
    },

    /// An atom record lists itself as one of its bonded neighbors.
    #[error("atom with ID {atom_id} is bonded to itself")]
    SelfBondingAtom {
        /// Identifier of the atom that incorrectly lists a self-bond.
        atom_id: usize,
    },
}

/// Errors raised while running the staged chemical perception pipeline.
///
/// Each variant corresponds to a logical section of perception so that downstream callers can
/// attribute failures to the relevant chemical heuristic.
#[derive(Debug, Error)]
pub enum PerceptionError {
    /// No valid Kekulé structure satisfied the aromatic subgraph constraints.
    #[error("failed to assign a valid Kekulé structure to the aromatic systems: {message}")]
    KekulizationFailed {
        /// Human-readable reason supplied by the Kekulé resolver.
        message: String,
    },

    /// Hybridization inference could not determine an sp/sp2/sp3 class for an atom.
    #[error(
        "could not infer hybridization for atom ID {atom_id}: unhandled steric number or chemical environment"
    )]
    HybridizationInference {
        /// Identifier of the atom whose steric number could not be mapped.
        atom_id: usize,
    },

    /// Conjugation perception failed within the upstream `pauling` library.
    #[error("conjugation perception failed via the 'pauling' library")]
    PaulingError(#[from] pauling::PerceptionError),

    /// Catch-all variant for perception failures that do not fit the other buckets.
    #[error("an unexpected perception error occurred: {0}")]
    Other(String),
}

/// Error reported when the typing engine stalls before all atoms receive types.
///
/// This typically indicates that the ruleset lacks coverage for the perceived environments or that
/// earlier perception output was incomplete.
#[derive(Debug, Error)]
#[error("engine stalled after {rounds_completed} rounds with {untyped_atom_ids:?} still untyped")]
pub struct AssignmentError {
    /// Unique identifiers of atoms that never converged to a final type.
    pub untyped_atom_ids: Vec<usize>,
    /// Total number of engine rounds completed before stalling.
    pub rounds_completed: u32,
}
