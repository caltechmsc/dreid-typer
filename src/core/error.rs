use thiserror::Error;

/// The primary error type for all operations within the `dreid-typer` library.
#[derive(Debug, Error)]
pub enum TyperError {
    /// An error occurred during the validation of the input `MolecularGraph`.
    #[error("invalid input graph")]
    InvalidInput(#[from] GraphValidationError),

    /// An error occurred while parsing the DREIDING typing rules from a TOML string.
    #[error("failed to parse typing rules")]
    RuleParse(#[from] toml::de::Error),

    /// A failure occurred during a specific step of the chemical perception pipeline.
    #[error("chemical perception failed during '{step}' step")]
    PerceptionFailed {
        step: String,
        #[source]
        source: PerceptionError,
    },

    /// The typing engine failed to assign a type to one or more atoms.
    #[error("atom typing failed")]
    AssignmentFailed(#[from] AssignmentError),
}

/// Errors related to the structural or logical integrity of the input `MolecularGraph`.
#[derive(Debug, Error)]
pub enum GraphValidationError {
    /// A bond references an atom ID that does not exist.
    #[error("bond references a non-existent atom with ID {atom_id}")]
    MissingAtom { atom_id: usize },

    /// An atom is defined to be bonded to itself.
    #[error("atom with ID {atom_id} is bonded to itself")]
    SelfBondingAtom { atom_id: usize },
}

/// Errors that can occur during the multi-step chemical perception phase.
#[derive(Debug, Error)]
pub enum PerceptionError {
    /// Failed to find a valid Kekulé structure for the aromatic systems in the molecule.
    #[error("failed to assign a valid Kekulé structure to the aromatic systems: {message}")]
    KekulizationFailed { message: String },

    /// Failed to determine the hybridization state for an atom.
    #[error(
        "could not infer hybridization for atom ID {atom_id}: unhandled steric number or chemical environment"
    )]
    HybridizationInference { atom_id: usize },

    /// An error originated from the external `pauling` library during conjugation perception.
    #[error("conjugation perception failed via the 'pauling' library")]
    PaulingError(#[from] pauling::PerceptionError),

    /// A generic error for other, less common perception failures.
    #[error("an unexpected perception error occurred: {0}")]
    Other(String),
}

/// An error indicating that the typing engine could not assign types to all atoms.
///
/// This typically occurs when the rule set is incomplete for the given molecule.
#[derive(Debug, Error)]
#[error("engine stalled after {rounds_completed} rounds with {untyped_atom_ids:?} still untyped")]
pub struct AssignmentError {
    /// A list of unique identifiers of all atoms that remained untyped.
    pub untyped_atom_ids: Vec<usize>,
    /// The total number of iterative rounds completed before the engine stalled.
    pub rounds_completed: u32,
}
