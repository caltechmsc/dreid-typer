use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TyperError {
    /// An error occurred while parsing the user-provided or default rule set.
    /// This indicates a syntax error in the TOML rule file.
    RuleParse(String),

    /// The input `MolecularGraph` is structurally invalid or inconsistent.
    /// This error occurs before any typing logic is applied.
    InvalidInputGraph(GraphValidationError),

    /// The chemical feature perception pipeline failed for a specific reason.
    /// This indicates a logic error in our perception rules or unhandled chemical environments.
    AnnotationFailed(AnnotationError),

    /// The iterative typing engine failed to assign a DREIDING type to all atoms.
    /// This is the primary failure mode of the core algorithm, indicating that
    /// the provided rules are insufficient or ambiguous for the given molecule.
    AssignmentFailed(AssignmentError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphValidationError {
    /// An atom ID was referenced (e.g., in a bond) but not defined in the atom list.
    MissingAtom { id: usize },
    /// Two or more atoms were defined with the same unique ID.
    DuplicateAtomId { id: usize },
    /// An atom was defined to be bonded to itself.
    SelfBondingAtom { id: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssignmentError {
    /// A list of unique identifiers of all atoms that remained untyped.
    pub untyped_atom_ids: Vec<usize>,
    /// The total number of iterative rounds completed before the engine stalled.
    pub rounds_completed: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnotationError {
    /// The hybridization inference logic could not determine a state for an atom.
    HybridizationInference { atom_id: usize },
    /// A generic error message for other potential failures in the pipeline.
    Other(String),
}

impl fmt::Display for TyperError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RuleParse(msg) => write!(f, "Rule parsing error: {}", msg),
            Self::InvalidInputGraph(err) => write!(f, "Invalid input graph: {}", err),
            Self::AnnotationFailed(err) => write!(f, "Chemical annotation failed: {}", err),
            Self::AssignmentFailed(err) => write!(f, "Atom typing failed: {}", err),
        }
    }
}

impl fmt::Display for GraphValidationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingAtom { id } => {
                write!(f, "a bond references a non-existent atom with ID {}", id)
            }
            Self::DuplicateAtomId { id } => {
                write!(f, "found duplicate definition for atom with ID {}", id)
            }
            Self::SelfBondingAtom { id } => write!(f, "atom with ID {} is bonded to itself", id),
        }
    }
}

impl fmt::Display for AssignmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "engine stalled after {} rounds with {} untyped atoms remaining. IDs: {:?}",
            self.rounds_completed,
            self.untyped_atom_ids.len(),
            self.untyped_atom_ids
        )
    }
}

impl fmt::Display for AnnotationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::HybridizationInference { atom_id } => {
                write!(
                    f,
                    "could not infer hybridization for atom with ID {}",
                    atom_id
                )
            }
            Self::Other(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::error::Error for TyperError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidInputGraph(e) => Some(e),
            Self::AnnotationFailed(e) => Some(e),
            Self::AssignmentFailed(e) => Some(e),
            _ => None,
        }
    }
}
impl std::error::Error for GraphValidationError {}
impl std::error::Error for AssignmentError {}
impl std::error::Error for AnnotationError {}
