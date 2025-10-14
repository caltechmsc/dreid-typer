//! Error types for the dreid-typer library.
//!
//! This module defines all error types that can be returned by the dreid-typer
//! library functions. These errors cover failures in rule parsing, input validation,
//! chemical perception, and atom typing across the three-phase pipeline.

use std::fmt;

/// The primary error type returned by dreid-typer library functions.
///
/// This enum encompasses all possible failure modes of the library, from input
/// validation to the core typing algorithm. It provides detailed error information
/// to help users diagnose issues with their molecular graphs or rule configurations.
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

/// Errors related to the structural validation of input molecular graphs.
///
/// These errors are detected during the initial validation phase before any
/// chemical processing begins. They indicate fundamental issues with the
/// graph's connectivity or atom definitions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphValidationError {
    /// An atom ID was referenced (e.g., in a bond) but not defined in the atom list.
    MissingAtom {
        /// The ID of the missing atom.
        id: usize,
    },
    /// Two or more atoms were defined with the same unique ID.
    DuplicateAtomId {
        /// The duplicated atom ID.
        id: usize,
    },
    /// An atom was defined to be bonded to itself.
    SelfBondingAtom {
        /// The ID of the self-bonding atom.
        id: usize,
    },
}

/// Error indicating failure of the atom typing engine.
///
/// This error occurs when the iterative rule-based typing algorithm cannot
/// assign DREIDING types to all atoms in the molecule. This typically happens
/// when the rule set is incomplete or when the molecule contains chemical
/// environments not covered by the available rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssignmentError {
    /// A list of unique identifiers of all atoms that remained untyped.
    pub untyped_atom_ids: Vec<usize>,
    /// The total number of iterative rounds completed before the engine stalled.
    pub rounds_completed: u32,
}

/// Errors occurring during the chemical feature perception phase.
///
/// These errors indicate failures in the perception pipeline that annotates
/// the molecular graph with chemical properties like hybridization, aromaticity,
/// and electron counts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnotationError {
    /// The hybridization inference logic could not determine a state for an atom.
    HybridizationInference {
        /// The ID of the atom for which hybridization could not be inferred.
        atom_id: usize,
    },
    /// A generic error message for other potential failures in the pipeline.
    Other(String),
}

impl fmt::Display for TyperError {
    /// Formats the error for display to users.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::TyperError;
    ///
    /// let error = TyperError::RuleParse("invalid TOML".to_string());
    /// assert!(format!("{}", error).contains("Rule parsing error"));
    /// ```
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
    /// Formats the graph validation error for display.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::GraphValidationError;
    ///
    /// let error = GraphValidationError::MissingAtom { id: 5 };
    /// assert!(format!("{}", error).contains("non-existent atom with ID 5"));
    /// ```
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
    /// Formats the assignment error with details about untyped atoms.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::AssignmentError;
    ///
    /// let error = AssignmentError {
    ///     untyped_atom_ids: vec![1, 3],
    ///     rounds_completed: 5,
    /// };
    /// let msg = format!("{}", error);
    /// assert!(msg.contains("5 rounds"));
    /// assert!(msg.contains("2 untyped atoms"));
    /// ```
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
    /// Formats the annotation error for display.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::AnnotationError;
    ///
    /// let error = AnnotationError::HybridizationInference { atom_id: 2 };
    /// assert!(format!("{}", error).contains("atom with ID 2"));
    /// ```
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
