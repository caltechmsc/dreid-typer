use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TyperError {
    /// An error occurred while parsing the user-provided or default rule set.
    /// This indicates a syntax error in the TOML rule file.
    RuleParse(String),

    /// The input `MolecularGraph` is structurally invalid or inconsistent.
    /// This error occurs before any typing logic is applied.
    InvalidInputGraph(GraphValidationError),

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
