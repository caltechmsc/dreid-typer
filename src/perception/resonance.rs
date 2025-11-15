use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let conjugated_systems =
        pauling::find_resonance_systems(molecule).map_err(PerceptionError::PaulingError)?;

    for system in conjugated_systems {
        for atom_id in system.atoms {
            if let Some(atom) = molecule.atoms.get_mut(atom_id) {
                atom.is_in_conjugated_system = true;
            } else {
                return Err(PerceptionError::Other(format!(
                    "pauling library returned an invalid atom ID ({}) that is out of bounds",
                    atom_id
                )));
            }
        }
    }

    Ok(())
}
