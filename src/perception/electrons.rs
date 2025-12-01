//! Assigns formal charges and lone pairs by recognizing charged motifs before falling back to
//! element-based heuristics.
//!
//! The routines scan for well-known functional groups (nitro, sulfone, ammonium, etc.), mark their
//! atoms as processed to avoid double counting, and finally run a general valence-based pass for the
//! remaining atoms.

use super::model::AnnotatedMolecule;
use crate::core::error::PerceptionError;
use crate::core::properties::{Element, GraphBondOrder};

/// Runs the electron perception pipeline on the annotated molecule.
///
/// Specialized passes label recognizable charged fragments first so the generic valence pass can
/// operate only on leftover atoms without clobbering explicit assignments.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose atoms will receive charges and lone pairs.
///
/// # Returns
///
/// `Ok(())` once all perception passes finish successfully.
///
/// # Errors
///
/// Propagates any [`PerceptionError`] emitted by helper routines (currently only `assign_general`
/// can error if an element lacks valence data).
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let mut processed = vec![false; molecule.atoms.len()];

    assign_nitrone_groups(molecule, &mut processed)?;
    assign_nitro_groups(molecule, &mut processed)?;
    assign_sulfur_oxides(molecule, &mut processed)?;
    assign_halogen_oxyanions(molecule, &mut processed)?;
    assign_phosphorus_oxides(molecule, &mut processed)?;
    assign_carboxylate_anions(molecule, &mut processed)?;
    assign_ammonium_and_iminium(molecule, &mut processed)?;
    assign_onium_ions(molecule, &mut processed)?;
    assign_phosphonium_ions(molecule, &mut processed)?;
    assign_enolate_phenate_anions(molecule, &mut processed)?;

    assign_general(molecule, &processed)?;

    Ok(())
}

/// Detects nitrones and applies the canonical charge distribution.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to inspect and mutate.
/// * `processed` - Scratch mask indicating atoms already assigned by previous passes.
///
/// # Returns
///
/// `Ok(())`; this routine never emits an error.
fn assign_nitrone_groups(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx]
            || molecule.atoms[n_idx].element != Element::N
            || molecule.atoms[n_idx].degree != 3
        {
            continue;
        }

        let mut double_bond_c_idx = None;
        let mut single_bond_o_idx = None;
        let mut single_bond_c_idx = None;

        for &(neighbor_id, order) in &molecule.adjacency[n_idx] {
            match (molecule.atoms[neighbor_id].element, order) {
                (Element::C, GraphBondOrder::Double) => double_bond_c_idx = Some(neighbor_id),
                (Element::O, GraphBondOrder::Single) => single_bond_o_idx = Some(neighbor_id),
                (Element::C, GraphBondOrder::Single) => single_bond_c_idx = Some(neighbor_id),
                _ => {}
            }
        }

        match (double_bond_c_idx, single_bond_o_idx, single_bond_c_idx) {
            (Some(c1), Some(o), Some(c2)) if !processed[c1] && !processed[o] && !processed[c2] => {
                molecule.atoms[n_idx].formal_charge = 1;
                molecule.atoms[n_idx].lone_pairs = 0;
                molecule.atoms[o].formal_charge = -1;
                molecule.atoms[o].lone_pairs = 3;

                processed[n_idx] = true;
                processed[o] = true;
            }
            _ => {}
        }
    }
    Ok(())
}

/// Labels nitro groups with the expected +1/N and -1/O assignment.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to inspect.
/// * `processed` - Mask tracking atoms already handled.
///
/// # Returns
///
/// Always returns `Ok(())`.
fn assign_nitro_groups(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx]
            || molecule.atoms[n_idx].element != Element::N
            || molecule.atoms[n_idx].degree != 3
        {
            continue;
        }

        let mut double_bond_o_idx = None;
        let mut single_bond_o_idx = None;
        let mut other_neighbor_count = 0;

        for &(neighbor_id, order) in &molecule.adjacency[n_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                match order {
                    GraphBondOrder::Double if double_bond_o_idx.is_none() => {
                        double_bond_o_idx = Some(neighbor_id)
                    }
                    GraphBondOrder::Single if single_bond_o_idx.is_none() => {
                        single_bond_o_idx = Some(neighbor_id)
                    }
                    _ => continue,
                }
            } else {
                other_neighbor_count += 1;
            }
        }

        match (double_bond_o_idx, single_bond_o_idx) {
            (Some(o1), Some(o2))
                if other_neighbor_count == 1 && !processed[o1] && !processed[o2] =>
            {
                molecule.atoms[n_idx].formal_charge = 1;
                molecule.atoms[n_idx].lone_pairs = 0;

                molecule.atoms[o1].formal_charge = 0;
                molecule.atoms[o1].lone_pairs = 2;

                molecule.atoms[o2].formal_charge = -1;
                molecule.atoms[o2].lone_pairs = 3;

                processed[n_idx] = true;
                processed[o1] = true;
                processed[o2] = true;
            }
            _ => {}
        }
    }
    Ok(())
}

/// Handles sulfoxide and sulfone motifs, setting sulfur and oxygen charges.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule being mutated.
/// * `processed` - Mask updated for sulfur/oxygen atoms already visited.
///
/// # Returns
///
/// `Ok(())` after all sulfur patterns are applied.
fn assign_sulfur_oxides(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for s_idx in 0..molecule.atoms.len() {
        if processed[s_idx] || molecule.atoms[s_idx].element != Element::S {
            continue;
        }

        let oxygen_neighbors: Vec<(usize, GraphBondOrder)> = molecule.adjacency[s_idx]
            .iter()
            .filter(|&&(id, _)| molecule.atoms[id].element == Element::O)
            .cloned()
            .collect();

        let double_bonded_oxygens: Vec<usize> = oxygen_neighbors
            .iter()
            .filter(|&&(_, order)| order == GraphBondOrder::Double)
            .map(|&(id, _)| id)
            .collect();

        if molecule.atoms[s_idx].degree == 3 && oxygen_neighbors.len() == 1 {
            let (o_idx, _) = oxygen_neighbors[0];
            if !processed[o_idx] {
                molecule.atoms[s_idx].formal_charge = 0;
                molecule.atoms[s_idx].lone_pairs = 1;
                molecule.atoms[o_idx].formal_charge = 0;
                molecule.atoms[o_idx].lone_pairs = 2;
                processed[s_idx] = true;
                processed[o_idx] = true;
            }
        } else if molecule.atoms[s_idx].degree == 4 && double_bonded_oxygens.len() == 2 {
            let o1_idx = double_bonded_oxygens[0];
            let o2_idx = double_bonded_oxygens[1];
            if !processed[o1_idx] && !processed[o2_idx] {
                molecule.atoms[s_idx].formal_charge = 2;
                molecule.atoms[s_idx].lone_pairs = 0;
                molecule.atoms[o1_idx].formal_charge = -1;
                molecule.atoms[o1_idx].lone_pairs = 3;
                molecule.atoms[o2_idx].formal_charge = -1;
                molecule.atoms[o2_idx].lone_pairs = 3;
                processed[s_idx] = true;
                processed[o1_idx] = true;
                processed[o2_idx] = true;
            }
        }
    }
    Ok(())
}

/// Assigns charges to halogen oxyanion oxygens, distinguishing single vs. double bonds.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule under inspection.
/// * `processed` - Mutable mask used to avoid reassigning oxygens.
///
/// # Returns
///
/// `Ok(())` when the pass completes.
fn assign_halogen_oxyanions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for center_idx in 0..molecule.atoms.len() {
        if processed[center_idx] {
            continue;
        }

        if !matches!(
            molecule.atoms[center_idx].element,
            Element::Cl | Element::Br | Element::I
        ) {
            continue;
        }

        let oxygen_neighbors: Vec<(usize, GraphBondOrder)> = molecule.adjacency[center_idx]
            .iter()
            .filter(|&&(neighbor_id, _)| molecule.atoms[neighbor_id].element == Element::O)
            .map(|&(neighbor_id, order)| (neighbor_id, order))
            .collect();

        if oxygen_neighbors.len() < 3 {
            continue;
        }

        for &(oxygen_idx, order) in &oxygen_neighbors {
            if processed[oxygen_idx] {
                continue;
            }

            let oxygen = &mut molecule.atoms[oxygen_idx];
            oxygen.lone_pairs = 2;
            oxygen.formal_charge = if order == GraphBondOrder::Single {
                -1
            } else {
                0
            };

            processed[oxygen_idx] = true;
        }
    }

    Ok(())
}

/// Captures phosphoryl species with a single P=O double bond and applies +1/-1 charges.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to mutate.
/// * `processed` - Mask marking atoms that should not be revisited later.
///
/// # Returns
///
/// Always returns `Ok(())`.
fn assign_phosphorus_oxides(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for p_idx in 0..molecule.atoms.len() {
        if processed[p_idx]
            || molecule.atoms[p_idx].element != Element::P
            || molecule.atoms[p_idx].degree != 4
        {
            continue;
        }

        let double_bonded_oxygens: Vec<usize> = molecule.adjacency[p_idx]
            .iter()
            .filter(|&&(id, order)| {
                molecule.atoms[id].element == Element::O && order == GraphBondOrder::Double
            })
            .map(|&(id, _)| id)
            .collect();

        if double_bonded_oxygens.len() == 1 {
            let o_idx = double_bonded_oxygens[0];
            if !processed[o_idx] {
                molecule.atoms[p_idx].formal_charge = 1;
                molecule.atoms[p_idx].lone_pairs = 0;
                molecule.atoms[o_idx].formal_charge = -1;
                molecule.atoms[o_idx].lone_pairs = 3;
                processed[p_idx] = true;
                processed[o_idx] = true;
            }
        }
    }
    Ok(())
}

/// Detects carboxylate anions and assigns the single-bonded oxygen a -1 charge.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to examine.
/// * `processed` - Mask updated for atoms within the recognized carboxylate.
///
/// # Returns
///
/// `Ok(())` when complete.
fn assign_carboxylate_anions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx]
            || molecule.atoms[c_idx].element != Element::C
            || molecule.atoms[c_idx].degree != 3
        {
            continue;
        }

        let mut double_bond_o_idx = None;
        let mut single_bond_o_idx = None;

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                match order {
                    GraphBondOrder::Double if double_bond_o_idx.is_none() => {
                        double_bond_o_idx = Some(neighbor_id)
                    }
                    GraphBondOrder::Single if single_bond_o_idx.is_none() => {
                        single_bond_o_idx = Some(neighbor_id)
                    }
                    _ => continue,
                }
            }
        }

        match (double_bond_o_idx, single_bond_o_idx) {
            (Some(o1), Some(o2))
                if !processed[o1] && !processed[o2] && molecule.atoms[o2].degree == 1 =>
            {
                molecule.atoms[c_idx].formal_charge = 0;
                molecule.atoms[c_idx].lone_pairs = 0;

                molecule.atoms[o1].formal_charge = 0;
                molecule.atoms[o1].lone_pairs = 2;

                molecule.atoms[o2].formal_charge = -1;
                molecule.atoms[o2].lone_pairs = 3;

                processed[c_idx] = true;
                processed[o1] = true;
                processed[o2] = true;
            }
            _ => {}
        }
    }
    Ok(())
}

/// Marks ammonium/iminium nitrogens as positively charged when geometry criteria match.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule containing candidate nitrogens.
/// * `processed` - Mask tracking which atoms have already been finalized.
///
/// # Returns
///
/// `Ok(())`; the function never errors.
fn assign_ammonium_and_iminium(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (n_idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag || molecule.atoms[n_idx].element != Element::N {
            continue;
        }

        let degree = molecule.atoms[n_idx].degree;
        let has_double_bond = molecule.adjacency[n_idx]
            .iter()
            .any(|&(_, order)| order == GraphBondOrder::Double);

        if degree == 4 || (degree == 3 && has_double_bond) {
            let atom = &mut molecule.atoms[n_idx];
            atom.formal_charge = 1;
            atom.lone_pairs = 0;
            *processed_flag = true;
        }
    }
    Ok(())
}

/// Handles oxonium/sulfonium (onium) ions by enforcing +1 charge and lone pairs.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule.
/// * `processed` - Mask preventing repeated assignment.
///
/// # Returns
///
/// `Ok(())` after the pass finishes.
fn assign_onium_ions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag {
            continue;
        }

        let element = molecule.atoms[idx].element;
        let degree = molecule.atoms[idx].degree;
        let has_pi_bond = molecule.adjacency[idx]
            .iter()
            .any(|&(_, order)| order != GraphBondOrder::Single);

        if (element == Element::O || element == Element::S) && degree == 3 && !has_pi_bond {
            let atom_mut = &mut molecule.atoms[idx];
            atom_mut.formal_charge = 1;
            atom_mut.lone_pairs = 1;
            *processed_flag = true;
        }
    }
    Ok(())
}

/// Detects tetravalent phosphorus centers lacking P=O and labels phosphonium ions.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule being mutated.
/// * `processed` - Mask ensuring atoms are only assigned once.
///
/// # Returns
///
/// Always returns `Ok(())`.
fn assign_phosphonium_ions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (p_idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag
            || molecule.atoms[p_idx].element != Element::P
            || molecule.atoms[p_idx].degree != 4
        {
            continue;
        }

        let has_double_bond_o = molecule.adjacency[p_idx].iter().any(|&(id, order)| {
            molecule.atoms[id].element == Element::O && order == GraphBondOrder::Double
        });

        if !has_double_bond_o {
            let atom = &mut molecule.atoms[p_idx];
            atom.formal_charge = 1;
            atom.lone_pairs = 0;
            *processed_flag = true;
        }
    }
    Ok(())
}

/// Marks enolates and phenates by giving terminal oxygens a -1 charge when attached to sp² carbon.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule being inspected.
/// * `processed` - Mask recording atoms that no longer need processing.
///
/// # Returns
///
/// `Ok(())` when done.
fn assign_enolate_phenate_anions(
    molecule: &mut AnnotatedMolecule,
    processed: &mut [bool],
) -> Result<(), PerceptionError> {
    for (o_idx, processed_flag) in processed.iter_mut().enumerate() {
        if *processed_flag
            || molecule.atoms[o_idx].element != Element::O
            || molecule.atoms[o_idx].degree != 1
        {
            continue;
        }

        let (neighbor_id, order) = molecule.adjacency[o_idx][0];

        if order != GraphBondOrder::Single {
            continue;
        }

        let neighbor = &molecule.atoms[neighbor_id];
        if neighbor.element == Element::C {
            let neighbor_is_sp2 = molecule.adjacency[neighbor_id]
                .iter()
                .any(|&(_, order)| order == GraphBondOrder::Double);

            if neighbor_is_sp2 {
                let atom = &mut molecule.atoms[o_idx];
                atom.formal_charge = -1;
                atom.lone_pairs = 3;
                *processed_flag = true;
            }
        }
    }
    Ok(())
}

/// Fallback valence-based assignment for atoms not matched by specialized rules.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose remaining atoms get charges/lone pairs.
/// * `processed` - Mask denoting which atoms should be skipped (already solved).
///
/// # Returns
///
/// `Ok(())` once complete.
///
/// # Errors
///
/// Returns [`PerceptionError::Other`] if an element lacks a `valence_electrons` definition.
fn assign_general(
    molecule: &mut AnnotatedMolecule,
    processed: &[bool],
) -> Result<(), PerceptionError> {
    for (i, processed_flag) in processed.iter().enumerate() {
        if *processed_flag {
            continue;
        }
        let element = molecule.atoms[i].element;

        let valence = element.valence_electrons().ok_or_else(|| {
            PerceptionError::Other(format!(
                "valence electrons not defined for element {:?}",
                element
            ))
        })?;

        let bonding_electrons: u8 = molecule.adjacency[i]
            .iter()
            .map(|&(_, order)| bond_order_to_valence(order))
            .sum();

        let mut lone_pairs = 0;

        let is_second_period = matches!(
            element,
            Element::B | Element::C | Element::N | Element::O | Element::F
        );

        if element == Element::H {
            let bonded_electrons = bonding_electrons.saturating_mul(2);
            if bonded_electrons <= 2 {
                lone_pairs = (2 - bonded_electrons) / 2;
            }
        } else if is_second_period {
            let bonded_electrons = bonding_electrons.saturating_mul(2);
            if bonded_electrons <= 8 {
                lone_pairs = (8 - bonded_electrons) / 2;
            }
        } else if valence >= bonding_electrons {
            lone_pairs = (valence - bonding_electrons) / 2;
        }

        let formal_charge = valence as i8 - bonding_electrons as i8 - (lone_pairs * 2) as i8;

        let atom_mut = &mut molecule.atoms[i];
        atom_mut.lone_pairs = lone_pairs;
        atom_mut.formal_charge = formal_charge;
    }
    Ok(())
}

/// Converts a bond order into its valence contribution for generic bookkeeping.
///
/// # Arguments
///
/// * `order` - Bond order encountered when summing bonding electrons.
///
/// # Returns
///
/// Valence contribution counted toward an atom's bonding electron total.
fn bond_order_to_valence(order: GraphBondOrder) -> u8 {
    match order {
        GraphBondOrder::Single => 1,
        GraphBondOrder::Double => 2,
        GraphBondOrder::Triple => 3,
        GraphBondOrder::Aromatic => panic!("Aromatic bonds should have been kekulized by now."),
    }
}
