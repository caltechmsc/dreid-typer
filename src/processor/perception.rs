//! Molecular perception algorithms for the DREIDING force field.
//!
//! This module implements the perception phase of the dreid-typer pipeline, converting a basic
//! `MolecularGraph` into a richly annotated `ProcessingGraph` with electron counts, hybridization,
//! ring membership, and aromaticity information. It combines general chemical algorithms with
//! functional group template matching to prepare molecules for atom typing.

use super::graph::{PerceptionSource, ProcessingGraph, RingInfo};
use crate::core::error::{AnnotationError, TyperError};
use crate::core::graph::MolecularGraph;
use crate::core::{BondOrder, Element, Hybridization};
use std::collections::{HashSet, VecDeque};

/// Calculates valence electrons, bonding electrons, and lone pairs for each atom.
///
/// This function initializes the `ProcessingGraph` with basic electron distribution information
/// required for subsequent perception steps. It computes valence electrons based on element type,
/// bonding electrons from bond orders, and derives lone pairs from the difference.
///
/// # Arguments
///
/// * `molecular_graph` - The input molecular graph containing atoms and bonds.
///
/// # Returns
///
/// A `ProcessingGraph` with electron counts annotated on each atom.
///
/// # Errors
///
/// Returns `TyperError::InvalidInputGraph` if the input graph cannot be converted to a processing graph.
pub(crate) fn perceive_electron_counts(
    molecular_graph: &MolecularGraph,
) -> Result<ProcessingGraph, TyperError> {
    let mut graph = ProcessingGraph::new(molecular_graph).map_err(TyperError::InvalidInputGraph)?;

    for atom in &mut graph.atoms {
        let valence = get_valence_electrons(atom.element).unwrap_or(0);

        atom.valence_electrons = valence;

        let bonding = graph.adjacency[atom.id]
            .iter()
            .map(|(_, order)| bond_order_contribution(*order))
            .sum::<u8>();
        atom.bonding_electrons = bonding;

        // Calculate available electrons after accounting for bonding and formal charge,
        // then assign lone pairs as half of the remaining electrons (octet rule approximation).
        let available = valence as i16 - bonding as i16 - atom.formal_charge as i16;
        let adjusted = available.max(0);
        let lone_pairs = (adjusted / 2) as u8;
        atom.lone_pairs = lone_pairs;
        atom.steric_number = 0;
        atom.hybridization = Hybridization::Unknown;
        atom.is_aromatic = false;
        atom.is_in_ring = false;
        atom.smallest_ring_size = None;
        atom.perception_source = None;
    }

    Ok(graph)
}

/// Identifies all rings in the molecular graph using cycle detection.
///
/// This function employs the Johnson cycle finding algorithm to detect all unique cycles
/// (rings) in the graph, which is essential for aromaticity perception and ring-based properties.
///
/// # Arguments
///
/// * `graph` - The processing graph to analyze for rings.
///
/// # Returns
///
/// A `RingInfo` structure containing all detected rings as sets of atom indices.
pub(crate) fn perceive_rings(graph: &ProcessingGraph) -> RingInfo {
    if graph.atoms.is_empty() {
        return RingInfo::default();
    }

    let mut finder = JohnsonCycleFinder::new(graph);
    let sorted_vec_cycles = finder.find_cycles_internal();
    RingInfo(sorted_vec_cycles)
}

/// Annotates atoms with ring membership and smallest ring size information.
///
/// This function processes the detected rings to mark atoms as being in rings and determines
/// the smallest ring each atom participates in, which is crucial for hybridization and typing rules.
///
/// # Arguments
///
/// * `graph` - The processing graph to annotate with ring information.
/// * `ring_info` - The detected rings from `perceive_rings`.
pub(crate) fn apply_ring_annotations(graph: &mut ProcessingGraph, ring_info: &RingInfo) {
    let mut atom_ring_sizes: Vec<Vec<u8>> = vec![vec![]; graph.atoms.len()];
    for ring in &ring_info.0 {
        let ring_len = ring.len() as u8;
        for &atom_id in ring {
            atom_ring_sizes[atom_id].push(ring_len);
        }
    }

    for (atom_id, atom) in graph.atoms.iter_mut().enumerate() {
        if !atom_ring_sizes[atom_id].is_empty() {
            atom.is_in_ring = true;
            atom.smallest_ring_size = atom_ring_sizes[atom_id].iter().min().copied();
        }
    }
}

/// Applies generic perception algorithms for aromaticity and hybridization.
///
/// This is a convenience function that orchestrates the application of general chemical rules
/// for aromaticity detection and hybridization assignment, skipping atoms already handled by templates.
///
/// # Arguments
///
/// * `graph` - The processing graph to annotate.
/// * `ring_info` - The detected rings for aromaticity analysis.
///
/// # Returns
///
/// An empty result on success.
///
/// # Errors
///
/// Returns `AnnotationError` if hybridization inference fails for any atom.
pub(crate) fn perceive_generic_properties(
    graph: &mut ProcessingGraph,
    ring_info: &RingInfo,
) -> Result<(), AnnotationError> {
    perceive_generic_aromaticity(graph, ring_info)?;
    perceive_generic_hybridization(graph)?;
    Ok(())
}

/// Detects aromatic rings using Hückel's rule and marks aromatic atoms.
///
/// This function analyzes each detected ring for aromaticity based on electron count and bond types,
/// respecting atoms already marked by functional group templates.
///
/// # Arguments
///
/// * `graph` - The processing graph to annotate with aromaticity.
/// * `ring_info` - The detected rings to analyze.
///
/// # Returns
///
/// An empty result on success.
///
/// # Errors
///
/// Returns `AnnotationError` if aromaticity detection encounters invalid electron contributions.
pub(crate) fn perceive_generic_aromaticity(
    graph: &mut ProcessingGraph,
    ring_info: &RingInfo,
) -> Result<(), AnnotationError> {
    let mut aromatic_atoms = HashSet::new();

    for ring_atom_ids in &ring_info.0 {
        if is_ring_aromatic(ring_atom_ids, graph) {
            aromatic_atoms.extend(ring_atom_ids.iter().copied());
        }
    }

    for atom_id in aromatic_atoms {
        if graph.atoms[atom_id].perception_source == Some(PerceptionSource::Template) {
            continue;
        }
        graph.atoms[atom_id].is_aromatic = true;
    }

    Ok(())
}

/// Assigns hybridization states to atoms based on steric number and special cases.
///
/// This function determines hybridization for atoms not already handled by templates, using
/// steric number (degree + lone pairs) and considering aromaticity and special element rules.
///
/// # Arguments
///
/// * `graph` - The processing graph to annotate with hybridization.
///
/// # Returns
///
/// An empty result on success.
///
/// # Errors
///
/// Returns `AnnotationError::HybridizationInference` if an atom's hybridization cannot be determined.
pub(crate) fn perceive_generic_hybridization(
    graph: &mut ProcessingGraph,
) -> Result<(), AnnotationError> {
    for atom in &mut graph.atoms {
        if atom.perception_source.is_some() {
            continue;
        }

        atom.steric_number = atom.degree + atom.lone_pairs;

        let hyb = if atom.is_aromatic {
            Hybridization::Resonant
        } else if is_special_non_hybridized(atom.element) || atom.degree == 0 {
            Hybridization::None
        } else {
            match atom.steric_number {
                4 => Hybridization::SP3,
                3 => Hybridization::SP2,
                2 => Hybridization::SP,
                _ => Hybridization::Unknown,
            }
        };

        if hyb == Hybridization::Unknown {
            return Err(AnnotationError::HybridizationInference { atom_id: atom.id });
        }

        atom.hybridization = hyb;
        atom.perception_source = Some(PerceptionSource::Generic);

        // Adjust steric number based on hybridization for consistency with DREIDING conventions.
        atom.steric_number = match hyb {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            Hybridization::None => atom.degree + atom.lone_pairs,
            Hybridization::Unknown => atom.steric_number,
        };
    }

    Ok(())
}

/// Determines if a ring is aromatic based on Hückel's rule.
///
/// Checks for explicit aromatic bonds first, then falls back to pi electron counting
/// following the 4n+2 rule for cyclic conjugated systems.
///
/// # Arguments
///
/// * `ring_atom_ids` - The atom indices forming the ring.
/// * `graph` - The processing graph containing the atoms.
///
/// # Returns
///
/// `true` if the ring is aromatic, `false` otherwise.
fn is_ring_aromatic(ring_atom_ids: &[usize], graph: &ProcessingGraph) -> bool {
    if ring_atom_ids.len() < 3 {
        return false;
    }

    let ring_atoms: HashSet<usize> = ring_atom_ids.iter().copied().collect();

    let all_bonds_aromatic = ring_atom_ids.iter().all(|&atom_id| {
        graph.adjacency[atom_id].iter().all(|(neighbor, order)| {
            if ring_atoms.contains(neighbor) {
                matches!(order, BondOrder::Aromatic)
            } else {
                true
            }
        })
    });

    if all_bonds_aromatic {
        return true;
    }

    // Count pi electrons contributed by each atom in the ring.
    let mut total_pi_electrons = 0u8;

    for &atom_id in ring_atom_ids {
        match count_atom_pi_contribution(atom_id, &ring_atoms, graph) {
            Some(contribution) => total_pi_electrons += contribution,
            None => return false,
        }
    }

    if total_pi_electrons < 2 {
        return false;
    }

    // Hückel's rule: aromatic if 4n+2 pi electrons.
    (total_pi_electrons - 2).is_multiple_of(4)
}

/// Calculates the pi electron contribution of an atom to its ring's aromaticity.
///
/// Considers the atom's steric number, lone pairs, and exocyclic pi bonds to determine
/// how many pi electrons it contributes to the ring's conjugated system.
///
/// # Arguments
///
/// * `atom_id` - The index of the atom to analyze.
/// * `ring_atoms` - The set of atoms in the ring.
/// * `graph` - The processing graph containing the atom.
///
/// # Returns
///
/// The number of pi electrons contributed, or `None` if the atom cannot participate in aromaticity.
fn count_atom_pi_contribution(
    atom_id: usize,
    ring_atoms: &HashSet<usize>,
    graph: &ProcessingGraph,
) -> Option<u8> {
    let atom = &graph.atoms[atom_id];
    let steric = match atom.perception_source {
        Some(PerceptionSource::Template) => atom.steric_number,
        _ => atom.degree + atom.lone_pairs,
    };

    if steric >= 4 && atom.lone_pairs == 0 {
        return None;
    }

    let has_exocyclic_pi_bond = graph.adjacency[atom_id].iter().any(|(neighbor, order)| {
        !ring_atoms.contains(neighbor) && matches!(order, BondOrder::Double | BondOrder::Triple)
    });

    if has_exocyclic_pi_bond {
        return Some(1);
    }

    match steric {
        2 | 3 => Some(1),
        4 if atom.lone_pairs > 0 => Some(2),
        _ => None,
    }
}

/// Converts a bond order to its electron contribution for bonding calculations.
///
/// Maps bond orders to the number of electrons involved in the bond.
///
/// # Arguments
///
/// * `order` - The bond order to convert.
///
/// # Returns
///
/// The number of bonding electrons contributed by this bond order.
fn bond_order_contribution(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 1,
    }
}

/// Checks if an element typically does not participate in hybridization.
///
/// Certain elements (halogens, noble gases, alkali/alkaline earth metals, transition metals)
/// are treated as non-hybridized in the DREIDING force field context.
///
/// # Arguments
///
/// * `element` - The element to check.
///
/// # Returns
///
/// `true` if the element is considered non-hybridized, `false` otherwise.
fn is_special_non_hybridized(element: Element) -> bool {
    matches!(
        element,
        Element::H
            | Element::F
            | Element::Cl
            | Element::Br
            | Element::I
            | Element::He
            | Element::Ne
            | Element::Ar
            | Element::Kr
            | Element::Xe
            | Element::Rn
            | Element::Li
            | Element::Na
            | Element::K
            | Element::Rb
            | Element::Cs
            | Element::Fr
            | Element::Be
            | Element::Mg
            | Element::Ca
            | Element::Sr
            | Element::Ba
            | Element::Ra
            | Element::Fe
            | Element::Zn
    )
}

/// Returns the number of valence electrons for a given element.
///
/// Uses standard periodic table valence electron counts for common elements.
///
/// # Arguments
///
/// * `element` - The element to query.
///
/// # Returns
///
/// The number of valence electrons, or `None` for unknown elements.
fn get_valence_electrons(element: Element) -> Option<u8> {
    use Element::*;
    match element {
        H => Some(1),
        He => Some(2),
        Li | Na | K | Rb | Cs | Fr => Some(1),
        Be | Mg | Ca | Sr | Ba | Ra => Some(2),
        B | Al | Ga | In | Tl => Some(3),
        C | Si | Ge | Sn | Pb => Some(4),
        N | P | As | Sb | Bi => Some(5),
        O | S | Se | Te | Po => Some(6),
        F | Cl | Br | I | At => Some(7),
        Ne | Ar | Kr | Xe | Rn => Some(8),
        _ => None,
    }
}

/// Cycle finder implementation using Johnson's algorithm.
///
/// This struct encapsulates the state for finding all unique cycles in a molecular graph,
/// ensuring each cycle is reported only once regardless of starting point.
struct JohnsonCycleFinder<'a> {
    graph: &'a ProcessingGraph,
    all_cycles: HashSet<Vec<usize>>,
}

impl<'a> JohnsonCycleFinder<'a> {
    /// Creates a new cycle finder for the given graph.
    ///
    /// # Arguments
    ///
    /// * `graph` - The processing graph to analyze for cycles.
    fn new(graph: &'a ProcessingGraph) -> Self {
        Self {
            graph,
            all_cycles: HashSet::new(),
        }
    }

    /// Finds all unique cycles in the graph.
    ///
    /// Iterates through all possible starting nodes to ensure complete cycle detection.
    ///
    /// # Returns
    ///
    /// A set of all detected cycles, each represented as a sorted vector of atom indices.
    fn find_cycles_internal(&mut self) -> HashSet<Vec<usize>> {
        let num_atoms = self.graph.atoms.len();
        for i in 0..num_atoms {
            self.find_cycles_from_node(i);
        }
        self.all_cycles.clone()
    }

    /// Finds cycles starting from a specific node using BFS.
    ///
    /// Uses a queue-based approach to explore paths from the start node, detecting
    /// cycles when returning to the start with sufficient length.
    ///
    /// # Arguments
    ///
    /// * `start_node` - The atom index to start cycle detection from.
    fn find_cycles_from_node(&mut self, start_node: usize) {
        let mut queue = VecDeque::new();
        queue.push_back(vec![start_node]);

        while let Some(path) = queue.pop_front() {
            let last_node = *path.last().unwrap();

            for (neighbor, _) in &self.graph.adjacency[last_node] {
                if *neighbor < start_node {
                    continue;
                }
                if path.contains(neighbor) {
                    if *neighbor == start_node && path.len() > 2 {
                        let mut cycle = path.clone();
                        cycle.sort_unstable();
                        self.all_cycles.insert(cycle);
                    }
                    continue;
                }

                let mut new_path = path.clone();
                new_path.push(*neighbor);
                queue.push_back(new_path);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::MolecularGraph;
    use crate::core::{BondOrder, Element};
    use crate::processor::templates;
    use std::collections::HashSet;

    #[test]
    fn electron_counts_assign_lone_pairs() {
        let mut mg = MolecularGraph::new();
        let o = mg.add_atom(Element::O, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);

        mg.add_bond(o, h1, BondOrder::Single).unwrap();
        mg.add_bond(o, h2, BondOrder::Single).unwrap();

        let graph = perceive_electron_counts(&mg).unwrap();
        let oxygen = &graph.atoms[o];
        assert_eq!(oxygen.valence_electrons, 6);
        assert_eq!(oxygen.bonding_electrons, 2);
        assert_eq!(oxygen.lone_pairs, 2);
        assert_eq!(oxygen.steric_number, 0);
    }

    #[test]
    fn ring_detection_identifies_simple_cycle() {
        let mut mg = MolecularGraph::new();
        let a = mg.add_atom(Element::C, 0);
        let b = mg.add_atom(Element::C, 0);
        let c = mg.add_atom(Element::C, 0);

        mg.add_bond(a, b, BondOrder::Single).unwrap();
        mg.add_bond(b, c, BondOrder::Single).unwrap();
        mg.add_bond(c, a, BondOrder::Single).unwrap();

        let graph = perceive_electron_counts(&mg).unwrap();
        let rings = perceive_rings(&graph);
        assert_eq!(rings.0.len(), 1);
    }

    #[test]
    fn aromaticity_flags_benzene() {
        let mut mg = MolecularGraph::new();
        let mut atoms = vec![];
        for _ in 0..6 {
            atoms.push(mg.add_atom(Element::C, 0));
        }
        for i in 0..6 {
            mg.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let rings = perceive_rings(&graph);
        apply_ring_annotations(&mut graph, &rings);
        perceive_generic_aromaticity(&mut graph, &rings).unwrap();

        assert!(graph.atoms.iter().all(|a| a.is_aromatic));
    }

    #[test]
    fn hybridization_assigns_resonant_for_aromatic_atoms() {
        let mut mg = MolecularGraph::new();
        let mut atoms = vec![];
        for _ in 0..6 {
            atoms.push(mg.add_atom(Element::C, 0));
        }
        for i in 0..6 {
            mg.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let rings = perceive_rings(&graph);
        apply_ring_annotations(&mut graph, &rings);
        perceive_generic_aromaticity(&mut graph, &rings).unwrap();
        perceive_generic_hybridization(&mut graph).unwrap();

        assert!(
            graph
                .atoms
                .iter()
                .all(|a| a.hybridization == Hybridization::Resonant)
        );
    }

    #[test]
    fn test_ammonium_ion_perception() {
        let mut mg = MolecularGraph::new();
        let n = mg.add_atom(Element::N, 1);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);

        mg.add_bond(n, h1, BondOrder::Single).unwrap();
        mg.add_bond(n, h2, BondOrder::Single).unwrap();
        mg.add_bond(n, h3, BondOrder::Single).unwrap();
        mg.add_bond(n, h4, BondOrder::Single).unwrap();

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let nitrogen = &graph.atoms[n];
        assert_eq!(nitrogen.formal_charge, 1);
        assert_eq!(nitrogen.lone_pairs, 0);
        assert_eq!(nitrogen.steric_number, 0);

        perceive_generic_hybridization(&mut graph).unwrap();
        let nitrogen = &graph.atoms[n];
        assert_eq!(nitrogen.steric_number, 4);
        assert_eq!(graph.atoms[n].hybridization, Hybridization::SP3);
    }

    #[test]
    fn test_acetate_ion_perception() {
        let mut mg = MolecularGraph::new();
        let c_methyl = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let c_carboxyl = mg.add_atom(Element::C, 0);
        let o_double = mg.add_atom(Element::O, 0);
        let o_single = mg.add_atom(Element::O, -1);

        mg.add_bond(c_methyl, h1, BondOrder::Single).unwrap();
        mg.add_bond(c_methyl, h2, BondOrder::Single).unwrap();
        mg.add_bond(c_methyl, h3, BondOrder::Single).unwrap();
        mg.add_bond(c_methyl, c_carboxyl, BondOrder::Single)
            .unwrap();
        mg.add_bond(c_carboxyl, o_double, BondOrder::Double)
            .unwrap();
        mg.add_bond(c_carboxyl, o_single, BondOrder::Single)
            .unwrap();

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let carbonyl_oxygen = &graph.atoms[o_double];
        assert_eq!(carbonyl_oxygen.formal_charge, 0);
        assert_eq!(carbonyl_oxygen.lone_pairs, 2);
        assert_eq!(carbonyl_oxygen.steric_number, 0);

        let alkoxide_oxygen = &graph.atoms[o_single];
        assert_eq!(alkoxide_oxygen.formal_charge, -1);
        assert_eq!(alkoxide_oxygen.lone_pairs, 3);
        assert_eq!(alkoxide_oxygen.steric_number, 0);

        templates::apply_functional_group_templates(&mut graph).unwrap();
        assert_eq!(graph.atoms[o_double].steric_number, 3);
        assert_eq!(graph.atoms[o_single].steric_number, 3);

        perceive_generic_hybridization(&mut graph).unwrap();
        assert_eq!(graph.atoms[o_double].hybridization, Hybridization::SP2);
        assert_eq!(graph.atoms[o_single].hybridization, Hybridization::SP2);
    }

    #[test]
    fn test_boron_trifluoride_hybridization() {
        let mut mg = MolecularGraph::new();
        let b = mg.add_atom(Element::B, 0);
        let f1 = mg.add_atom(Element::F, 0);
        let f2 = mg.add_atom(Element::F, 0);
        let f3 = mg.add_atom(Element::F, 0);

        mg.add_bond(b, f1, BondOrder::Single).unwrap();
        mg.add_bond(b, f2, BondOrder::Single).unwrap();
        mg.add_bond(b, f3, BondOrder::Single).unwrap();

        let mut graph = perceive_electron_counts(&mg).unwrap();
        perceive_generic_hybridization(&mut graph).unwrap();
        assert_eq!(graph.atoms[b].steric_number, 3);
        assert_eq!(graph.atoms[b].hybridization, Hybridization::SP2);
    }

    #[test]
    fn test_pyrrole_aromatic_pi_contribution() {
        let mut mg = MolecularGraph::new();
        let n = mg.add_atom(Element::N, 0);
        let c1 = mg.add_atom(Element::C, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let c3 = mg.add_atom(Element::C, 0);
        let c4 = mg.add_atom(Element::C, 0);
        let h_n = mg.add_atom(Element::H, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);

        mg.add_bond(n, c1, BondOrder::Single).unwrap();
        mg.add_bond(c1, c2, BondOrder::Double).unwrap();
        mg.add_bond(c2, c3, BondOrder::Single).unwrap();
        mg.add_bond(c3, c4, BondOrder::Double).unwrap();
        mg.add_bond(c4, n, BondOrder::Single).unwrap();
        mg.add_bond(n, h_n, BondOrder::Single).unwrap();
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();
        mg.add_bond(c2, h2, BondOrder::Single).unwrap();
        mg.add_bond(c3, h3, BondOrder::Single).unwrap();
        mg.add_bond(c4, h4, BondOrder::Single).unwrap();

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let rings = perceive_rings(&graph);
        apply_ring_annotations(&mut graph, &rings);
        perceive_generic_properties(&mut graph, &rings).unwrap();

        let ring_atoms = rings.0.iter().next().unwrap();
        let ring_set: HashSet<usize> = ring_atoms.iter().copied().collect();

        assert!(ring_atoms.iter().all(|&idx| graph.atoms[idx].is_aromatic));
        assert_eq!(count_atom_pi_contribution(n, &ring_set, &graph), Some(2));
    }

    #[test]
    fn test_furan_aromatic_pi_contribution() {
        let mut mg = MolecularGraph::new();
        let o = mg.add_atom(Element::O, 0);
        let c1 = mg.add_atom(Element::C, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let c3 = mg.add_atom(Element::C, 0);
        let c4 = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);

        mg.add_bond(o, c1, BondOrder::Single).unwrap();
        mg.add_bond(c1, c2, BondOrder::Double).unwrap();
        mg.add_bond(c2, c3, BondOrder::Single).unwrap();
        mg.add_bond(c3, c4, BondOrder::Double).unwrap();
        mg.add_bond(c4, o, BondOrder::Single).unwrap();
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();
        mg.add_bond(c2, h2, BondOrder::Single).unwrap();
        mg.add_bond(c3, h3, BondOrder::Single).unwrap();
        mg.add_bond(c4, h4, BondOrder::Single).unwrap();

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let rings = perceive_rings(&graph);
        apply_ring_annotations(&mut graph, &rings);
        perceive_generic_properties(&mut graph, &rings).unwrap();

        let ring_atoms = rings.0.iter().next().unwrap();
        let ring_set: HashSet<usize> = ring_atoms.iter().copied().collect();

        assert!(ring_atoms.iter().all(|&idx| graph.atoms[idx].is_aromatic));
        assert_eq!(count_atom_pi_contribution(o, &ring_set, &graph), Some(2));
    }

    #[test]
    fn test_pyridine_aromatic_pi_contribution() {
        let mut mg = MolecularGraph::new();
        let n = mg.add_atom(Element::N, 0);
        let c1 = mg.add_atom(Element::C, 0);
        let c2 = mg.add_atom(Element::C, 0);
        let c3 = mg.add_atom(Element::C, 0);
        let c4 = mg.add_atom(Element::C, 0);
        let c5 = mg.add_atom(Element::C, 0);
        let h1 = mg.add_atom(Element::H, 0);
        let h2 = mg.add_atom(Element::H, 0);
        let h3 = mg.add_atom(Element::H, 0);
        let h4 = mg.add_atom(Element::H, 0);
        let h5 = mg.add_atom(Element::H, 0);

        mg.add_bond(n, c1, BondOrder::Double).unwrap();
        mg.add_bond(c1, c2, BondOrder::Single).unwrap();
        mg.add_bond(c2, c3, BondOrder::Double).unwrap();
        mg.add_bond(c3, c4, BondOrder::Single).unwrap();
        mg.add_bond(c4, c5, BondOrder::Double).unwrap();
        mg.add_bond(c5, n, BondOrder::Single).unwrap();
        mg.add_bond(c1, h1, BondOrder::Single).unwrap();
        mg.add_bond(c2, h2, BondOrder::Single).unwrap();
        mg.add_bond(c3, h3, BondOrder::Single).unwrap();
        mg.add_bond(c4, h4, BondOrder::Single).unwrap();
        mg.add_bond(c5, h5, BondOrder::Single).unwrap();

        let mut graph = perceive_electron_counts(&mg).unwrap();
        let rings = perceive_rings(&graph);
        apply_ring_annotations(&mut graph, &rings);
        perceive_generic_properties(&mut graph, &rings).unwrap();

        let ring_atoms = rings.0.iter().next().unwrap();
        let ring_set: HashSet<usize> = ring_atoms.iter().copied().collect();

        assert!(ring_atoms.iter().all(|&idx| graph.atoms[idx].is_aromatic));
        assert_eq!(count_atom_pi_contribution(n, &ring_set, &graph), Some(1));
    }
}
