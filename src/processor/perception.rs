use super::graph::{PerceptionSource, ProcessingGraph, RingInfo};
use crate::core::error::{AnnotationError, TyperError};
use crate::core::graph::MolecularGraph;
use crate::core::{BondOrder, Element, Hybridization};
use std::collections::{HashSet, VecDeque};

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

pub(crate) fn perceive_rings(graph: &ProcessingGraph) -> RingInfo {
    if graph.atoms.is_empty() {
        return RingInfo::default();
    }

    let mut finder = JohnsonCycleFinder::new(graph);
    let sorted_vec_cycles = finder.find_cycles_internal();
    RingInfo(sorted_vec_cycles)
}

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

pub(crate) fn perceive_generic_properties(
    graph: &mut ProcessingGraph,
    ring_info: &RingInfo,
) -> Result<(), AnnotationError> {
    perceive_generic_aromaticity(graph, ring_info)?;
    perceive_generic_hybridization(graph)?;
    Ok(())
}

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

pub(crate) fn perceive_generic_hybridization(
    graph: &mut ProcessingGraph,
) -> Result<(), AnnotationError> {
    let mut initial_hybs = Vec::with_capacity(graph.atoms.len());
    for i in 0..graph.atoms.len() {
        if graph.atoms[i].perception_source == Some(PerceptionSource::Template) {
            initial_hybs.push(graph.atoms[i].hybridization);
            continue;
        }

        let atom = &mut graph.atoms[i];
        let steric = atom.degree + atom.lone_pairs;
        atom.steric_number = steric;

        let hyb = if atom.is_aromatic {
            Hybridization::Resonant
        } else if is_special_non_hybridized(atom.element) || atom.degree == 0 {
            Hybridization::None
        } else {
            match steric {
                4 => Hybridization::SP3,
                3 => Hybridization::SP2,
                2 => Hybridization::SP,
                _ => Hybridization::Unknown,
            }
        };
        initial_hybs.push(hyb);
    }

    for i in 0..graph.atoms.len() {
        if graph.atoms[i].perception_source == Some(PerceptionSource::Template) {
            continue;
        }

        let mut final_hyb = initial_hybs[i];

        if final_hyb == Hybridization::SP3 && graph.atoms[i].lone_pairs > 0 {
            let atom_view = &graph.atoms[i];
            let hydrogen_neighbors = graph.adjacency[atom_view.id]
                .iter()
                .filter(|(neighbor_id, _)| graph.atoms[*neighbor_id].element == Element::H)
                .count();
            let has_protective_hydrogen = atom_view.element == Element::O && hydrogen_neighbors > 0;
            let is_multihydrogen = hydrogen_neighbors >= 2;
            let is_anion = atom_view.formal_charge < 0;

            if !has_protective_hydrogen && !is_anion && !is_multihydrogen {
                let has_pi_system_neighbor =
                    graph.adjacency[atom_view.id]
                        .iter()
                        .any(|(neighbor_id, _)| {
                            let neighbor_hyb = initial_hybs[*neighbor_id];
                            matches!(
                                neighbor_hyb,
                                Hybridization::SP2 | Hybridization::SP | Hybridization::Resonant
                            )
                        });

                if has_pi_system_neighbor {
                    final_hyb = Hybridization::SP2;
                }
            }
        }

        if final_hyb == Hybridization::Unknown {
            return Err(AnnotationError::HybridizationInference {
                atom_id: graph.atoms[i].id,
            });
        }

        graph.atoms[i].hybridization = final_hyb;
        graph.atoms[i].perception_source = Some(PerceptionSource::Generic);

        graph.atoms[i].steric_number = match final_hyb {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            Hybridization::None => graph.atoms[i].degree + graph.atoms[i].lone_pairs,
            Hybridization::Unknown => graph.atoms[i].steric_number,
        };
    }

    Ok(())
}

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

    (total_pi_electrons - 2) % 4 == 0
}

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

fn bond_order_contribution(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 1,
    }
}

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

struct JohnsonCycleFinder<'a> {
    graph: &'a ProcessingGraph,
    all_cycles: HashSet<Vec<usize>>,
}

impl<'a> JohnsonCycleFinder<'a> {
    fn new(graph: &'a ProcessingGraph) -> Self {
        Self {
            graph,
            all_cycles: HashSet::new(),
        }
    }

    fn find_cycles_internal(&mut self) -> HashSet<Vec<usize>> {
        let num_atoms = self.graph.atoms.len();
        for i in 0..num_atoms {
            self.find_cycles_from_node(i);
        }
        self.all_cycles.clone()
    }

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
