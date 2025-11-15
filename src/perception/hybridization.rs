use super::model::{AnnotatedAtom, AnnotatedMolecule};
use crate::core::error::PerceptionError;
use crate::core::properties::{Element, Hybridization};

pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    for i in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[i];

        let steric_number = atom.degree + atom.lone_pairs;

        let hybridization = determine_hybridization(atom, steric_number)?;

        let atom_mut = &mut molecule.atoms[i];
        atom_mut.hybridization = hybridization;

        atom_mut.steric_number = match hybridization {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            _ => steric_number,
        };
    }
    Ok(())
}

fn determine_hybridization(
    atom: &AnnotatedAtom,
    steric_number: u8,
) -> Result<Hybridization, PerceptionError> {
    if is_non_hybridized_element(atom.element) {
        return Ok(Hybridization::None);
    }

    if atom.is_in_conjugated_system && !atom.is_anti_aromatic {
        if steric_number <= 3 || (steric_number == 4 && atom.lone_pairs > 0) {
            return Ok(Hybridization::Resonant);
        }
    }

    if atom.is_aromatic {
        return Ok(Hybridization::SP2);
    }

    match steric_number {
        4 => Ok(Hybridization::SP3),
        3 => Ok(Hybridization::SP2),
        2 => Ok(Hybridization::SP),
        0 | 1 => Ok(Hybridization::None),
        _ => Err(PerceptionError::HybridizationInference { atom_id: atom.id }),
    }
}

fn is_non_hybridized_element(element: Element) -> bool {
    matches!(
        element,
        Element::Li
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
            | Element::F
            | Element::Cl
            | Element::Br
            | Element::I
            | Element::At
            | Element::He
            | Element::Ne
            | Element::Ar
            | Element::Kr
            | Element::Xe
            | Element::Rn
            | Element::Sc
            | Element::Ti
            | Element::V
            | Element::Cr
            | Element::Mn
            | Element::Fe
            | Element::Co
            | Element::Ni
            | Element::Cu
            | Element::Zn
            | Element::Y
            | Element::Zr
            | Element::Nb
            | Element::Mo
            | Element::Tc
            | Element::Ru
            | Element::Rh
            | Element::Pd
            | Element::Ag
            | Element::Cd
            | Element::Hf
            | Element::Ta
            | Element::W
            | Element::Re
            | Element::Os
            | Element::Ir
            | Element::Pt
            | Element::Au
            | Element::Hg
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::properties::{BondOrder, Element, Hybridization};

    #[derive(Clone, Copy)]
    struct AtomSpec {
        element: Element,
        degree: u8,
        lone_pairs: u8,
        is_conjugated: bool,
        is_aromatic: bool,
        is_anti_aromatic: bool,
    }

    impl AtomSpec {
        fn new(element: Element) -> Self {
            Self {
                element,
                degree: 0,
                lone_pairs: 0,
                is_conjugated: false,
                is_aromatic: false,
                is_anti_aromatic: false,
            }
        }

        fn with_degree(mut self, degree: u8) -> Self {
            self.degree = degree;
            self
        }

        fn with_lone_pairs(mut self, lone_pairs: u8) -> Self {
            self.lone_pairs = lone_pairs;
            self
        }

        fn conjugated(mut self) -> Self {
            self.is_conjugated = true;
            self
        }

        fn aromatic(mut self) -> Self {
            self.is_aromatic = true;
            self
        }

        fn anti_aromatic(mut self) -> Self {
            self.is_anti_aromatic = true;
            self
        }
    }

    fn build_molecule(specs: &[AtomSpec]) -> AnnotatedMolecule {
        let atoms = specs
            .iter()
            .enumerate()
            .map(|(idx, spec)| AnnotatedAtom {
                id: idx,
                element: spec.element,
                formal_charge: 0,
                lone_pairs: spec.lone_pairs,
                degree: spec.degree,
                is_in_ring: spec.is_aromatic,
                smallest_ring_size: None,
                is_aromatic: spec.is_aromatic,
                is_anti_aromatic: spec.is_anti_aromatic,
                is_in_conjugated_system: spec.is_conjugated,
                is_resonant: false,
                steric_number: spec.degree + spec.lone_pairs,
                hybridization: Hybridization::Unknown,
            })
            .collect();

        let mut adjacency: Vec<Vec<(usize, BondOrder)>> = Vec::with_capacity(specs.len());
        for _ in 0..specs.len() {
            adjacency.push(Vec::new());
        }

        AnnotatedMolecule {
            atoms,
            bonds: Vec::new(),
            adjacency,
            rings: Vec::new(),
        }
    }

    fn perceive_specs(specs: &[AtomSpec]) -> AnnotatedMolecule {
        let mut molecule = build_molecule(specs);
        perceive(&mut molecule).expect("hybridization inference should succeed");
        molecule
    }

    #[test]
    fn non_hybridized_elements_remain_none() {
        let molecule = perceive_specs(&[AtomSpec::new(Element::Na).with_degree(1)]);
        let atom = &molecule.atoms[0];
        assert_eq!(atom.hybridization, Hybridization::None);
        assert_eq!(atom.steric_number, 1);
    }

    #[test]
    fn conjugated_atoms_with_steric_three_become_resonant() {
        let molecule = perceive_specs(&[AtomSpec::new(Element::C).with_degree(3).conjugated()]);
        let atom = &molecule.atoms[0];
        assert_eq!(atom.hybridization, Hybridization::Resonant);
        assert_eq!(atom.steric_number, 3);
    }

    #[test]
    fn conjugated_atoms_with_lone_pair_rehybridize_from_four() {
        let molecule = perceive_specs(&[AtomSpec::new(Element::N)
            .with_degree(3)
            .with_lone_pairs(1)
            .conjugated()]);
        let atom = &molecule.atoms[0];
        assert_eq!(atom.hybridization, Hybridization::Resonant);
        assert_eq!(atom.steric_number, 3, "steric number should collapse to 3");
    }

    #[test]
    fn anti_aromatic_atoms_skip_resonant_assignment() {
        let molecule = perceive_specs(&[AtomSpec::new(Element::C)
            .with_degree(3)
            .conjugated()
            .anti_aromatic()]);
        let atom = &molecule.atoms[0];
        assert_eq!(atom.hybridization, Hybridization::SP2);
        assert_eq!(atom.steric_number, 3);
    }

    #[test]
    fn aromatic_atoms_default_to_sp2_planarity() {
        let molecule = perceive_specs(&[AtomSpec::new(Element::C).with_degree(4).aromatic()]);
        let atom = &molecule.atoms[0];
        assert_eq!(atom.hybridization, Hybridization::SP2);
        assert_eq!(atom.steric_number, 3);
    }

    #[test]
    fn vsepr_rules_assign_expected_hybridizations() {
        let molecule = perceive_specs(&[
            AtomSpec::new(Element::C).with_degree(4),
            AtomSpec::new(Element::C).with_degree(3),
            AtomSpec::new(Element::C).with_degree(2),
            AtomSpec::new(Element::H).with_degree(1),
        ]);

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP3);
        assert_eq!(molecule.atoms[0].steric_number, 4);

        assert_eq!(molecule.atoms[1].hybridization, Hybridization::SP2);
        assert_eq!(molecule.atoms[1].steric_number, 3);

        assert_eq!(molecule.atoms[2].hybridization, Hybridization::SP);
        assert_eq!(molecule.atoms[2].steric_number, 2);

        assert_eq!(molecule.atoms[3].hybridization, Hybridization::None);
        assert_eq!(molecule.atoms[3].steric_number, 1);
    }

    #[test]
    fn steric_numbers_above_four_raise_an_error() {
        let mut molecule = build_molecule(&[AtomSpec::new(Element::C).with_degree(5)]);
        let err = perceive(&mut molecule).expect_err("steric 5 should fail");

        match err {
            PerceptionError::HybridizationInference { atom_id } => assert_eq!(atom_id, 0),
            other => panic!("unexpected error returned: {other:?}"),
        }
    }
}
