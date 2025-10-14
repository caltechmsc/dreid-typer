//! Core types and enumerations for chemical elements and molecular properties.
//!
//! This module defines fundamental types used throughout the dreid-typer library,
//! including chemical elements, bond orders, and hybridization states. These types
//! provide the basic building blocks for representing molecular structures and
//! properties in the perception, typing, and building phases of the pipeline.

use std::fmt;
use std::str::FromStr;

pub mod error;
pub mod graph;

/// Represents a chemical element with its atomic number.
///
/// This enum covers all known chemical elements, organized by their periodic table
/// groups. Elements are represented by their standard atomic symbols and can be
/// parsed from strings or displayed as symbols.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Element {
    // --- Non-metals ---
    /// Hydrogen (atomic number 1)
    H = 1,
    /// Helium (atomic number 2)
    He,
    /// Boron (atomic number 5)
    B = 5,
    /// Carbon (atomic number 6)
    C,
    /// Nitrogen (atomic number 7)
    N,
    /// Oxygen (atomic number 8)
    O,
    /// Fluorine (atomic number 9)
    F,
    /// Neon (atomic number 10)
    Ne,
    /// Silicon (atomic number 14)
    Si = 14,
    /// Phosphorus (atomic number 15)
    P,
    /// Sulfur (atomic number 16)
    S,
    /// Chlorine (atomic number 17)
    Cl,
    /// Argon (atomic number 18)
    Ar,
    /// Germanium (atomic number 32)
    Ge = 32,
    /// Arsenic (atomic number 33)
    As,
    /// Selenium (atomic number 34)
    Se,
    /// Bromine (atomic number 35)
    Br,
    /// Krypton (atomic number 36)
    Kr,
    /// Antimony (atomic number 51)
    Sb = 51,
    /// Tellurium (atomic number 52)
    Te,
    /// Iodine (atomic number 53)
    I,
    /// Xenon (atomic number 54)
    Xe,
    /// Astatine (atomic number 85)
    At = 85,
    /// Radon (atomic number 86)
    Rn,

    // --- Alkali Metals ---
    /// Lithium (atomic number 3)
    Li = 3,
    /// Sodium (atomic number 11)
    Na = 11,
    /// Potassium (atomic number 19)
    K = 19,
    /// Rubidium (atomic number 37)
    Rb = 37,
    /// Caesium (atomic number 55)
    Cs = 55,
    /// Francium (atomic number 87)
    Fr = 87,

    // --- Alkaline Earth Metals ---
    /// Beryllium (atomic number 4)
    Be = 4,
    /// Magnesium (atomic number 12)
    Mg = 12,
    /// Calcium (atomic number 20)
    Ca = 20,
    /// Strontium (atomic number 38)
    Sr = 38,
    /// Barium (atomic number 56)
    Ba = 56,
    /// Radium (atomic number 88)
    Ra = 88,

    // --- Transition Metals ---
    /// Scandium (atomic number 21)
    Sc = 21,
    /// Titanium (atomic number 22)
    Ti,
    /// Vanadium (atomic number 23)
    V,
    /// Chromium (atomic number 24)
    Cr,
    /// Manganese (atomic number 25)
    Mn,
    /// Iron (atomic number 26)
    Fe,
    /// Cobalt (atomic number 27)
    Co,
    /// Nickel (atomic number 28)
    Ni,
    /// Copper (atomic number 29)
    Cu,
    /// Zinc (atomic number 30)
    Zn,
    /// Yttrium (atomic number 39)
    Y = 39,
    /// Zirconium (atomic number 40)
    Zr,
    /// Niobium (atomic number 41)
    Nb,
    /// Molybdenum (atomic number 42)
    Mo,
    /// Technetium (atomic number 43)
    Tc,
    /// Ruthenium (atomic number 44)
    Ru,
    /// Rhodium (atomic number 45)
    Rh,
    /// Palladium (atomic number 46)
    Pd,
    /// Silver (atomic number 47)
    Ag,
    /// Cadmium (atomic number 48)
    Cd,
    /// Hafnium (atomic number 72)
    Hf = 72,
    /// Tantalum (atomic number 73)
    Ta,
    /// Tungsten (atomic number 74)
    W,
    /// Rhenium (atomic number 75)
    Re,
    /// Osmium (atomic number 76)
    Os,
    /// Iridium (atomic number 77)
    Ir,
    /// Platinum (atomic number 78)
    Pt,
    /// Gold (atomic number 79)
    Au,
    /// Mercury (atomic number 80)
    Hg,

    // --- Post-transition Metals ---
    /// Aluminium (atomic number 13)
    Al = 13,
    /// Gallium (atomic number 31)
    Ga = 31,
    /// Indium (atomic number 49)
    In = 49,
    /// Thallium (atomic number 81)
    Tl = 81,
    /// Tin (atomic number 50)
    Sn = 50,
    /// Lead (atomic number 82)
    Pb = 82,
    /// Bismuth (atomic number 83)
    Bi = 83,
    /// Polonium (atomic number 84)
    Po = 84,

    // --- Lanthanides ---
    /// Lanthanum (atomic number 57)
    La = 57,
    /// Cerium (atomic number 58)
    Ce,
    /// Praseodymium (atomic number 59)
    Pr,
    /// Neodymium (atomic number 60)
    Nd,
    /// Promethium (atomic number 61)
    Pm,
    /// Samarium (atomic number 62)
    Sm,
    /// Europium (atomic number 63)
    Eu,
    /// Gadolinium (atomic number 64)
    Gd,
    /// Terbium (atomic number 65)
    Tb,
    /// Dysprosium (atomic number 66)
    Dy,
    /// Holmium (atomic number 67)
    Ho,
    /// Erbium (atomic number 68)
    Er,
    /// Thulium (atomic number 69)
    Tm,
    /// Ytterbium (atomic number 70)
    Yb,
    /// Lutetium (atomic number 71)
    Lu = 71,

    // --- Actinides ---
    /// Actinium (atomic number 89)
    Ac = 89,
    /// Thorium (atomic number 90)
    Th,
    /// Protactinium (atomic number 91)
    Pa,
    /// Uranium (atomic number 92)
    U,
    /// Neptunium (atomic number 93)
    Np,
    /// Plutonium (atomic number 94)
    Pu,
    /// Americium (atomic number 95)
    Am,
    /// Curium (atomic number 96)
    Cm,
    /// Berkelium (atomic number 97)
    Bk,
    /// Californium (atomic number 98)
    Cf,
    /// Einsteinium (atomic number 99)
    Es,
    /// Fermium (atomic number 100)
    Fm,
    /// Mendelevium (atomic number 101)
    Md,
    /// Nobelium (atomic number 102)
    No,
    /// Lawrencium (atomic number 103)
    Lr = 103,
}

/// Error type for failed element symbol parsing.
///
/// This error is returned when attempting to parse an invalid string
/// as a chemical element symbol.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseElementError {
    /// The invalid string that could not be parsed.
    invalid_string: String,
}

impl fmt::Display for ParseElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid element symbol: '{}'", self.invalid_string)
    }
}
impl std::error::Error for ParseElementError {}

impl FromStr for Element {
    type Err = ParseElementError;

    /// Parses a chemical element from its standard symbol string.
    ///
    /// # Arguments
    ///
    /// * `s` - The element symbol string (e.g., "C" for carbon, "Na" for sodium).
    ///
    /// # Returns
    ///
    /// A `Result` containing the parsed `Element` on success.
    ///
    /// # Errors
    ///
    /// Returns `ParseElementError` if the string is not a valid element symbol.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::str::FromStr;
    /// use dreid_typer::Element;
    ///
    /// let carbon = Element::from_str("C").unwrap();
    /// assert_eq!(carbon, Element::C);
    ///
    /// let sodium = Element::from_str("Na").unwrap();
    /// assert_eq!(sodium, Element::Na);
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Element::H),
            "He" => Ok(Element::He),
            "Li" => Ok(Element::Li),
            "Be" => Ok(Element::Be),
            "B" => Ok(Element::B),
            "C" => Ok(Element::C),
            "N" => Ok(Element::N),
            "O" => Ok(Element::O),
            "F" => Ok(Element::F),
            "Ne" => Ok(Element::Ne),
            "Na" => Ok(Element::Na),
            "Mg" => Ok(Element::Mg),
            "Al" => Ok(Element::Al),
            "Si" => Ok(Element::Si),
            "P" => Ok(Element::P),
            "S" => Ok(Element::S),
            "Cl" => Ok(Element::Cl),
            "Ar" => Ok(Element::Ar),
            "K" => Ok(Element::K),
            "Ca" => Ok(Element::Ca),
            "Sc" => Ok(Element::Sc),
            "Ti" => Ok(Element::Ti),
            "V" => Ok(Element::V),
            "Cr" => Ok(Element::Cr),
            "Mn" => Ok(Element::Mn),
            "Fe" => Ok(Element::Fe),
            "Co" => Ok(Element::Co),
            "Ni" => Ok(Element::Ni),
            "Cu" => Ok(Element::Cu),
            "Zn" => Ok(Element::Zn),
            "Ga" => Ok(Element::Ga),
            "Ge" => Ok(Element::Ge),
            "As" => Ok(Element::As),
            "Se" => Ok(Element::Se),
            "Br" => Ok(Element::Br),
            "Kr" => Ok(Element::Kr),
            "Rb" => Ok(Element::Rb),
            "Sr" => Ok(Element::Sr),
            "Y" => Ok(Element::Y),
            "Zr" => Ok(Element::Zr),
            "Nb" => Ok(Element::Nb),
            "Mo" => Ok(Element::Mo),
            "Tc" => Ok(Element::Tc),
            "Ru" => Ok(Element::Ru),
            "Rh" => Ok(Element::Rh),
            "Pd" => Ok(Element::Pd),
            "Ag" => Ok(Element::Ag),
            "Cd" => Ok(Element::Cd),
            "In" => Ok(Element::In),
            "Sn" => Ok(Element::Sn),
            "Sb" => Ok(Element::Sb),
            "Te" => Ok(Element::Te),
            "I" => Ok(Element::I),
            "Xe" => Ok(Element::Xe),
            "Cs" => Ok(Element::Cs),
            "Ba" => Ok(Element::Ba),
            "La" => Ok(Element::La),
            "Ce" => Ok(Element::Ce),
            "Pr" => Ok(Element::Pr),
            "Nd" => Ok(Element::Nd),
            "Pm" => Ok(Element::Pm),
            "Sm" => Ok(Element::Sm),
            "Eu" => Ok(Element::Eu),
            "Gd" => Ok(Element::Gd),
            "Tb" => Ok(Element::Tb),
            "Dy" => Ok(Element::Dy),
            "Ho" => Ok(Element::Ho),
            "Er" => Ok(Element::Er),
            "Tm" => Ok(Element::Tm),
            "Yb" => Ok(Element::Yb),
            "Lu" => Ok(Element::Lu),
            "Hf" => Ok(Element::Hf),
            "Ta" => Ok(Element::Ta),
            "W" => Ok(Element::W),
            "Re" => Ok(Element::Re),
            "Os" => Ok(Element::Os),
            "Ir" => Ok(Element::Ir),
            "Pt" => Ok(Element::Pt),
            "Au" => Ok(Element::Au),
            "Hg" => Ok(Element::Hg),
            "Tl" => Ok(Element::Tl),
            "Pb" => Ok(Element::Pb),
            "Bi" => Ok(Element::Bi),
            "Po" => Ok(Element::Po),
            "At" => Ok(Element::At),
            "Rn" => Ok(Element::Rn),
            "Fr" => Ok(Element::Fr),
            "Ra" => Ok(Element::Ra),
            "Ac" => Ok(Element::Ac),
            "Th" => Ok(Element::Th),
            "Pa" => Ok(Element::Pa),
            "U" => Ok(Element::U),
            "Np" => Ok(Element::Np),
            "Pu" => Ok(Element::Pu),
            "Am" => Ok(Element::Am),
            "Cm" => Ok(Element::Cm),
            "Bk" => Ok(Element::Bk),
            "Cf" => Ok(Element::Cf),
            "Es" => Ok(Element::Es),
            "Fm" => Ok(Element::Fm),
            "Md" => Ok(Element::Md),
            "No" => Ok(Element::No),
            "Lr" => Ok(Element::Lr),
            _ => Err(ParseElementError {
                invalid_string: s.to_owned(),
            }),
        }
    }
}

impl fmt::Display for Element {
    /// Displays the element as its standard chemical symbol.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::Element;
    ///
    /// let carbon = Element::C;
    /// assert_eq!(format!("{}", carbon), "C");
    ///
    /// let sodium = Element::Na;
    /// assert_eq!(format!("{}", sodium), "Na");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let symbol = match self {
            Element::H => "H",
            Element::He => "He",
            Element::Li => "Li",
            Element::Be => "Be",
            Element::B => "B",
            Element::C => "C",
            Element::N => "N",
            Element::O => "O",
            Element::F => "F",
            Element::Ne => "Ne",
            Element::Na => "Na",
            Element::Mg => "Mg",
            Element::Al => "Al",
            Element::Si => "Si",
            Element::P => "P",
            Element::S => "S",
            Element::Cl => "Cl",
            Element::Ar => "Ar",
            Element::K => "K",
            Element::Ca => "Ca",
            Element::Sc => "Sc",
            Element::Ti => "Ti",
            Element::V => "V",
            Element::Cr => "Cr",
            Element::Mn => "Mn",
            Element::Fe => "Fe",
            Element::Co => "Co",
            Element::Ni => "Ni",
            Element::Cu => "Cu",
            Element::Zn => "Zn",
            Element::Ga => "Ga",
            Element::Ge => "Ge",
            Element::As => "As",
            Element::Se => "Se",
            Element::Br => "Br",
            Element::Kr => "Kr",
            Element::Rb => "Rb",
            Element::Sr => "Sr",
            Element::Y => "Y",
            Element::Zr => "Zr",
            Element::Nb => "Nb",
            Element::Mo => "Mo",
            Element::Tc => "Tc",
            Element::Ru => "Ru",
            Element::Rh => "Rh",
            Element::Pd => "Pd",
            Element::Ag => "Ag",
            Element::Cd => "Cd",
            Element::In => "In",
            Element::Sn => "Sn",
            Element::Sb => "Sb",
            Element::Te => "Te",
            Element::I => "I",
            Element::Xe => "Xe",
            Element::Cs => "Cs",
            Element::Ba => "Ba",
            Element::La => "La",
            Element::Ce => "Ce",
            Element::Pr => "Pr",
            Element::Nd => "Nd",
            Element::Pm => "Pm",
            Element::Sm => "Sm",
            Element::Eu => "Eu",
            Element::Gd => "Gd",
            Element::Tb => "Tb",
            Element::Dy => "Dy",
            Element::Ho => "Ho",
            Element::Er => "Er",
            Element::Tm => "Tm",
            Element::Yb => "Yb",
            Element::Lu => "Lu",
            Element::Hf => "Hf",
            Element::Ta => "Ta",
            Element::W => "W",
            Element::Re => "Re",
            Element::Os => "Os",
            Element::Ir => "Ir",
            Element::Pt => "Pt",
            Element::Au => "Au",
            Element::Hg => "Hg",
            Element::Tl => "Tl",
            Element::Pb => "Pb",
            Element::Bi => "Bi",
            Element::Po => "Po",
            Element::At => "At",
            Element::Rn => "Rn",
            Element::Fr => "Fr",
            Element::Ra => "Ra",
            Element::Ac => "Ac",
            Element::Th => "Th",
            Element::Pa => "Pa",
            Element::U => "U",
            Element::Np => "Np",
            Element::Pu => "Pu",
            Element::Am => "Am",
            Element::Cm => "Cm",
            Element::Bk => "Bk",
            Element::Cf => "Cf",
            Element::Es => "Es",
            Element::Fm => "Fm",
            Element::Md => "Md",
            Element::No => "No",
            Element::Lr => "Lr",
        };
        write!(f, "{}", symbol)
    }
}

/// Represents the order of a chemical bond.
///
/// This enum defines the possible bond orders used in molecular graphs,
/// including aromatic bonds which are treated specially in force field calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum BondOrder {
    /// Single bond (order 1)
    Single = 1,
    /// Double bond (order 2)
    Double = 2,
    /// Triple bond (order 3)
    Triple = 3,
    /// Aromatic bond (special case for conjugated systems)
    Aromatic = 4,
}

/// Error type for failed bond order parsing.
///
/// This error is returned when attempting to parse an invalid string
/// as a bond order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseBondOrderError {
    /// The invalid string that could not be parsed.
    invalid_string: String,
}

impl fmt::Display for ParseBondOrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid bond order string: '{}'", self.invalid_string)
    }
}
impl std::error::Error for ParseBondOrderError {}

impl FromStr for BondOrder {
    type Err = ParseBondOrderError;

    /// Parses a bond order from its string representation.
    ///
    /// # Arguments
    ///
    /// * `s` - The bond order string ("Single", "Double", "Triple", or "Aromatic").
    ///
    /// # Returns
    ///
    /// A `Result` containing the parsed `BondOrder` on success.
    ///
    /// # Errors
    ///
    /// Returns `ParseBondOrderError` if the string is not a valid bond order.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::str::FromStr;
    /// use dreid_typer::BondOrder;
    ///
    /// let single = BondOrder::from_str("Single").unwrap();
    /// assert_eq!(single, BondOrder::Single);
    ///
    /// let double = BondOrder::from_str("Double").unwrap();
    /// assert_eq!(double, BondOrder::Double);
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Single" => Ok(BondOrder::Single),
            "Double" => Ok(BondOrder::Double),
            "Triple" => Ok(BondOrder::Triple),
            "Aromatic" => Ok(BondOrder::Aromatic),
            _ => Err(ParseBondOrderError {
                invalid_string: s.to_owned(),
            }),
        }
    }
}

impl fmt::Display for BondOrder {
    /// Displays the bond order as its string representation.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::BondOrder;
    ///
    /// let single = BondOrder::Single;
    /// assert_eq!(format!("{}", single), "Single");
    ///
    /// let aromatic = BondOrder::Aromatic;
    /// assert_eq!(format!("{}", aromatic), "Aromatic");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let order_str = match self {
            BondOrder::Single => "Single",
            BondOrder::Double => "Double",
            BondOrder::Triple => "Triple",
            BondOrder::Aromatic => "Aromatic",
        };
        write!(f, "{}", order_str)
    }
}

/// Represents the hybridization state of an atom.
///
/// Hybridization describes the orbital configuration of atoms in molecules,
/// which is crucial for determining molecular geometry and force field parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Hybridization {
    /// sp hybridization (linear geometry)
    SP,
    /// sp² hybridization (trigonal planar geometry)
    SP2,
    /// sp³ hybridization (tetrahedral geometry)
    SP3,
    /// Resonant hybridization (for aromatic systems)
    Resonant,
    /// No hybridization (for certain special cases)
    None,
    /// Unknown hybridization (fallback state)
    Unknown,
}

/// Error type for failed hybridization parsing.
///
/// This error is returned when attempting to parse an invalid string
/// as a hybridization state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseHybridizationError {
    /// The invalid string that could not be parsed.
    invalid_string: String,
}

impl fmt::Display for ParseHybridizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid hybridization string: '{}'", self.invalid_string)
    }
}
impl std::error::Error for ParseHybridizationError {}

impl FromStr for Hybridization {
    type Err = ParseHybridizationError;

    /// Parses a hybridization state from its string representation.
    ///
    /// # Arguments
    ///
    /// * `s` - The hybridization string ("SP", "SP2", "SP3", "Resonant", or "None").
    ///
    /// # Returns
    ///
    /// A `Result` containing the parsed `Hybridization` on success.
    ///
    /// # Errors
    ///
    /// Returns `ParseHybridizationError` if the string is not a valid hybridization.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::str::FromStr;
    /// use dreid_typer::Hybridization;
    ///
    /// let sp3 = Hybridization::from_str("SP3").unwrap();
    /// assert_eq!(sp3, Hybridization::SP3);
    ///
    /// let resonant = Hybridization::from_str("Resonant").unwrap();
    /// assert_eq!(resonant, Hybridization::Resonant);
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SP" => Ok(Hybridization::SP),
            "SP2" => Ok(Hybridization::SP2),
            "SP3" => Ok(Hybridization::SP3),
            "Resonant" => Ok(Hybridization::Resonant),
            "None" => Ok(Hybridization::None),
            _ => Err(ParseHybridizationError {
                invalid_string: s.to_owned(),
            }),
        }
    }
}

impl fmt::Display for Hybridization {
    /// Displays the hybridization state as its string representation.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::Hybridization;
    ///
    /// let sp3 = Hybridization::SP3;
    /// assert_eq!(format!("{}", sp3), "SP3");
    ///
    /// let resonant = Hybridization::Resonant;
    /// assert_eq!(format!("{}", resonant), "Resonant");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let hyb_str = match self {
            Hybridization::SP => "SP",
            Hybridization::SP2 => "SP2",
            Hybridization::SP3 => "SP3",
            Hybridization::Resonant => "Resonant",
            Hybridization::None => "None",
            Hybridization::Unknown => "Unknown",
        };
        write!(f, "{}", hyb_str)
    }
}
