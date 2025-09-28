use std::fmt;
use std::str::FromStr;

pub mod graph;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Element {
    // --- Non-metals ---
    H = 1,
    He,
    B = 5,
    C,
    N,
    O,
    F,
    Ne,
    Si = 14,
    P,
    S,
    Cl,
    Ar,
    Ge = 32,
    As,
    Se,
    Br,
    Kr,
    Sb = 51,
    Te,
    I,
    Xe,
    At = 85,
    Rn,

    // --- Alkali Metals ---
    Li = 3,
    Na = 11,
    K = 19,
    Rb = 37,
    Cs = 55,
    Fr = 87,

    // --- Alkaline Earth Metals ---
    Be = 4,
    Mg = 12,
    Ca = 20,
    Sr = 38,
    Ba = 56,
    Ra = 88,

    // --- Transition Metals ---
    Sc = 21,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Y = 39,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    Hf = 72,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,

    // --- Post-transition Metals ---
    Al = 13,
    Ga = 31,
    In = 49,
    Tl = 81,
    Sn = 50,
    Pb = 82,
    Bi = 83,
    Po = 84,

    // --- Lanthanides ---
    La = 57,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu = 71,

    // --- Actinides ---
    Ac = 89,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr = 103,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseElementError {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum BondOrder {
    Single = 1,
    Double = 2,
    Triple = 3,
    Aromatic = 4,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseBondOrderError {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Hybridization {
    SP,
    SP2,
    SP3,
    Resonant,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseHybridizationError {
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

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SP" => Ok(Hybridization::SP),
            "SP2" => Ok(Hybridization::SP2),
            "SP3" => Ok(Hybridization::SP3),
            "Resonant" => Ok(Hybridization::Resonant),
            _ => Err(ParseHybridizationError {
                invalid_string: s.to_owned(),
            }),
        }
    }
}

impl fmt::Display for Hybridization {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let hyb_str = match self {
            Hybridization::SP => "SP",
            Hybridization::SP2 => "SP2",
            Hybridization::SP3 => "SP3",
            Hybridization::Resonant => "Resonant",
        };
        write!(f, "{}", hyb_str)
    }
}
