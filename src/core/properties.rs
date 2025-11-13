use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u16)]
pub enum Element {
    // --- Non-metals ---
    H = 1,
    He = 2,
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

    // --- Post-transition Metals ---
    Al = 13,
    Ga = 31,
    In = 49,
    Sn = 50,
    Tl = 81,
    Pb,
    Bi,
    Po,

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
    Lu,

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
    Lr,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseElementError(String);

impl fmt::Display for ParseElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid element symbol: '{}'", self.0)
    }
}
impl std::error::Error for ParseElementError {}

impl FromStr for Element {
    type Err = ParseElementError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Self::H),
            "He" => Ok(Self::He),
            "Li" => Ok(Self::Li),
            "Be" => Ok(Self::Be),
            "B" => Ok(Self::B),
            "C" => Ok(Self::C),
            "N" => Ok(Self::N),
            "O" => Ok(Self::O),
            "F" => Ok(Self::F),
            "Ne" => Ok(Self::Ne),
            "Na" => Ok(Self::Na),
            "Mg" => Ok(Self::Mg),
            "Al" => Ok(Self::Al),
            "Si" => Ok(Self::Si),
            "P" => Ok(Self::P),
            "S" => Ok(Self::S),
            "Cl" => Ok(Self::Cl),
            "Ar" => Ok(Self::Ar),
            "K" => Ok(Self::K),
            "Ca" => Ok(Self::Ca),
            "Sc" => Ok(Self::Sc),
            "Ti" => Ok(Self::Ti),
            "V" => Ok(Self::V),
            "Cr" => Ok(Self::Cr),
            "Mn" => Ok(Self::Mn),
            "Fe" => Ok(Self::Fe),
            "Co" => Ok(Self::Co),
            "Ni" => Ok(Self::Ni),
            "Cu" => Ok(Self::Cu),
            "Zn" => Ok(Self::Zn),
            "Ga" => Ok(Self::Ga),
            "Ge" => Ok(Self::Ge),
            "As" => Ok(Self::As),
            "Se" => Ok(Self::Se),
            "Br" => Ok(Self::Br),
            "Kr" => Ok(Self::Kr),
            "Rb" => Ok(Self::Rb),
            "Sr" => Ok(Self::Sr),
            "Y" => Ok(Self::Y),
            "Zr" => Ok(Self::Zr),
            "Nb" => Ok(Self::Nb),
            "Mo" => Ok(Self::Mo),
            "Tc" => Ok(Self::Tc),
            "Ru" => Ok(Self::Ru),
            "Rh" => Ok(Self::Rh),
            "Pd" => Ok(Self::Pd),
            "Ag" => Ok(Self::Ag),
            "Cd" => Ok(Self::Cd),
            "In" => Ok(Self::In),
            "Sn" => Ok(Self::Sn),
            "Sb" => Ok(Self::Sb),
            "Te" => Ok(Self::Te),
            "I" => Ok(Self::I),
            "Xe" => Ok(Self::Xe),
            "Cs" => Ok(Self::Cs),
            "Ba" => Ok(Self::Ba),
            "La" => Ok(Self::La),
            "Ce" => Ok(Self::Ce),
            "Pr" => Ok(Self::Pr),
            "Nd" => Ok(Self::Nd),
            "Pm" => Ok(Self::Pm),
            "Sm" => Ok(Self::Sm),
            "Eu" => Ok(Self::Eu),
            "Gd" => Ok(Self::Gd),
            "Tb" => Ok(Self::Tb),
            "Dy" => Ok(Self::Dy),
            "Ho" => Ok(Self::Ho),
            "Er" => Ok(Self::Er),
            "Tm" => Ok(Self::Tm),
            "Yb" => Ok(Self::Yb),
            "Lu" => Ok(Self::Lu),
            "Hf" => Ok(Self::Hf),
            "Ta" => Ok(Self::Ta),
            "W" => Ok(Self::W),
            "Re" => Ok(Self::Re),
            "Os" => Ok(Self::Os),
            "Ir" => Ok(Self::Ir),
            "Pt" => Ok(Self::Pt),
            "Au" => Ok(Self::Au),
            "Hg" => Ok(Self::Hg),
            "Tl" => Ok(Self::Tl),
            "Pb" => Ok(Self::Pb),
            "Bi" => Ok(Self::Bi),
            "Po" => Ok(Self::Po),
            "At" => Ok(Self::At),
            "Rn" => Ok(Self::Rn),
            "Fr" => Ok(Self::Fr),
            "Ra" => Ok(Self::Ra),
            "Ac" => Ok(Self::Ac),
            "Th" => Ok(Self::Th),
            "Pa" => Ok(Self::Pa),
            "U" => Ok(Self::U),
            "Np" => Ok(Self::Np),
            "Pu" => Ok(Self::Pu),
            "Am" => Ok(Self::Am),
            "Cm" => Ok(Self::Cm),
            "Bk" => Ok(Self::Bk),
            "Cf" => Ok(Self::Cf),
            "Es" => Ok(Self::Es),
            "Fm" => Ok(Self::Fm),
            "Md" => Ok(Self::Md),
            "No" => Ok(Self::No),
            "Lr" => Ok(Self::Lr),
            _ => Err(ParseElementError(s.to_string())),
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
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
pub struct ParseBondOrderError(String);

impl fmt::Display for ParseBondOrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid bond order: '{}'", self.0)
    }
}
impl std::error::Error for ParseBondOrderError {}

impl FromStr for BondOrder {
    type Err = ParseBondOrderError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Single" => Ok(Self::Single),
            "Double" => Ok(Self::Double),
            "Triple" => Ok(Self::Triple),
            "Aromatic" => Ok(Self::Aromatic),
            _ => Err(ParseBondOrderError(s.to_string())),
        }
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
