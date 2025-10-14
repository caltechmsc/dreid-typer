//! Rule-based typing system for DREIDING force field atom types.
//!
//! This module implements the core rule engine that assigns DREIDING atom types
//! during the typing phase of the dreid-typer pipeline. It provides structures
//! for defining typing rules in TOML format and functions for parsing and
//! applying these rules to molecular graphs.

use crate::core::{Element, Hybridization, error::TyperError};
use serde::{Deserialize, Deserializer, de};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

mod default;

/// Deserializes a hash map with string keys that need to be parsed into a specific type.
///
/// This deserializer converts TOML tables with string keys into `HashMap<K, V>`
/// where the keys are parsed from strings using `FromStr`.
fn hash_map_from_str_keys<'de, K, V, D>(deserializer: D) -> Result<HashMap<K, V>, D::Error>
where
    K: FromStr + Eq + std::hash::Hash,
    K::Err: fmt::Display,
    V: Deserialize<'de>,
    D: Deserializer<'de>,
{
    let string_map = HashMap::<String, V>::deserialize(deserializer)?;
    string_map
        .into_iter()
        .map(|(k, v)| {
            K::from_str(&k)
                .map(|key| (key, v))
                .map_err(de::Error::custom)
        })
        .collect()
}

/// Conditions that must be met for a typing rule to apply to an atom.
///
/// This struct defines all the chemical and structural criteria that an atom
/// must satisfy for a DREIDING typing rule to be applicable. It supports both
/// simple property checks and complex neighbor-based conditions.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Default)]
#[serde(deny_unknown_fields)]
pub struct Conditions {
    /// The chemical element of the atom must match this value.
    #[serde(default, deserialize_with = "from_str_optional")]
    pub element: Option<Element>,
    /// The formal charge of the atom must match this value.
    #[serde(default)]
    pub formal_charge: Option<i8>,
    /// The degree (number of bonded neighbors) must match this value.
    #[serde(default)]
    pub degree: Option<u8>,
    /// The number of lone pairs on the atom must match this value.
    #[serde(default)]
    pub lone_pairs: Option<u8>,
    /// The steric number (degree + lone pairs) must match this value.
    #[serde(default)]
    pub steric_number: Option<u8>,
    /// The hybridization state must match this value.
    #[serde(default, deserialize_with = "from_str_optional")]
    pub hybridization: Option<Hybridization>,
    /// Whether the atom must be part of a ring structure.
    #[serde(default)]
    pub is_in_ring: Option<bool>,
    /// Whether the atom must be aromatic.
    #[serde(default)]
    pub is_aromatic: Option<bool>,
    /// The size of the smallest ring containing this atom must match this value.
    #[serde(default)]
    pub smallest_ring_size: Option<u8>,

    /// Counts of neighboring atoms by element type.
    #[serde(default, deserialize_with = "hash_map_from_str_keys")]
    pub neighbor_elements: HashMap<Element, u8>,
    /// Counts of neighboring atoms by their assigned DREIDING type.
    #[serde(default)]
    pub neighbor_types: HashMap<String, u8>,
}

/// Deserializes an optional string into an optional type that implements `FromStr`.
///
/// This handles TOML fields that may be absent or present as strings to be parsed.
fn from_str_optional<'de, T, D>(deserializer: D) -> Result<Option<T>, D::Error>
where
    T: FromStr,
    T::Err: fmt::Display,
    D: Deserializer<'de>,
{
    let opt_s = Option::<String>::deserialize(deserializer)?;
    match opt_s {
        Some(s) => T::from_str(&s).map(Some).map_err(de::Error::custom),
        None => Ok(None),
    }
}

/// A single rule for assigning DREIDING atom types.
///
/// Each rule defines a set of conditions that an atom must meet and the
/// resulting DREIDING type to assign when those conditions are satisfied.
/// Rules are prioritized and applied iteratively during the typing phase.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    /// A descriptive name for the rule (for debugging and documentation).
    pub name: String,
    /// The priority of this rule (higher values take precedence).
    pub priority: i32,
    /// The DREIDING atom type to assign when conditions are met.
    #[serde(rename = "type")]
    pub result_type: String,
    /// The conditions that must be satisfied for this rule to apply.
    pub conditions: Conditions,
}

/// Internal representation of a complete ruleset loaded from TOML.
#[derive(Debug, Clone, Deserialize)]
struct Ruleset {
    /// The list of rules in the ruleset.
    #[serde(rename = "rule")]
    rules: Vec<Rule>,
}

use std::sync::LazyLock;

/// Lazily-loaded default DREIDING typing rules.
///
/// The default rules are compiled into the binary and parsed once on first access.
/// This ensures efficient loading while providing comprehensive coverage of
/// common chemical environments.
static DEFAULT_RULES: LazyLock<Result<Vec<Rule>, TyperError>> =
    LazyLock::new(|| parse_rules(default::DEFAULT_RULES_TOML));

/// Parses a TOML string into a vector of typing rules.
///
/// This function deserializes TOML content containing rule definitions
/// and validates the structure and field values.
///
/// # Arguments
///
/// * `content` - The TOML string containing rule definitions.
///
/// # Returns
///
/// A `Result` containing the parsed rules on success.
///
/// # Errors
///
/// Returns `TyperError::RuleParse` if the TOML is malformed or contains invalid values.
///
/// # Examples
///
/// ```
/// use dreid_typer::rules::parse_rules;
///
/// let toml_content = r#"
///     [[rule]]
///     name = "Carbon_SP3"
///     priority = 100
///     type = "C_3"
///     conditions = { element = "C", hybridization = "SP3" }
/// "#;
/// let rules = parse_rules(toml_content).unwrap();
/// assert_eq!(rules.len(), 1);
/// assert_eq!(rules[0].name, "Carbon_SP3");
/// ```
pub fn parse_rules(content: &str) -> Result<Vec<Rule>, TyperError> {
    let ruleset: Ruleset =
        toml::from_str(content).map_err(|e| TyperError::RuleParse(e.to_string()))?;
    Ok(ruleset.rules)
}

/// Returns the default DREIDING typing rules compiled into the library.
///
/// These rules provide comprehensive coverage for common organic and
/// inorganic molecules. The rules are loaded lazily on first access.
pub(crate) fn get_default_rules() -> Result<&'static [Rule], TyperError> {
    DEFAULT_RULES
        .as_ref()
        .map(|vec| vec.as_slice())
        .map_err(|e| e.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn parse_rules_with_valid_content() {
        let toml_content = r#"
            [[rule]]
            name = "C_Aromatic"
            priority = 150
            type = "C_R"
            conditions = { element = "C", is_aromatic = true }
        "#;
        let rules = parse_rules(toml_content).unwrap();
        assert_eq!(rules.len(), 1);
        assert_eq!(
            rules[0],
            Rule {
                name: "C_Aromatic".to_string(),
                priority: 150,
                result_type: "C_R".to_string(),
                conditions: Conditions {
                    element: Some(Element::C),
                    is_aromatic: Some(true),
                    ..Default::default()
                }
            }
        );
    }

    #[test]
    fn parse_rules_with_invalid_toml_syntax() {
        let toml_content = "this is not valid toml";
        let result = parse_rules(toml_content);
        assert!(matches!(result, Err(TyperError::RuleParse(_))));
    }

    #[test]
    fn parse_rules_with_unknown_field_in_rule() {
        let toml_content = r#"
            [[rule]]
            name = "test"
            priority = 1
            type = "t"
            unknown_field = "value"
            conditions = { element = "C" }
        "#;
        let result = parse_rules(toml_content);
        assert!(matches!(result, Err(TyperError::RuleParse(_))));
    }

    #[test]
    fn parse_rules_with_unknown_field_in_conditions() {
        let toml_content = r#"
            [[rule]]
            name = "test"
            priority = 1
            type = "t"
            conditions = { element = "C", unknown_condition = true }
        "#;
        let result = parse_rules(toml_content);
        assert!(matches!(result, Err(TyperError::RuleParse(_))));
    }

    #[test]
    fn parse_rules_with_all_condition_fields() {
        let toml_content = r#"
            [[rule]]
            name = "test_all"
            priority = 100
            type = "X_0"
            [rule.conditions]
            element = "C"
            formal_charge = 0
            degree = 4
            lone_pairs = 0
            steric_number = 4
            hybridization = "SP3"
            is_in_ring = true
            is_aromatic = false
            smallest_ring_size = 6
            neighbor_elements = { H = 3, C = 1 }
            neighbor_types = { "C_3" = 1 }
        "#;
        let rules = parse_rules(toml_content).unwrap();
        let mut expected_neighbor_elements = HashMap::new();
        expected_neighbor_elements.insert(Element::H, 3);
        expected_neighbor_elements.insert(Element::C, 1);
        let mut expected_neighbor_types = HashMap::new();
        expected_neighbor_types.insert("C_3".to_string(), 1);

        assert_eq!(
            rules[0].conditions,
            Conditions {
                element: Some(Element::C),
                formal_charge: Some(0),
                degree: Some(4),
                lone_pairs: Some(0),
                steric_number: Some(4),
                hybridization: Some(Hybridization::SP3),
                is_in_ring: Some(true),
                is_aromatic: Some(false),
                smallest_ring_size: Some(6),
                neighbor_elements: expected_neighbor_elements,
                neighbor_types: expected_neighbor_types,
            }
        );
    }

    #[test]
    fn parse_rules_with_invalid_element() {
        let toml_content = r#"
            [[rule]]
            name = "test"
            priority = 1
            type = "t"
            conditions = { element = "InvalidElement" }
        "#;
        let result = parse_rules(toml_content);
        assert!(matches!(result, Err(TyperError::RuleParse(_))));
    }

    #[test]
    fn parse_rules_with_invalid_hybridization() {
        let toml_content = r#"
            [[rule]]
            name = "test"
            priority = 1
            type = "t"
            conditions = { hybridization = "SP4" }
        "#;
        let result = parse_rules(toml_content);
        assert!(matches!(result, Err(TyperError::RuleParse(_))));
    }

    #[test]
    fn parse_rules_with_invalid_neighbor_element_key() {
        let toml_content = r#"
            [[rule]]
            name = "test"
            priority = 1
            type = "t"
            conditions = { neighbor_elements = { "Xx" = 1 } }
        "#;
        let result = parse_rules(toml_content);
        assert!(matches!(result, Err(TyperError::RuleParse(_))));
    }

    #[test]
    fn get_default_rules_succeeds_and_is_not_empty() {
        let rules = get_default_rules().unwrap();
        assert!(!rules.is_empty());
        assert_eq!(rules.len(), 43);
    }
}
