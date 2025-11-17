//! Defines the serialized rule language that drives the DREIDING typing engine.
//!
//! This module owns the `Rule` schema, parsing helpers, cache for the built-in rule deck, and
//! serde utilities that allow rules to reference elements, hybridizations, and neighbor counts by
//! symbolic keys.

use crate::core::properties::{Element, Hybridization};
use serde::{Deserialize, Deserializer, de};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use std::sync::OnceLock;

/// User-facing rule that assigns a DREIDING type when its conditions match an atom.
///
/// Rules are ordered by priority and evaluated by the typing engine until the first rule whose
/// `conditions` filters match succeeds, at which point the `result_type` label is emitted.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    /// Human-readable identifier that makes debugging rule decks easier.
    pub name: String,
    /// Ordering value where larger priorities are evaluated before smaller ones.
    pub priority: i32,
    /// Resulting atom type string written into the topology output.
    #[serde(rename = "type")]
    pub result_type: String,
    /// Set of optional property filters that determine when the rule fires.
    pub conditions: Conditions,
}

/// Optional property filters attached to a [`Rule`].
///
/// Each field defaults to `None` or an empty map, meaning the corresponding property is not
/// considered when deciding whether a rule applies. When populated, the field describes an exact
/// value that must match the target atom.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Default)]
#[serde(deny_unknown_fields)]
pub struct Conditions {
    /// Required element identity.
    #[serde(default)]
    pub element: Option<Element>,
    /// Required formal charge assigned during perception.
    #[serde(default)]
    pub formal_charge: Option<i8>,

    /// Required degree (number of σ bonds) for the atom.
    #[serde(default)]
    pub degree: Option<u8>,
    /// Whether the atom must belong (or not belong) to a ring system.
    #[serde(default)]
    pub is_in_ring: Option<bool>,

    /// Required lone-pair count after electron perception.
    #[serde(default)]
    pub lone_pairs: Option<u8>,
    /// Required hybridization assignment (SP, SP2, SP3, etc.).
    #[serde(default)]
    pub hybridization: Option<Hybridization>,

    /// Whether the atom must be aromatic.
    #[serde(default)]
    pub is_aromatic: Option<bool>,
    /// Whether the atom must be anti-aromatic.
    #[serde(default)]
    pub is_anti_aromatic: Option<bool>,
    /// Whether the atom must participate in any resonance system.
    #[serde(default)]
    pub is_resonant: Option<bool>,

    /// Minimum counts for neighbor elements keyed by element symbol strings.
    #[serde(default, deserialize_with = "deserialize_str_keyed_map")]
    pub neighbor_elements: HashMap<Element, u8>,
    /// Minimum counts for neighbor atom types identified by their DREIDING labels.
    #[serde(default)]
    pub neighbor_types: HashMap<String, u8>,
}

/// Helper struct that mirrors the `[ [rule] ]` array in the TOML file.
#[derive(Deserialize)]
struct Ruleset {
    #[serde(rename = "rule")]
    rules: Vec<Rule>,
}

/// Parses a TOML ruleset string into a list of [`Rule`] values.
///
/// The input must contain a sequence of `[[rule]]` tables whose fields match the schema defined
/// by [`Rule`] and [`Conditions`]. Any unknown keys or invalid enum values will cause parsing to
/// fail.
///
/// # Arguments
///
/// * `content` - TOML source text containing zero or more `[[rule]]` tables.
///
/// # Returns
///
/// Vector of fully materialized [`Rule`] structs in the order they appear in the file.
///
/// # Errors
///
/// Returns [`toml::de::Error`] when the document is not valid TOML or when one of the rule fields
/// fails validation (e.g., unknown element symbol, missing required property, or wrong type).
///
/// # Examples
///
/// ```
/// use dreid_typer::rules::parse_rules;
///
/// let toml = r#"
///     [[rule]]
///     name = "C_sp2"
///     priority = 10
///     type = "C_R"
///     [rule.conditions]
///     element = "C"
///     degree = 3
/// "#;
///
/// let rules = parse_rules(toml).unwrap();
/// assert_eq!(rules.len(), 1);
/// assert_eq!(rules[0].name, "C_sp2");
/// assert_eq!(rules[0].result_type, "C_R");
/// ```
pub fn parse_rules(content: &str) -> Result<Vec<Rule>, toml::de::Error> {
    let ruleset: Ruleset = toml::from_str(content)?;
    Ok(ruleset.rules)
}

static DEFAULT_RULES: OnceLock<Vec<Rule>> = OnceLock::new();

const DEFAULT_RULES_TOML: &str = include_str!("../../resources/default.rules.toml");

/// Returns the lazily parsed, embedded DREIDING rule deck.
///
/// This function parses the bundled `resources/default.rules.toml` file the first time it is
/// called and caches the resulting `Vec<Rule>` behind a [`OnceLock`]. Subsequent calls reuse the
/// same slice reference.
///
/// # Returns
///
/// Immutable slice containing all default rules in priority order.
///
/// # Panics
///
/// Panics if the embedded TOML blob cannot be parsed. Such a failure indicates that the library
/// was built with an invalid `default.rules.toml` file.
pub fn get_default_rules() -> &'static [Rule] {
    DEFAULT_RULES.get_or_init(|| {
        parse_rules(DEFAULT_RULES_TOML)
            .expect("Failed to parse embedded default DREIDING rules. This is a library bug.")
    })
}

/// Generic helper that deserializes TOML maps with string keys into typed keys via `FromStr`.
///
/// # Arguments
///
/// * `deserializer` - Serde deserializer for the field annotated with `deserialize_with`.
///
/// # Errors
///
/// Returns a serde error when either the nested map cannot be read or any key fails to parse via
/// `FromStr`.
fn deserialize_str_keyed_map<'de, K, V, D>(deserializer: D) -> Result<HashMap<K, V>, D::Error>
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

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;
    use std::ptr;

    const SAMPLE_RULES: &str = r#"
        [[rule]]
        name = "C_sp2"
        priority = 10
        type = "C_R"
        [rule.conditions]
        element = "C"
        degree = 3
        hybridization = "SP2"
        is_aromatic = true
        neighbor_elements = { N = 1 }
        neighbor_types = { "N_R" = 1 }

        [[rule]]
        name = "H_sp"
        priority = 5
        type = "H_"
        [rule.conditions]
        element = "H"
    "#;

    #[test]
    fn parse_rules_parses_multiple_entries() {
        let rules = parse_rules(SAMPLE_RULES).expect("sample TOML should parse");
        assert_eq!(rules.len(), 2);

        let first = &rules[0];
        assert_eq!(first.name, "C_sp2");
        assert_eq!(first.priority, 10);
        assert_eq!(first.result_type, "C_R");
        assert_eq!(first.conditions.element, Some(Element::C));
        assert_eq!(first.conditions.degree, Some(3));
        assert_eq!(first.conditions.hybridization, Some(Hybridization::SP2));
        assert_eq!(first.conditions.is_aromatic, Some(true));
        assert_eq!(
            first.conditions.neighbor_elements.get(&Element::N),
            Some(&1)
        );
        assert_eq!(first.conditions.neighbor_types.get("N_R"), Some(&1));

        let second = &rules[1];
        assert_eq!(second.name, "H_sp");
        assert_eq!(second.conditions.element, Some(Element::H));
        assert!(second.conditions.neighbor_elements.is_empty());
    }

    #[test]
    fn parse_rules_rejects_missing_required_fields() {
        let invalid = r#"
            [[rule]]
            name = "Invalid"
            priority = 1
            [rule.conditions]
            element = "C"
        "#;
        let err = parse_rules(invalid).expect_err("rule missing 'type' must fail");
        assert!(err.to_string().contains("missing field"));
    }

    #[test]
    fn parse_rules_reports_invalid_neighbor_element_key() {
        let invalid = r#"
            [[rule]]
            name = "InvalidElement"
            priority = 1
            type = "C_R"
            [rule.conditions]
            neighbor_elements = { Xx = 1 }
        "#;
        let err = parse_rules(invalid).expect_err("unknown element symbol should fail");
        assert!(
            err.to_string().contains("Xx"),
            "error should mention the problematic key"
        );
    }

    #[test]
    fn get_default_rules_is_cached_and_non_empty() {
        let first = get_default_rules();
        let second = get_default_rules();
        assert!(!first.is_empty(), "embedded rules must not be empty");
        assert!(ptr::eq(first, second), "rules slice should be cached");
    }

    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    enum SmallEnum {
        Foo,
        Bar,
    }

    impl FromStr for SmallEnum {
        type Err = &'static str;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s {
                "Foo" => Ok(SmallEnum::Foo),
                "Bar" => Ok(SmallEnum::Bar),
                _ => Err("unknown key"),
            }
        }
    }

    #[derive(Debug, Deserialize)]
    struct SmallMap {
        #[serde(deserialize_with = "deserialize_str_keyed_map")]
        map: HashMap<SmallEnum, u8>,
    }

    #[test]
    fn deserialize_str_keyed_map_supports_custom_enums() {
        let value: SmallMap =
            toml::from_str("map = { Foo = 1, Bar = 2 }").expect("enum keys should deserialize");
        assert_eq!(value.map.get(&SmallEnum::Foo), Some(&1));
        assert_eq!(value.map.get(&SmallEnum::Bar), Some(&2));
    }

    #[test]
    fn deserialize_str_keyed_map_reports_invalid_enum_values() {
        let err = toml::from_str::<SmallMap>("map = { Baz = 3 }")
            .expect_err("unknown enum key should fail");
        assert!(err.to_string().contains("unknown key"));
    }
}
