use crate::core::properties::{Element, Hybridization};
use serde::{Deserialize, Deserializer, de};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use std::sync::OnceLock;

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    pub name: String,
    pub priority: i32,
    #[serde(rename = "type")]
    pub result_type: String,
    pub conditions: Conditions,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Default)]
#[serde(deny_unknown_fields)]
pub struct Conditions {
    #[serde(default)]
    pub element: Option<Element>,
    #[serde(default)]
    pub formal_charge: Option<i8>,

    #[serde(default)]
    pub degree: Option<u8>,
    #[serde(default)]
    pub is_in_ring: Option<bool>,

    #[serde(default)]
    pub lone_pairs: Option<u8>,
    #[serde(default)]
    pub hybridization: Option<Hybridization>,

    #[serde(default)]
    pub is_aromatic: Option<bool>,
    #[serde(default)]
    pub is_anti_aromatic: Option<bool>,
    #[serde(default)]
    pub is_resonant: Option<bool>,

    #[serde(default, deserialize_with = "deserialize_str_keyed_map")]
    pub neighbor_elements: HashMap<Element, u8>,
    #[serde(default)]
    pub neighbor_types: HashMap<String, u8>,
}

#[derive(Deserialize)]
struct Ruleset {
    #[serde(rename = "rule")]
    rules: Vec<Rule>,
}

pub fn parse_rules(content: &str) -> Result<Vec<Rule>, toml::de::Error> {
    let ruleset: Ruleset = toml::from_str(content)?;
    Ok(ruleset.rules)
}

static DEFAULT_RULES: OnceLock<Vec<Rule>> = OnceLock::new();

const DEFAULT_RULES_TOML: &str = include_str!("../../resources/default.rules.toml");

pub fn get_default_rules() -> &'static [Rule] {
    DEFAULT_RULES.get_or_init(|| {
        parse_rules(DEFAULT_RULES_TOML)
            .expect("Failed to parse embedded default DREIDING rules. This is a library bug.")
    })
}

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
