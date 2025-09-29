use crate::core::{Element, Hybridization, error::TyperError};
use serde::{Deserialize, Deserializer, de};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

mod default;

fn from_str<'de, T, D>(deserializer: D) -> Result<T, D::Error>
where
    T: FromStr,
    T::Err: fmt::Display,
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    T::from_str(&s).map_err(de::Error::custom)
}

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

#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Default)]
#[serde(deny_unknown_fields)]
pub struct Conditions {
    #[serde(default, deserialize_with = "from_str_optional")]
    pub element: Option<Element>,
    #[serde(default)]
    pub degree: Option<u8>,
    #[serde(default, deserialize_with = "from_str_optional")]
    pub hybridization: Option<Hybridization>,
    #[serde(default)]
    pub is_in_ring: Option<bool>,
    #[serde(default)]
    pub is_aromatic: Option<bool>,
    #[serde(default)]
    pub smallest_ring_size: Option<u8>,

    #[serde(default, deserialize_with = "hash_map_from_str_keys")]
    pub neighbor_elements: HashMap<Element, u8>,
    #[serde(default)]
    pub neighbor_types: HashMap<String, u8>,
}

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

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Rule {
    pub name: String,
    pub priority: i32,
    #[serde(rename = "type")]
    pub result_type: String,
    pub conditions: Conditions,
}

#[derive(Debug, Clone, Deserialize)]
struct Ruleset {
    #[serde(rename = "rule")]
    rules: Vec<Rule>,
}

use once_cell::sync::Lazy;

static DEFAULT_RULES: Lazy<Result<Vec<Rule>, TyperError>> =
    Lazy::new(|| parse_rules(default::DEFAULT_RULES_TOML));

pub fn parse_rules(content: &str) -> Result<Vec<Rule>, TyperError> {
    let ruleset: Ruleset =
        toml::from_str(content).map_err(|e| TyperError::RuleParse(e.to_string()))?;
    Ok(ruleset.rules)
}

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
            degree = 4
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
                degree: Some(4),
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
        assert_eq!(rules.len(), 44);
    }
}
