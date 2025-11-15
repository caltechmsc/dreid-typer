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
