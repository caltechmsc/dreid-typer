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
