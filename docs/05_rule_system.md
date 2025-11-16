# Reference: The DREIDING Rule System

After perception, every atom in an `AnnotatedMolecule` carries rich metadata: element, degree, lone pairs, hybridization, aromaticity, resonance state, and the smallest ring it participates in. The typing engine does not hard-code chemistry; instead, it evaluates TOML rules that describe how those annotations map to DREIDING atom types. This document is the complete guide to that rule layer.

## Rule Structure

Rules are declared as `[[rule]]` tables inside a TOML file. Each rule must provide four keys:

- `name` (`string`): descriptive identifier surfaced in diagnostics.
- `priority` (`integer`): conflict resolver; **larger values win** if multiple rules match an atom during the same iteration.
- `type` (`string`): the DREIDING atom type emitted when the rule fires.
- `conditions` (`table`): property checks an atom must satisfy. All listed checks must pass.

Example:

```toml
[[rule]]
name = "N_Trigonal_SP2"
priority = 200
type = "N_2"
conditions = { element = "N", steric_number = 3, is_aromatic = false }
```

At runtime, `typing::rules::parse_rules` converts the TOML into strongly typed `Rule` structures. `typing::rules::get_default_rules` lazily parses the embedded `resources/default.rules.toml`, so applications can either use the canonical ruleset directly or append their own entries before starting the typing engine.

## Available Conditions

Conditions operate on the immutable snapshot of an `AnnotatedAtom`. Because perception already computed lone pairs, hybridization, resonance, and ring membership, rules can simply read fields and avoid bespoke chemistry code. Every key is optional; omitting a key turns it into a wildcard.

The following table details every valid key that can be used inside the `conditions` table.

| Key                           | Type    | Description                                                                                                                                                      |
| ----------------------------- | ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Atom-Intrinsic Properties** |         | Properties derived from the atom itself.                                                                                                                         |
| `element`                     | String  | The atom's element symbol (e.g., `"C"`, `"Na"`). Must be a valid symbol.                                                                                         |
| `formal_charge`               | Integer | The formal charge of the atom (e.g., `1`, `0`, `-1`).                                                                                                            |
| `degree`                      | Integer | The number of directly bonded neighbor atoms.                                                                                                                    |
| `lone_pairs`                  | Integer | The number of lone electron pairs, as calculated during the Perception Phase.                                                                                    |
| `steric_number`               | Integer | The sum of `degree` and `lone_pairs`. A primary indicator of geometry.                                                                                           |
| `hybridization`               | String  | The perceived hybridization state. Valid values: `"SP"`, `"SP2"`, `"SP3"`, `"Resonant"`, `"None"`.                                                               |
| `is_in_ring`                  | Boolean | `true` if the atom is part of any detected ring system.                                                                                                          |
| `is_aromatic`                 | Boolean | `true` if the atom is part of a perceived aromatic system.                                                                                                       |
| `is_anti_aromatic`            | Boolean | `true` if perception tagged the atom as belonging to an anti-aromatic ring.                                                                                      |
| `is_resonant`                 | Boolean | `true` if resonance analysis marked the atom as delocalized (e.g., phenoxide oxygen).                                                                            |
| `smallest_ring_size`          | Integer | The size of the smallest ring the atom belongs to (e.g., `5` for furan).                                                                                         |
| **Neighbor-Based Properties** |         | Properties derived from the atom's immediate neighbors.                                                                                                          |
| `neighbor_elements`           | Table   | Specifies the **exact counts** of neighboring elements. Atoms not listed are assumed to be zero.                                                                 |
| `neighbor_types`              | Table   | Specifies the **exact counts** of the **final assigned types** of neighboring atoms. This is the key condition that enables context-dependent, iterative typing. |

**Example of `neighbor_elements`:**
The following condition matches a hydrogen atom bonded to exactly two boron atoms (as in diborane).

```toml
conditions = { element = "H", degree = 2, neighbor_elements = { B = 2 } }
```

**Example of `neighbor_types`:**
This condition would match a carbon atom bonded to exactly one `C_3` atom and three `H_` atoms.

```toml
conditions = { element = "C", neighbor_types = { "C_3" = 1, "H_" = 3 } }
```

### The Role of `priority` and `neighbor_types`

- **Priority:** The `priority` key is the sole mechanism for resolving conflicts. When an atom matches multiple rules, the one with the highest `priority` value is definitively chosen in that iteration.
- **Iteration trigger:** `neighbor_types` refers to already-assigned neighbor atom types. Early rounds may skip these rules while neighbors are still untyped. The engine keeps iterating, seeding newly determined types back into the graph, until every atom is stable. See [Typing Engine](./03_typing_engine.md) for the convergence strategy.

## Default Ruleset Philosophy and Key Atom Types

`resources/default.rules.toml` tracks the original DREIDING priorities while embracing the richer perception data. The layout is intentionally layered:

1. **500+** – Exotic safeties and overrides (e.g., diborane bridging hydrogens).
2. **400s** – Delocalized or aromatic atoms (`*_R`, resonance-stabilized heteroatoms) that must outrank geometry-only rules.
3. **100–300** – VSEPR-driven workhorses keyed off steric number and hybridization.
4. **<100** – Simple fallbacks such as halogens, alkali/alkaline-earth metals, and default hydrogens.

Representative entries are summarized below (table retained for quick reference):

| Atom Type         | DREIDING Description          | Key Rule Condition(s) in `dreiding.rules.toml`                    | Priority |
| :---------------- | :---------------------------- | :---------------------------------------------------------------- | :------: |
| `H_`              | Standard Hydrogen             | `{ element = "H" }`                                               |    1     |
| `H_HB`            | Hydrogen-Bonding Hydrogen     | `{ element = "H", neighbor_elements = { O = 1 } }` or `{ N = 1 }` |   ~250   |
| `H_b`             | Bridging Hydrogen (Diborane)  | `{ degree = 2, neighbor_elements = { B = 2 } }`                   |   500    |
| `C_3`             | sp³ Tetrahedral Carbon        | `{ element = "C", steric_number = 4 }`                            |   100    |
| `C_2`             | sp² Trigonal Carbon           | `{ element = "C", steric_number = 3, is_aromatic = false }`       |   200    |
| `C_1`             | sp Linear Carbon              | `{ element = "C", steric_number = 2 }`                            |   300    |
| `C_R`             | Resonant/Aromatic Carbon      | `{ element = "C", is_aromatic = true }`                           |   400    |
| `N_3`             | sp³ Nitrogen (Amine/Ammonium) | `{ element = "N", steric_number = 4 }`                            |   100    |
| `N_2`             | sp² Nitrogen (Imine/Amide)    | `{ element = "N", steric_number = 3, is_aromatic = false }`       |   200    |
| `N_R`             | Resonant/Aromatic Nitrogen    | `{ element = "N", is_aromatic = true }`                           |   400    |
| `O_3`             | sp³ Oxygen (Ether/Alcohol)    | `{ element = "O", steric_number = 4 }`                            |   100    |
| `O_2`             | sp² Oxygen (Carbonyl)         | `{ element = "O", steric_number = 3 }`                            |   200    |
| `O_R`             | Resonant Oxygen (Phenol)      | `{ element = "O", hybridization = "Resonant" }`                   |   401    |
| `S_R`             | Resonant Sulfur (Thiophene)   | `{ element = "S", hybridization = "Resonant" }`                   |   400    |
| `P_3`             | sp³ Phosphorus (Phosphate)    | `{ element = "P", steric_number = 4 }`                            |   100    |
| `S_3`             | sp³ Sulfur (Thiol/Sulfide)    | `{ element = "S", hybridization = "SP3" }`                        |   100    |
| `F_`, `Cl_`, etc. | Halogens                      | `{ element = "F" }`, etc.                                         |    50    |
| `Na`, `Ca`, etc.  | Metal Ions                    | `{ element = "Na" }`, etc.                                        |    20    |

## How to Extend the Rule System

Customizing typing means editing TOML, not Rust. Typical workflow:

1. **Author a TOML snippet** (e.g., `my_copper_rules.toml`).
2. **Pick priorities carefully.** Choose values that let your rules coexist with (or outrank) the defaults.
3. **Load rules at runtime.** Parse the TOML with `dreid_typer::rules::parse_rules` and pass the resulting slice into `assign_topology_with_rules`. (If you want to extend the canonical DREIDING file, copy `resources/default.rules.toml` into your project, edit it, and parse that content before appending your custom entries.)

```toml
# my_copper_rules.toml
[[rule]]
name = "Ion_Cu_Divalent"
priority = 20
type = "Cu+2"
conditions = { element = "Cu", formal_charge = 2 }
```

```rust
use dreid_typer::{assign_topology_with_rules, rules, MolecularGraph, MolecularTopology, TyperError};

fn type_with_custom_rules(graph: &MolecularGraph) -> Result<MolecularTopology, TyperError> {
   // Load the default ruleset.
   let mut all_rules = rules::parse_rules(include_str!("dreiding.rules.toml"))?;

   // Append or override entries programmatically.
   all_rules.extend(rules::parse_rules(include_str!("my_copper_rules.toml"))?);

   assign_topology_with_rules(graph, &all_rules)
}
```

Because the engine merely consumes structured data, you can version-control TOML files, generate them from other toolchains, or even ship different rulesets for different force fields—all without recompiling `dreid-typer`.
