# Reference: The DREIDING Rule System

The `dreid-typer` engine's power and flexibility stem from its rule-based design. The entire atom typing logic is defined in human-readable TOML files, allowing for easy inspection, modification, and extension without altering the core library code. This document serves as the definitive technical reference for the rule system.

## Rule Structure

All rules are defined within a TOML file as an array of tables, with each rule prefixed by `[[rule]]`. Each rule is an independent object containing four mandatory keys.

- `name` (`string`): A unique, descriptive name for the rule (e.g., `"C_Aromatic"`). This is used for debugging and identification.
- `priority` (`integer`): A number that determines the rule's precedence in cases of conflict. **Higher numbers correspond to higher priority** and are evaluated first.
- `type` (`string`): The DREIDING atom type string to be assigned if all conditions are met (e.g., `"C_R"`).
- `conditions` (`table`): A table of one or more conditions that an atom must satisfy for the rule to be considered a match.

**Example of a complete rule:**

```toml
[[rule]]
# A human-readable identifier for the rule.
name = "N_Trigonal_SP2"

# This rule will be chosen over any rule with a priority < 200.
priority = 200

# The atom type to be assigned on a successful match.
type = "N_2"

# A collection of properties the atom must have.
conditions = { element = "N", steric_number = 3, is_aromatic = false }
```

## Available Conditions

The `conditions` table is the heart of a rule. It defines the specific chemical context that an atom must match. An atom must satisfy **all** specified conditions within a single rule for that rule to apply.

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
- **Iteration Trigger:** The `neighbor_types` condition is special. Because it depends on the _final_ types of other atoms, rules containing it may not match in the first round of typing. The iterative engine will continue to run, propagating new type information, until these context-dependent rules can be satisfied and the system reaches a stable state. For a detailed explanation, see the [Typing Engine documentation](./03_typing_engine.md).

## Default Ruleset Philosophy and Key Atom Types

The default ruleset (`dreiding.rules.toml`) is designed to be a faithful and robust implementation of the original DREIDING philosophy. Its structure follows a clear hierarchy:

1. **Highest Priority (500+):** Extremely specific and rare cases (e.g., bridging hydrogens).
2. **High Priority (400s):** Strong, non-local chemical features like aromaticity that should override simpler geometric definitions.
3. **Medium Priority (100-300):** The main workhorse rules based on geometry (steric number and hybridization) from VSEPR theory.
4. **Low Priority (1-99):** General fallback rules for common cases like halogens and a default for standard hydrogen.

The table below summarizes some of the most common atom types and the key conditions that define them in the default ruleset.

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
| `P_3`             | sp³ Phosphorus (Phosphate)    | `{ element = "P", steric_number = 4 }`                            |   100    |
| `S_3`             | sp³ Sulfur (Thiol/Sulfide)    | `{ element = "S", hybridization = "SP3" }`                        |   100    |
| `F_`, `Cl_`, etc. | Halogens                      | `{ element = "F" }`, etc.                                         |    50    |
| `Na`, `Ca`, etc.  | Metal Ions                    | `{ element = "Na" }`, etc.                                        |    20    |

## How to Extend the Rule System

The system is designed for easy extension. To add support for new elements or define custom types:

1. **Create a custom TOML file** (e.g., `my_copper_rules.toml`).
2. **Define your new rules.** Ensure you choose a priority that correctly interacts with existing rules. For a new metal ion, a low priority is appropriate.

   ```toml
   # my_copper_rules.toml
   [[rule]]
   name = "Ion_Cu_Divalent"
   priority = 20
   type = "Cu+2"
   conditions = { element = "Cu", formal_charge = 2 }
   ```

3. **Load and use the rules in your code.** You can either combine them with the default rules or use them exclusively.

   ```rust
   use dreid_typer::{rules, assign_topology_with_rules, MolecularGraph};

   fn type_molecule_with_custom_rules(graph: &MolecularGraph) {
       // Load default rules
       let mut all_rules = rules::get_default_rules().unwrap().to_vec();

       // Load and append custom rules
       let custom_rules_str = include_str!("my_copper_rules.toml");
       let custom_rules = rules::parse_rules(custom_rules_str).unwrap();
       all_rules.extend(custom_rules);

       // Run the typer with the combined ruleset
       let topology = assign_topology_with_rules(graph, &all_rules).unwrap();
       // ...
   }
   ```
