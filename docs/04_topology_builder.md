# Phase 3: The Topology Builder

With atom types resolved, the builder translates the annotated molecule into a `MolecularTopology`. This stage is pure graph traversal — no additional chemistry is inferred — but it's where canonical force-field terms emerge.

`builder::build_topology` takes two inputs:

1. The immutable `AnnotatedMolecule` output of perception.
2. The `Vec<String>` of atom types returned by the typing engine.

It produces `MolecularTopology { atoms, bonds, angles, torsions, inversions }`, all deduplicated and ready for downstream MD engines.

```mermaid
graph TD
    A(<b>AnnotatedMolecule</b>) --> B{builder::build_topology}
    T(<b>Atom Types</b>) --> B
    B --> C(<b>MolecularTopology</b>)
```

## Atom Table

`build_atoms` walks the annotated atoms and copies their element, hybridization, and ID while splicing in the final type string (`atom_types[ann_atom.id]`). This produces the topology's `atoms` vector.

## Connectivity Terms

Every interaction term uses the molecule's adjacency lists and bond table, which already reflect Kekulé-expanded bond orders.

### Bonds

`build_bonds` maps each `BondEdge` to a `Bond`, sorting the atom indices so that `(i, j)` and `(j, i)` collapse to the same representation. The resulting set is stored in a `HashSet<Bond>` to prevent duplicates before being collected into a `Vec`.

### Angles (`build_angles`)

For each atom `j` (the angle center), consider all unordered pairs of neighbors `(i, k)` taken from `adjacency[j]`. Each pair yields `Angle::new(i, j, k)`, which internally sorts the outer atoms to maintain canonical order. Because combinations are generated without repetition, every unique `i-j-k` angle appears exactly once.

Pseudocode:

```text
for center in atoms:
    neighbors = adjacency[center]
    for each unordered pair (i, k) in neighbors:
        angles.insert(Angle::new(i, center, k))
```

### Torsions (`build_torsions`)

Torsions are enumerated around each bond `j-k`:

1. Iterate over every stored bond.
2. For each neighbor `i` of `j` (excluding `k`) and each neighbor `l` of `k` (excluding `j` and `i`), emit `Torsion::new(i, j, k, l)`.
3. The constructor compares `(i, j, k, l)` to its reverse `(l, k, j, i)` and keeps the lexicographically smaller tuple to guarantee uniqueness.

This approach naturally covers both directions (i.e., `i-j-k-l` and `l-k-j-i`) without generating duplicates.

### Inversions (`build_inversions`)

Inversions enforce planarity at trigonal centers. The builder scans every atom and checks two conditions:

1. Degree equals 3.
2. Hybridization equals `Hybridization::SP2` or `Hybridization::Resonant`.

Per the DREIDING paper, **each planar center generates three inversion terms**, with each neighbor taking turn as the "axis":

For center I with neighbors {J, K, L}:

- Inversion(center=I, axis=J, plane={K, L})
- Inversion(center=I, axis=K, plane={J, L})
- Inversion(center=I, axis=L, plane={J, K})

The constructor `Inversion::new(center, axis, plane1, plane2)` sorts only the two plane atoms (not the axis), ensuring the three terms per center remain distinct.

## Why Canonical Forms Matter

- **Deduplication:** All intermediate collections are `HashSet`s, so deterministic ordering of atom IDs is required to detect duplicates.
- **Stable output:** Simulation pipelines downstream can diff or cache topologies knowing that rerunning the builder yields identical ordering.
- **Serialization friendliness:** Canonical tuples simplify hashing/serialization and reduce noise in reference files.

## Error Handling

The builder assumes its inputs are valid. By the time `build_topology` runs, perception has already verified graph connectivity, and typing has succeeded. Therefore, the functions are infallible and cannot produce errors on their own.

## Summary

The topology builder does not invent new chemistry; it formalizes the perceived molecule into the geometric primitives expected by DREIDING-compatible engines. Canonical ordering and set-based generation ensure reproducible, duplicate-free interaction lists every time.
