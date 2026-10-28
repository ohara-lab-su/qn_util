# FNC (Fixed Neighbour Constraint)

## What Is FNC?

FNC (Fixed Neighbour Constraint) is an RMC_POT feature for constraining local structures.

Conventional Reverse Monte Carlo methods evaluate agreement with experimental data, potentially producing structures in which chemical bonds break, coordination numbers change substantially, or networks collapse.

FNC introduces constraints into the Monte Carlo acceptance decision to keep specified neighboring atoms within prescribed distance ranges.

---

# Role of FNC in RMC_POT

FNC is used with:

```text
switch = 1
switch = 2
switch = 3
```

In contrast, `switch = 4` uses bonded potentials based on a GROMACS-style Topology file. FNC and Topology therefore overlap in purpose but differ fundamentally in their constraint mechanisms.

---

# Concept of FNC

FNC constrains **distances only**. For example, the Si–O distance may be required to remain within:

```text
1.55 Å – 1.75 Å
```

Moves outside this range are rejected or have a reduced acceptance probability. FNC therefore specifies an **allowed distance range**, not a bond potential.

---

# FNC Generation in This Library

FNC files are generated from Central → Bridge bond information obtained by Qn network analysis. The analysis retains only:

```text
(atom_i,
 atom_j,
 type_id)
```

Here, `type_id` identifies a distance range. The writer converts these identifiers into RMC_POT FNC definitions.

---

# FNC Types

In RMC_POT, distance ranges are assigned type identifiers, for example:

```text
Type 1 → 1.55–1.75 Å
Type 2 → 2.45–2.70 Å
```

The library reads `fnc_ranges` from JSON and maps each `type_id` to `(rmin, rmax)`. Pairs with the same distance constraint can share one type.

---

# FNC Pairs

The analysis produces pairs in the form:

```text
(atom_i,
 atom_j,
 type_id)
```

The writer reorganizes them into the RMC_POT format, which specifies, for each atom, its constrained neighbors and the corresponding constraint types. Thus, output is organized as an **atom list**, rather than a simple pair list.

---

# FNC File Structure

The generated FNC file has approximately the following layout:

```text
Number of FNC types
Minimum distance for each type
Maximum distance for each type
Number of atoms
Atom 1
Number of constraints
Partner atom
Constraint type
Atom 2
...
```

For each atom, the number of constraints is followed by the constrained partners and their types. This corresponds to the fixed-neighbor constraint input expected by RMC_POT.

---

# Pair Symmetry

A constraint between Si and O is stored in both directions, Si → O and O → Si, because RMC_POT expects a list of constraints for each atom. Pairs are therefore symmetrized internally before output.

---

# Type Management

The writer checks that type identifiers form a consecutive sequence `1, 2, 3, ...`. RMC_POT reads the type-to-distance-range mapping in numerical order. Missing identifiers are treated as configuration errors.

---

# Limitations of FNC

FNC constrains only distances; it cannot directly represent bond angles, dihedrals, or stiffness.

For tetrahedral networks such as SiO₂, constraining Si–O distances alone does not preserve O–Si–O or Si–O–Si angles. When stricter local-structure preservation is needed, `switch = 4` and Topology-based bonded potentials are more appropriate.

---

# Position in the Library

`QnNetworkAnalysis` generates abstract network information, and `RmcConstraintWriter` converts that information into FNC input. The analysis algorithm does not depend on FNC syntax, and the writer does not perform Qn classification. This separation allows both components to be maintained independently.

The next chapter describes generation of Topology (`*.top`) files containing bonds, angles, and molecules.
