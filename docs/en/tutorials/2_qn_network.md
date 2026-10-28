# Qn Network Analysis

## Purpose

To create FNC or Topology input for RMC_POT, it is necessary to determine which atoms form each local structure. This library represents a local structure using three layers:

```text
Central
      │
   Bridge
      │
Second
```

For SiO₂, this corresponds to:

```text
Si
 │
 O
 │
Si
```

The library analyzes this network to obtain Qn values and bond information.

---

# Network Model

The general network representation is:

```text
Central ----- Bridge ----- Second
```

Atomic species are specified in JSON. For SiO₂:

```text
Central : Si
Bridge  : O
Second  : Si
```

This design makes the analysis algorithm itself independent of SiO₂.

---

# Definition of Qn

Qn is defined as the number of Bridge atoms bonded to a Central atom that are also bonded to another Central atom. For example:

```text
      Si
      │
O─────┼─────O
│           │
Si         Si
```

This corresponds to `Q2`. A Bridge atom that is not connected to another Central atom is treated as a non-bridging oxygen (NBO) in the SiO₂ case.

---

# Analysis Procedure

The analysis proceeds approximately as follows:

```text
POSCAR
  ↓
Extract atomic species
  ↓
Select Central atoms
  ↓
Select Bridge atoms
  ↓
Select Second atoms
  ↓
Search Central–Bridge neighbors
  ↓
Search Bridge–Second neighbors
  ↓
Identify bridging atoms
  ↓
Determine Qn
  ↓
Generate bond information
  ↓
Generate angle information
  ↓
Generate FNC pairs
```

---

# 1. Extracting Atomic Species

Atomic species are read from POSCAR and classified according to the JSON entries `central_types`, `bridge_types`, and `second_types`. No distance-based selection is performed at this stage.

---

# 2. Central–Bridge Search

Neighbors between Central and Bridge atoms are identified using one of two methods: `cutoff` or `nearest`.

### cutoff

All atoms whose distances lie between `first_min` and `first_max` are selected.

### nearest

Only the specified number `N` of candidates satisfying the distance conditions are selected. For SiO₂, `nearest = 4` is commonly used to obtain SiO₄ tetrahedra consistently.

---

# 3. Bridge–Second Search

Second atoms are identified from the perspective of each Bridge atom. Two methods are available.

## second_use_first_network = false

A new distance-based search is performed around each Bridge atom to find Second atoms, for example O → Si. This is the simplest method.

## second_use_first_network = true

No separate Bridge–Second search is performed. Instead, the previously determined Central → Bridge network is reversed to construct Bridge → Second connections.

This approach preserves consistency between the Central–Bridge and Bridge–Second relationships and is recommended for SiO₂.

---

# 4. Identifying Bridging Atoms

For each Bridge atom, the number of connected Second atoms is examined. The original Central atom is excluded; only connections to **another Central atom** qualify the atom as a bridge.

---

# 5. Determining Qn

For each Central atom, the number of Bridge atoms connected to another Central atom is counted:

```text
4 shared bridges → Q4
3 shared bridges → Q3
2 shared bridges → Q2
```

Qn does not count all Bridge neighbors. It counts only those shared with another Central atom.

---

# Analysis Results

The analysis generates:

- Qn values for each Central atom.
- Bridge information.
- Bond information.
- Angle information.
- FNC pairs.

These are abstract network data independent of RMC formatting. The subsequent writer uses only these results to generate RMC_POT FNC and Topology files:

```text
Analysis
   ↓
Abstract network
   ↓
RMC files
```

This three-stage design keeps the analysis algorithm and output format independent.
