# JSON Configuration File

## Overview

Most aspects of the analysis algorithm can be configured through a JSON file. This allows users to change the target elements, neighbor-selection rules, bond parameters, angle parameters, and FNC types without modifying the program.

The current configuration targets SiO₂, but the design allows adaptation to other materials by changing JSON settings.

---

# Configuration Flow

```text
POSCAR
  ↓
JSON
  ↓
QnNetworkBuilder
  ↓
QnNetworkAnalysis
  ↓
RmcConstraintWriter
```

The builder passes JSON configuration to the relevant classes. The details of the analysis are determined by the configuration.

---

# network

The `network` configuration defines the Central → Bridge → Second network through `central_types`, `bridge_types`, and `second_types`.

For SiO₂:

```text
Central = Si
Bridge = O
Second = Si
```

These three atomic categories are the targets of network analysis.

---

# Neighbor Search

Neighbor searches are configured separately for `first` and `second`.

### First

The `first` settings control the Central → Bridge search, including the search method, minimum distance, maximum distance, and number of nearest neighbors.

### Second

The `second` settings control the Bridge → Second search, including minimum and maximum distances.

---

# second_use_first_network

This is an important configuration option.

- `false`: perform a new distance-based Bridge → Second search.
- `true`: reverse the previously constructed Central → Bridge network to identify Second atoms.

For SiO₂, `true` provides better consistency between the two sets of connections.

---

# Bond Parameters

Bond settings for Topology generally specify atom types, equilibrium distances, and force constants. The writer exports these settings to `[ bonds ]`.

---

# Angle Parameters

Angle settings generally specify angle types, equilibrium angles, and force constants. The writer exports these settings to `[ angles ]`.

---

# Atomic Information

The `[ atoms ]` section requires properties such as atom type, mass, and charge. These are also specified in JSON. The current implementation uses settings for Si and O.

---

# FNC Types

FNC types map identifiers to distance ranges. For example:

```text
Type 1 → 1.55–1.75 Å
```

The writer converts these settings into a `type_id → (rmin, rmax)` mapping and generates the FNC header.

---

# Role of JSON

Both the analysis and the writer use JSON, but they consume different settings.

The analysis uses neighbor-search rules, Qn classification settings, and atomic species. The writer uses bond and angle parameters, masses, charges, and FNC distance ranges.

---

# SiO₂-Specific Settings

The current JSON configuration assumes:

```text
Central = Si
Bridge = O
Second = Si
nearest = 4
```

This constructs SiO₄ tetrahedra. O–Si–O angle settings likewise assume tetrahedral geometry. Thus, although the configuration can be changed, its current values are specific to SiO₂.

---

# Applying the Library to Other Materials

The analysis algorithm uses the general Central → Bridge → Second model. It can therefore be adapted relatively easily to bridging-oxygen networks such as:

```text
P-O-P
B-O-B
Ge-O-Ge
Al-O-Si
```

Typically, the settings that need adjustment are atomic species, neighbor distances, bond parameters, angle parameters, and FNC ranges.

---

# Materials That Are Difficult to Handle

For metallic glasses containing Fe, Ni, or Cu, and for some liquids or simple crystals, the Central → Bridge → Second model may not be meaningful.

Changing JSON alone is insufficient for such systems; the analysis algorithm in `QnNetworkAnalysis` must be extended.

---

# Design Philosophy

The design is based on **algorithm + JSON**, not an algorithm wholly defined by JSON. JSON controls elements, distances, and parameters, while the program implements the Central–Bridge–Second network model.

This gives the library substantial reusability for SiO₂ and related network materials, but materials with fundamentally different network concepts require algorithmic extensions.

The next chapter examines how far the current SiO₂ implementation can be extended to other materials and discusses possible future development.
