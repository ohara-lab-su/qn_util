# SiO₂ Implementation and Extension to Other Materials

## Overview

The library is designed around a Central → Bridge → Second network model. Although the current implementation targets SiO₂, the underlying algorithm has been generalized to some extent.

However, the definition of Qn, the construction of Topology files, and bond parameters still strongly reflect SiO₂ structures. This chapter distinguishes general-purpose components from SiO₂-specific ones.

---

# SiO₂-Specific Components

The four main SiO₂-specific aspects are:

- The definition of Qn.
- SiO₄ tetrahedra.
- Bond parameters.
- Angle parameters.

These cannot necessarily be reused simply by changing element names.

---

# 1. Definition of Qn

Qn originally describes the number of bridging oxygen atoms associated with an SiO₄ tetrahedron in silicates.

For example, Q4 means that all four oxygens are bridging oxygens, while Q3 means that one oxygen is non-bridging.

The term Qn is therefore specific to silicate chemistry. Although the library uses this name internally, its underlying operation is more generally described as **counting bridging atoms**.

---

# 2. SiO₄ Tetrahedra

The default neighbor-search setting is `nearest = 4`, reflecting the assumption of SiO₄ tetrahedra. The same setting can be used for PO₄ tetrahedra.

For octahedral units such as TiO₆ or AlO₆, however, the setting must be changed to `nearest = 6`.

Thus, the neighbor-search algorithm is general, but its default values target SiO₂.

---

# 3. Bond Parameters

Topology bond parameters currently assume Si–O bonds, including their equilibrium distances. Other materials require parameters for bonds such as P–O, Ge–O, or Al–O. These changes can be made through JSON.

---

# 4. Angle Parameters

The current implementation focuses on O–Si–O and Si–O–Si angles. Materials containing PO₄, BO₃, or GeO₄ structural units require corresponding changes.

---

# General-Purpose Components

## Structure

`VaspStructure` is fully general-purpose. It depends neither on the elements nor on the network model; it handles POSCAR structural data.

## Neighbor Search

Neighbor searches are not intrinsically material-specific. They use distances and settings such as `nearest`.

## Network Model

The Central → Bridge → Second model is relatively general and can represent networks such as P–O–P, Ge–O–Ge, Al–O–Si, and B–O–B.

## Writer

The writer consumes analysis results without embedding knowledge of Si or O. It simply converts the results into RMC_POT input formats, making it the most general-purpose component.

---

# Applicable Materials

## Materials That Are Relatively Easy to Support

```text
SiO2
GeO2
P2O5
B2O3
AlPO4
Al-Si-O systems
```

These materials have the concept of bridging oxygen, so the bridge-network model can be reused relatively directly.

## Materials Requiring Minor Adjustments

```text
Li2SiO3
Li2Si2O5
Na2SiO3
```

These contain modifier cations such as Li or Na, which are not Bridge atoms. Qn classification should consider only Si → O → Si connections.

The neighbor-search algorithm can still be used, but the analysis conditions must be adjusted in JSON.

## Materials Requiring Algorithm Changes

```text
Li3PO4
Li3PS4
NASICON
```

Li is not a network-forming atom in these systems. If the analysis targets Li diffusion or Li sites, the bridge-network concept does not match the objective. Direct use of `QnNetworkAnalysis` is therefore difficult.

## Materials That Are Difficult to Support

Metallic glasses containing Fe, Ni, Cu, Pd, or Pt do not have bridging oxygen, and the Qn concept does not apply.

Such systems require substantially different approaches, for example neighbor graphs, Voronoi analysis, coordination-number (CN) analysis, or Delaunay analysis.

---

# Future Extensions

The current class is named `QnNetworkAnalysis`, although its underlying operation is better described as **Bridge Network Analysis**.

A future design could introduce an abstract `NetworkAnalysis` class with specialized implementations:

```text
NetworkAnalysis
  ├── QnNetworkAnalysis
  ├── PO4NetworkAnalysis
  ├── OctahedralNetworkAnalysis
  └── MetalNetworkAnalysis
```

The builder would select an analysis class based on JSON configuration. This would allow the analysis algorithm to be replaced without changing the writer.

---

# Position of This Library

Rather than merely being a library for Qn analysis, it is designed to **analyze local networks and convert the results into RMC_POT constraint input**.

Qn analysis is central to the current SiO₂ implementation, but the architecture follows three stages:

```text
Structural analysis
        ↓
Abstract network
        ↓
RMC_POT input
```

This separation of responsibilities makes it relatively straightforward to add new analysis algorithms, materials, and constraint formats.

This concludes the discussion of the library's design principles, analysis algorithms, RMC_POT FNC/Topology generation, and applicability to SiO₂ and related materials.
