# README0
# RMC Topology / FNC Generator

## Overview

This library analyzes local network structures from VASP/POSCAR structures and automatically generates **FNC (`*.fnc`)** and **Topology (`*.top`)** files, which provide structural constraints for RMC_POT.

Originally implemented as a single script named `qn_network.py`, it was converted into a library and divided into classes according to their responsibilities for integration into the RMC-DFT workflow.

The library is intended to:

- Extract local networks from VASP structures.
- Analyze Qn structures.
- Convert the results into constraint formats understood by RMC_POT.
- Support automatic generation within an RMC-DFT workflow.

**This is not a general-purpose structural analysis library.** It is specifically designed to generate RMC_POT input files and follows the concepts underlying its FNC and Topology formats.

---

# Design Principles

Rather than implementing all processing in one large class, the library separates responsibilities among classes:

```text
             POSCAR
                │
                ▼
      VaspStructure
                │
                ▼
     QnNetworkAnalysis
                │
                ▼
   RmcConstraintWriter
                │
                ▼
        *.fnc / *.top
```

Users only need to interact with `QnNetworkBuilder`, the entry point that coordinates the workflow. They do not need to manage the details of the analysis.

---

# Class Structure

## VaspStructure

This class manages structural information only. Its responsibilities are:

- Reading POSCAR files.
- Storing lattice information.
- Storing atomic information.
- Calculating distances under periodic boundary conditions (PBC).
- Calculating PBC displacement vectors.
- Calculating angles.

It has no knowledge of Qn, Topology, or FNC. It is implemented as a pure structural-information class.

---

## QnNetworkAnalysis

This is the core class responsible for network analysis. It performs:

- Atomic-species extraction.
- Neighbor searches.
- Qn classification.
- Bond-information generation.
- Angle-information generation.
- Pair generation for FNC.

Results are retained as Qn distributions, bond information, bridge information, and FNC pairs. These results are independent of the RMC file format.

In other words, this class handles:

```text
Structure
    ↓
Analysis
    ↓
Abstract network information
```

---

## RmcConstraintWriter

This class converts analysis results into RMC_POT input formats. Its responsibilities are FNC output and Topology output. It performs no structural analysis; its only input is `QnNetworkAnalysis`.

---

## QnNetworkBuilder

This is the public entry point to the library. The builder executes the following sequence:

```text
POSCAR
      ↓
Structure
      ↓
Analysis
      ↓
Writer
```

Users only need to call:

```python
builder = QnNetworkBuilder.from_poscar(...)
builder.write_all(...)
```

The builder itself contains no analysis logic.

---

# Library Features

The library clearly separates **structural analysis** from **RMC input generation**. Consequently, changing the analysis algorithm does not require changing the writer. Conversely, changes to the RMC_POT file format do not affect the analysis algorithm.

This separation of responsibilities also makes it easier to support constraint formats other than RMC_POT in the future.

---

# FNC and Topology in RMC_POT

RMC_POT provides two principal approaches to preserving local structure:

- FNC (Fixed Neighbour Constraint).
- Topology (Bonded Potential).

This library can generate input files for both. Although their objectives are similar, their constraint mechanisms differ substantially.

The following chapters explain Qn network analysis, FNC, and Topology in that order.
