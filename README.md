# Qn-Util — Qn Network Analysis and RMC Constraint Generation

[![en](https://img.shields.io/badge/lang-en-red.svg)](https://github.com/ohara-lab-su/qn_util/blob/main/README.md)
[![ja](https://img.shields.io/badge/lang-ja-yellow.svg)](https://github.com/ohara-lab-su/qn_util/blob/main/README.ja.md)

[![GitHub](https://img.shields.io/badge/GitHub-ohara--lab--su%2Fqn_util-181717?logo=github)](https://github.com/ohara-lab-su/qn_util/)
[![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-Documentation-222222?logo=githubpages&logoColor=white)](https://ohara-lab-su.github.io/qn_util/ja/)
[![Ohara Lab](https://img.shields.io/badge/Ohara%20Lab-Home-4B5563?logo=homeassistant&logoColor=white)](https://ohara-lab-su.github.io/index.ja.html)

Qn-Util is a Python library for constructing and analyzing Qn networks and generating constraint files for RMC_POT, using VASP POSCAR/CONTCAR structures and JSON configuration files as input. It outputs Qn distributions, network connectivity information, and analysis reports.

The library is primarily intended for generating RMC_POT constraint files and is used for research and study purposes.

## Main Components

- `VaspStructure`: Reads POSCAR files and calculates distances and angles under periodic boundary conditions.
- `QnNetworkAnalysis`: Builds networks from JSON configurations, counts Qn species, and outputs connectivity information and analysis results.
- `RmcConstraintWriter`: Generates FNC=3 constraint files, unit definitions, and network topology files for RMC_POT.
- `QnNetworkBuilder`: A facade that preserves existing calling conventions and provides combined access to analysis and output functions.
- `ConfigError`: An exception for configuration and input validation errors.

## Input and Output

Inputs are structural files in POSCAR/CONTCAR format and JSON configuration files. Outputs include `*_qn.csv`, `*_network_bonds.csv`, `*_report.txt`, `*_rmc_fnc3.fnc`, and `*_rmc_unit.top`. See the tutorials for required configuration parameters and the meaning of the constraints.

## Python API Example

```python
from qn_util.qn_network_builder import QnNetworkBuilder

builder = QnNetworkBuilder.from_poscar("POSCAR", "qn_config.json")
print(builder.get_qn_counts())
builder.write_analysis("sample")
builder.write_all("sample")
```

`write_all` generates analysis results and RMC input files. See `scripts/` for concrete configuration examples.

## Installation

```bash
pip install qn-util
```

## Installing Example Scripts

After installing the package, run the following command in your working directory to copy the files bundled under `scripts/` (Python scripts and `qn_config.json`) into the current directory.

```bash
qn_util install
```

If files with the same names already exist, the command exits with an error without overwriting them. Specify `--force` only when you want to overwrite existing files.

```bash
qn_util install --force
```

## Documentation

The tutorials cover Qn networks, FNC, topology, JSON configuration, and future development. The API Reference is generated from source docstrings.

[日本語](https://ohara-lab-su.github.io/qn_util/ja/) / [English](https://ohara-lab-su.github.io/qn_util/en/)
