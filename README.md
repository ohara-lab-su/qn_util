# Qn-Util — Qn Network Analysis and RMC Constraint Generation

[![en](https://img.shields.io/badge/lang-en-red.svg)](https://github.com/ohara-lab-su/qn_util/blob/main/README.md)
[![ja](https://img.shields.io/badge/lang-ja-yellow.svg)](https://github.com/ohara-lab-su/qn_util/blob/main/README.ja.md)
[![ja](https://img.shields.io/badge/lang-ja-yellow.svg)](https://ohara-lab-su.github.io/qn_util/ja/index.html)

[![GitHub](https://img.shields.io/badge/GitHub-ohara--lab--su%2Fqn_util-181717?logo=github)](https://github.com/ohara-lab-su/qn_util/)
[![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-Documentation-222222?logo=githubpages&logoColor=white)](https://ohara-lab-su.github.io/qn_util/en/)
[![Ohara Lab](https://img.shields.io/badge/Ohara%20Lab-Home-4B5563?logo=homeassistant&logoColor=white)](https://ohara-lab-su.github.io/index.en.html)

Qn-Util is a Python library for constructing and analyzing Qn networks from VASP POSCAR/CONTCAR structures and JSON configuration files, and generating constraints for RMC_POT. It outputs Qn distributions, network bond information, and analysis reports.

## Main components

- `VaspStructure`: reads POSCAR files and computes distances and angles under periodic boundary conditions.
- `QnNetworkAnalysis`: constructs networks from JSON configuration, counts Qn species, and exports network and analysis data.
- `RmcConstraintWriter`: writes RMC_POT FNC=3 constraints and unit/network topology files.
- `QnNetworkBuilder`: compatibility facade that exposes analysis and output operations through a single entry point.
- `ConfigError`: exception used for invalid configuration and input data.

## Inputs and outputs

Inputs are a POSCAR/CONTCAR structure and JSON configuration. Outputs include `*_qn.csv`, `*_network_bonds.csv`, `*_report.txt`, `*_rmc_fnc3.fnc`, and `*_rmc_unit.top`. Consult the tutorials for configuration details and constraint definitions.

## Python API example

```python
from qn_util.qn_network_builder import QnNetworkBuilder

builder = QnNetworkBuilder.from_poscar("POSCAR", "qn_config.json")
print(builder.get_qn_counts())
builder.write_analysis("sample")
builder.write_all("sample")
```

`write_all` generates analysis results and RMC input files. See `example/` for example scripts and configuration.

## Installation

```bash
pip install qn-util
```

## Documentation

The tutorials cover Qn networks, FNC constraints, topology, JSON configuration, and possible extensions. The API Reference is generated from source docstrings.


## Documentation

[English](https://ohara-lab-su.github.io/qn_util/en/) / [日本語](https://ohara-lab-su.github.io/qn_util/ja/)
