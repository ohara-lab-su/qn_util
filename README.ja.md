# Qn-Util — Qnネットワーク解析とRMC制約生成

[![en](https://img.shields.io/badge/lang-en-red.svg)](https://github.com/ohara-lab-su/qn_util/blob/main/README.md)
[![ja](https://img.shields.io/badge/lang-ja-yellow.svg)](https://github.com/ohara-lab-su/qn_util/blob/main/README.ja.md)

[![GitHub](https://img.shields.io/badge/GitHub-ohara--lab--su%2Fqn_util-181717?logo=github)](https://github.com/ohara-lab-su/qn_util/)
[![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-Documentation-222222?logo=githubpages&logoColor=white)](https://ohara-lab-su.github.io/qn_util/ja/)
[![Ohara Lab](https://img.shields.io/badge/Ohara%20Lab-Home-4B5563?logo=homeassistant&logoColor=white)](https://ohara-lab-su.github.io/index.ja.html)

Qn-Util は VASP の POSCAR/CONTCAR 構造と JSON 設定を入力として、Qn ネットワークの構築・解析および RMC_POT 用の制約ファイル生成を行う Python ライブラリです。Qn 分布、ネットワーク結合、解析レポートを出力します。

## 主な構成

- `VaspStructure`：POSCAR の読込、周期境界条件下の距離・角度計算。
- `QnNetworkAnalysis`：JSON 設定に基づくネットワーク構築、Qn 集計、結合情報と解析結果の出力。
- `RmcConstraintWriter`：RMC_POT 用 FNC=3 制約ファイル、ユニットおよびネットワーク topology の出力。
- `QnNetworkBuilder`：既存の呼び出しを維持する facade。解析と出力機能をまとめて利用できます。
- `ConfigError`：設定・入力の検証に関する例外。

## 入力と出力

入力は POSCAR/CONTCAR 形式の構造ファイルと JSON 設定です。出力には `*_qn.csv`、`*_network_bonds.csv`、`*_report.txt`、`*_rmc_fnc3.fnc`、`*_rmc_unit.top` などがあります。必要な設定値と制約の意味はチュートリアルを参照してください。

## Python API の例

```python
from qn_util.qn_network_builder import QnNetworkBuilder

builder = QnNetworkBuilder.from_poscar("POSCAR", "qn_config.json")
print(builder.get_qn_counts())
builder.write_analysis("sample")
builder.write_all("sample")
```

`write_all` は解析結果と RMC 入力ファイルを生成します。具体的な設定例は `example/` を参照してください。

## インストール

```bash
pip install qn-util
```

## ドキュメント

チュートリアルでは Qn ネットワーク、FNC、topology、JSON 設定、今後の展開を説明します。API Reference はソースの docstring から生成します。


## ドキュメント

[日本語](https://ohara-lab-su.github.io/qn_util/ja/) / [English](https://ohara-lab-su.github.io/qn_util/en/)
