#!/usr/bin/env python
"""
K.NAKADA, kengo.nakada@gmail.com
"""
from typing import Any, Dict, List, Optional, Sequence, Tuple
import itertools
import json
import math
import re
from collections import Counter, defaultdict
import csv

import numpy as np

from qn_util.config_error import ConfigError
from qn_util.vasp_structure import VaspStructure


class QnNetworkAnalysis:
    """Qn network の構築、解析、解析結果出力を行うクラス。"""

    def __init__(
        self,
        structure: Optional["VaspStructure"] = None,
        config: Optional[Dict[str, Any]] = None,
    ) -> None:

        self.poscar = structure
        self.structure = structure
        self.config = config
        self.first_map = defaultdict(list)
        self.bridge_second_map = defaultdict(list)
        self.qn_by_central = {}
        self.bridges_by_central = {}
        self.second_by_bridge_for_central = {}
        self.bond_pairs = set()
        self.fnc_pairs = []
        self.fnc_type_names = {}
        self.fnc_error = None
        self.angles = []

        if structure is not None and config is not None:
            data = QnNetworkAnalysis.build_network(structure, config)
            try:
                QnNetworkAnalysis.build_fnc3_pairs(structure, config, data)
            except ConfigError as exc:
                self.fnc_error = str(exc)
                data.fnc_pairs = []
                data.fnc_type_names = {}
            self.first_map = data.first_map
            self.bridge_second_map = data.bridge_second_map
            self.qn_by_central = data.qn_by_central
            self.bridges_by_central = data.bridges_by_central
            self.second_by_bridge_for_central = data.second_by_bridge_for_central
            self.bond_pairs = data.bond_pairs
            self.fnc_pairs = data.fnc_pairs
            self.fnc_type_names = data.fnc_type_names
            self.angles = data.angles

    @classmethod
    def from_poscar(
        cls,
        poscar_name: str,
        config_name: str,
    ) -> "QnNetworkAnalysis":
        """VASP 構造ファイルと JSON 設定から QnNetworkAnalysis を作る。"""
        config = cls.load_config(config_name)
        atom_types = cls.get_structure_atom_types(config)
        structure = VaspStructure.read_poscar(poscar_name, atom_types)
        return cls(structure, config)

    def get_qn_counts(self) -> Dict[int, int]:
        """Qn 分布を返す。"""
        counter = Counter(self.qn_by_central.values())
        result = {}
        for key in sorted(counter.keys()):
            result[int(key)] = int(counter[key])
        return result

    def write_qn_csv(self, path: str) -> None:
        QnNetworkAnalysis.write_qn_csv_file(path, self.structure, self)

    def write_network_bonds_csv(self, path: str) -> None:
        QnNetworkAnalysis.write_network_bonds_csv_file(path, self.structure, self)

    def write_report(self, path: str) -> None:
        QnNetworkAnalysis.write_report_file(path, self.structure, self.config, self)

    def write_analysis(
        self,
        prefix: str,
    ) -> None:
        self.write_qn_csv(prefix + "_qn.csv")
        self.write_network_bonds_csv(prefix + "_network_bonds.csv")
        self.write_report(prefix + "_report.txt")

    def write_fnc3(
        self,
        path: str,
        title: Optional[str] = None,
    ) -> None:
        from rmc_dft.rmc_pot.qn.rmc_constraint_writer import RmcConstraintWriter

        writer = RmcConstraintWriter(self)
        writer.write_fnc3(path, title=title)

    def write_unit_topology(
        self,
        path: str,
    ) -> None:
        from rmc_dft.rmc_pot.qn.rmc_constraint_writer import RmcConstraintWriter

        writer = RmcConstraintWriter(self)
        writer.write_unit_topology(path)

    def write_network_topology(
        self,
        path: str,
    ) -> None:
        from rmc_dft.rmc_pot.qn.rmc_constraint_writer import RmcConstraintWriter

        writer = RmcConstraintWriter(self)
        writer.write_network_topology(path)

    def write_all(
        self,
        prefix: str,
    ) -> None:
        from rmc_dft.rmc_pot.qn.rmc_constraint_writer import RmcConstraintWriter

        writer = RmcConstraintWriter(self)
        writer.write_rmc_inputs(prefix)
        self.write_analysis(prefix)

    @staticmethod
    def load_config(
        path: str,
    ) -> Dict[str, Any]:
        """Qn/FNC 生成用 JSON 設定を読む。

        設定ファイルには atom type 順、Qn 定義、FNC=3 の距離拘束、
        TOP 用の bond/angle パラメータを置く。
        """
        with open(path, "r", encoding="utf-8") as f:
            config = json.load(f)
        return config

    @staticmethod
    def get_required(
        config: Dict[str, Any],
        key: str,
    ) -> Any:
        """dict から必須キーを取得する。

        通常の ``dict.__getitem__`` でも取得できるが、エラー時にどの設定キーが
        欠けているかを明示するため、この小関数を通す。
        """
        if key not in config:
            raise ConfigError("config に %s がありません。" % key)
        return config[key]

    @staticmethod
    def get_structure_atom_types(
        config: Dict[str, Any],
    ) -> Optional[List[str]]:
        """POSCAR に元素名行が無い場合に使う atom type 順を取得する。"""
        structure_config = config.get("structure", {})
        if "atom_types" in structure_config:
            return list(structure_config["atom_types"])
        if "atom_types" in config:
            return list(config["atom_types"])
        return None

    @staticmethod
    def get_network_config(
        config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """VASP 構造から Qn network を作るための入力設定を取得する。"""
        if "network" in config:
            return config["network"]
        if "qn" in config:
            return config["qn"]
        raise ConfigError("config に network がありません。")

    @staticmethod
    def _strip_rmc_symbol(symbol: str) -> str:
        """RMC_POT の CHEMICAL-SYMBOLS を VASP の元素ラベルへ寄せる。

        例: ``Si4+`` -> ``Si``, ``O2-`` -> ``O``。
        元素記号として解釈できない文字列はそのまま返す。
        """
        match = re.match(r"^([A-Z][a-z]?)", str(symbol).strip())
        if match is None:
            return str(symbol).strip()
        return match.group(1)

    @staticmethod
    def _parse_free_format_sections(path: str) -> List[Tuple[str, Dict[str, List[List[str]]]]]:
        """RMC_POT free-format ``*.dat`` を section 単位で読む。

        同一 section の複数 instance と、同一 keyword の複数出現を保持する。
        keyword の値が次行へ継続する形式にも対応する。
        """
        sections = []
        current_tag = None
        current_data = None
        current_key = None
        current_values = None

        with open(path, "r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.split("!", 1)[0].strip()
                if not line or line.startswith("#"):
                    continue

                if line.startswith("[") and "]" in line:
                    tag = line[1:line.index("]")].strip().upper()
                    current_tag = tag
                    current_data = {}
                    sections.append((tag, current_data))
                    current_key = None
                    current_values = None
                    continue

                if current_data is None:
                    continue

                if "=" in line:
                    key, value = line.split("=", 1)
                    key = key.strip().upper()
                    values = value.split()
                    current_data.setdefault(key, []).append(values)
                    current_key = key
                    current_values = current_data[key][-1]
                    continue

                # RMC_POT free format では keyword の複数 parameter を次行へ継続できる。
                if current_key is not None and current_values is not None:
                    current_values.extend(line.split())

        return sections

    @staticmethod
    def network_config_from_rmc_dat(
        path: str,
        atom_types: Optional[Sequence[str]] = None,
    ) -> Dict[str, Any]:
        """RMC_POT ``*.dat`` の ``[ SNC ]`` から Qn 幾何条件を作る。

        Qn の目標値、fraction、sigma は「どの Qn をどれだけ要求するか」の
        constraint parameter であり、Qn の数え方そのものには使わない。
        同じ幾何条件で Q0-Q5 など複数の SNC が定義されている場合は一つにまとめる。
        異なる幾何条件の SNC が同じファイルに存在する場合は、自動的に混ぜずに停止する。
        """
        sections = QnNetworkAnalysis._parse_free_format_sections(path)

        chemical_symbols = None
        for tag, data in sections:
            if tag != "GENERAL":
                continue
            values_list = data.get("CHEMICAL-SYMBOLS", [])
            if values_list:
                chemical_symbols = [
                    QnNetworkAnalysis._strip_rmc_symbol(v) for v in values_list[-1]
                ]

        type_labels = chemical_symbols if chemical_symbols is not None else (list(atom_types) if atom_types is not None else None)
        if type_labels is None:
            raise ConfigError(
                "RMC_POT *.dat の atom type index を VASP 元素へ対応付けられません。"
                " CHEMICAL-SYMBOLS または atom_types が必要です。"
            )
        type_labels = [QnNetworkAnalysis._strip_rmc_symbol(v) for v in type_labels]

        geometries = []
        for tag, data in sections:
            if tag != "SNC":
                continue

            ctype_values = data.get("CTYPE_FTYPE_FROM_TO", [])
            stype_values = data.get("STYPE_FROM_TO", [])
            if not ctype_values or not stype_values:
                raise ConfigError(
                    "[ SNC ] には CTYPE_FTYPE_FROM_TO と STYPE_FROM_TO が必要です。"
                )

            cvals = ctype_values[-1]
            if len(cvals) != 4:
                raise ConfigError(
                    "CTYPE_FTYPE_FROM_TO は central_type, first_type, rmin, rmax の4値です。"
                )
            central_type_index = int(cvals[0])
            first_type_index = int(cvals[1])
            first_min = float(cvals[2])
            first_max = float(cvals[3])

            if central_type_index < 1 or central_type_index > len(type_labels):
                raise ConfigError("SNC central type index が atom type 数の範囲外です。")
            if first_type_index < 1 or first_type_index > len(type_labels):
                raise ConfigError("SNC first type index が atom type 数の範囲外です。")

            second_neighbors = []
            for svals in stype_values:
                if len(svals) != 3:
                    raise ConfigError(
                        "STYPE_FROM_TO は second_type, rmin, rmax の3値です。"
                    )
                second_type_index = int(svals[0])
                if second_type_index < 1 or second_type_index > len(type_labels):
                    raise ConfigError("SNC second type index が atom type 数の範囲外です。")
                second_neighbors.append(
                    {
                        "type": type_labels[second_type_index - 1],
                        "min": float(svals[1]),
                        "max": float(svals[2]),
                    }
                )

            # STYPE の順序は Qn 判定結果に影響しないので signature では整列する。
            signature_second = tuple(
                sorted(
                    (x["type"], float(x["min"]), float(x["max"]))
                    for x in second_neighbors
                )
            )
            signature = (
                type_labels[central_type_index - 1],
                type_labels[first_type_index - 1],
                first_min,
                first_max,
                signature_second,
            )
            geometries.append((signature, second_neighbors))

        if not geometries:
            raise ConfigError("RMC_POT *.dat に [ SNC ] がありません。")

        unique = {}
        for signature, second_neighbors in geometries:
            unique.setdefault(signature, second_neighbors)

        if len(unique) != 1:
            raise ConfigError(
                "RMC_POT *.dat に異なる Qn/SNC 幾何条件が複数あります。"
                " 自動解析では一意に選べません。"
            )

        signature, second_neighbors = next(iter(unique.items()))
        central_type, bridge_type, first_min, first_max, _ = signature
        return {
            "structure": {"atom_types": list(type_labels)},
            "network": {
                "central_types": [central_type],
                "bridge_types": [bridge_type],
                "first_min": first_min,
                "first_max": first_max,
                "second_neighbors": second_neighbors,
            },
        }

    @staticmethod
    def get_second_neighbor_specs(qn_config: Dict[str, Any]) -> List[Dict[str, Any]]:
        """RMC_POT SNC の second-neighbour type/range 定義を返す。

        推奨形式は ``second_neighbors`` で、RMC_POT の複数 ``STYPE_FROM_TO`` に
        一対一に対応する。従来設定 ``second_types`` + ``second_min/max`` も、
        同一距離範囲を各 second type に適用する表現として読み込める。
        """
        if "second_neighbors" in qn_config:
            specs = []
            for item in qn_config["second_neighbors"]:
                specs.append(
                    {
                        "type": str(QnNetworkAnalysis.get_required(item, "type")),
                        "min": float(QnNetworkAnalysis.get_required(item, "min")),
                        "max": float(QnNetworkAnalysis.get_required(item, "max")),
                    }
                )
            if not specs:
                raise ConfigError("network.second_neighbors が空です。")
            return specs

        second_types = list(QnNetworkAnalysis.get_required(qn_config, "second_types"))
        if "second_ranges" in qn_config:
            ranges = qn_config["second_ranges"]
            specs = []
            for atom_type in second_types:
                if atom_type not in ranges:
                    raise ConfigError("second_ranges に %s がありません。" % atom_type)
                values = ranges[atom_type]
                if len(values) != 2:
                    raise ConfigError("second_ranges[%s] は [min, max] です。" % atom_type)
                specs.append(
                    {"type": atom_type, "min": float(values[0]), "max": float(values[1])}
                )
            return specs

        second_min = float(QnNetworkAnalysis.get_required(qn_config, "second_min"))
        second_max = float(QnNetworkAnalysis.get_required(qn_config, "second_max"))
        return [
            {"type": atom_type, "min": second_min, "max": second_max}
            for atom_type in second_types
        ]

    @staticmethod
    def validate_rmcpot_qn_config(qn_config: Dict[str, Any]) -> None:
        """RMC_POT SNC に存在しない旧 Qn 判定オプションを拒否する。"""
        forbidden = [
            key
            for key in ("first_mode", "first_nearest_n", "second_use_first_network")
            if key in qn_config
        ]
        if forbidden:
            raise ConfigError(
                "RMC_POT SNC に存在しない Qn 判定オプションがあります: %s. "
                "これらは削除してください。" % ", ".join(forbidden)
            )

    @staticmethod
    def get_rmc_config(
        config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """RMC_POT 出力用設定を取得する。"""
        if "rmc" in config:
            return config["rmc"]
        return config

    @staticmethod
    def get_fnc_config(
        config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """RMC_POT .fnc 出力用設定を取得する。"""
        rmc_config = QnNetworkAnalysis.get_rmc_config(config)
        if "fnc" in rmc_config:
            return rmc_config["fnc"]
        if "fnc3" in config:
            return config["fnc3"]
        raise ConfigError("config に rmc.fnc がありません。")

    @staticmethod
    def get_top_root_config(
        config: Dict[str, Any],
    ) -> Dict[str, Any]:
        """RMC_POT .top 出力用設定の root を取得する。"""
        rmc_config = QnNetworkAnalysis.get_rmc_config(config)
        if "top" in rmc_config:
            return rmc_config["top"]
        return config

    @staticmethod
    def get_top_atom_parameter(
        config: Dict[str, Any],
        label: str,
    ) -> Tuple[float, float]:
        """[ atoms ] 行へ書く charge と mass を atom type ごとに取得する。

        この値は RMC_POT topology の [ atoms ] 列を埋めるための入力である。
        現状の exporter は non-bonded potential を有効化しないため、charge は
        Coulomb 計算の有効化フラグではない。
        """
        top_root = QnNetworkAnalysis.get_top_root_config(config)
        if "atoms" in top_root:
            atoms_config = top_root["atoms"]
            if label not in atoms_config:
                raise ConfigError("rmc.top.atoms に %s がありません。" % label)
            atom_config = atoms_config[label]
            charge = float(QnNetworkAnalysis.get_required(atom_config, "charge"))
            mass = float(QnNetworkAnalysis.get_required(atom_config, "mass"))
            return charge, mass

        masses = config.get("masses", {})
        charges = config.get("charges", {})
        charge = 0.0
        if label in charges:
            charge = float(charges[label])
        mass = 1.0
        if label in masses:
            mass = float(masses[label])
        return charge, mass

    @staticmethod
    def pbc_vector(i: int, j: int, frac: np.ndarray, cell: np.ndarray) -> np.ndarray:
        """PBC 下で i から j への最短イメージベクトルを返す。

        Args:
            i: 始点原子 index。0 始まり。
            j: 終点原子 index。0 始まり。
            frac: fractional coordinate。
            cell: 格子ベクトル行列。単位は Angstrom。
        """
        df = frac[j] - frac[i]
        df = df - np.rint(df)
        return df.dot(cell)

    @staticmethod
    def pbc_distance(i: int, j: int, frac: np.ndarray, cell: np.ndarray) -> float:
        """PBC 下で原子 i-j 間の最短距離を返す。単位は Angstrom。"""
        vec = QnNetworkAnalysis.pbc_vector(i, j, frac, cell)
        return float(np.linalg.norm(vec))

    @staticmethod
    def angle_degree(
        i: int,
        j: int,
        k: int,
        frac: np.ndarray,
        cell: np.ndarray,
    ) -> float:
        """PBC 下で angle i-j-k を degree 単位で返す。

        j が中心原子である。現在の出力では topology の angle パラメータは
        JSON 設定値を使うため、この関数は検査・拡張用である。
        """
        v1 = QnNetworkAnalysis.pbc_vector(j, i, frac, cell)
        v2 = QnNetworkAnalysis.pbc_vector(j, k, frac, cell)
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 <= 0.0 or n2 <= 0.0:
            return 0.0
        cos_value = float(np.dot(v1, v2) / (n1 * n2))
        if cos_value < -1.0:
            cos_value = -1.0
        if cos_value > 1.0:
            cos_value = 1.0
        return float(math.degrees(math.acos(cos_value)))

    @staticmethod
    def pair_key(
        i: int,
        j: int,
    ) -> Tuple[int, int]:
        """原子対を順序に依存しない key に正規化する。

        FNC=3 の距離拘束や topology の bond は i-j と j-i を同一視する。
        そのため、小さい index を先に置いた tuple を key として使う。
        """
        if i < j:
            return (i, j)
        return (j, i)

    @staticmethod
    def select_indices(
        labels: Sequence[str],
        wanted_types: Sequence[str],
    ) -> List[int]:
        """指定した atom type に一致する原子 index を抽出する。"""
        wanted = set(wanted_types)
        result = []
        for i, label in enumerate(labels):
            if label in wanted:
                result.append(i)
        return result

    @staticmethod
    def find_neighbors(
        centers: Sequence[int],
        candidates: Sequence[int],
        poscar: VaspStructure,
        min_distance: float,
        max_distance: float,
        exclude_self: bool,
    ) -> Dict[int, List[int]]:
        """距離範囲で neighbour list を作る。

        Args:
            centers: neighbour を探す中心原子 index 群。
            candidates: neighbour 候補原子 index 群。
            poscar: 構造データ。
            min_distance: 最小距離。単位は Angstrom。
            max_distance: 最大距離。単位は Angstrom。
            exclude_self: True の場合、同一 index は除外する。

        Returns:
            dict: center index -> neighbour index list。
        """
        result = {}
        for center in centers:
            hits = []
            for candidate in candidates:
                if exclude_self:
                    if center == candidate:
                        continue
                distance = QnNetworkAnalysis.pbc_distance(
                    center, candidate, poscar.frac, poscar.cell
                )
                if distance >= min_distance and distance <= max_distance:
                    hits.append(candidate)
            hits.sort()
            result[center] = hits
        return result

    @staticmethod
    def build_network(
        poscar: VaspStructure,
        config: Dict[str, Any],
    ):
        """RMC_POT の SNC/Qn 定義そのままで Qn network を構築する。

        central type ``i`` と first-neighbour type ``j`` の距離が
        ``dmin1 <= r <= dmax1`` に入る first neighbour をすべて対象にする。
        各 first neighbour について、``STYPE_FROM_TO`` に対応する second type
        ``k1..kN`` の少なくとも1原子が指定距離範囲内に存在すれば、その first
        neighbour を Qn に1つ加える。second atom が複数存在しても寄与は1である。
        second type が central type と同じ場合、central atom 自身は除外する。

        パラメータ取得方法（手動設定、RMC_POT *.dat から自動取得など）によって
        この計算ロジックは変化しない。
        """
        qn_config = QnNetworkAnalysis.get_network_config(config)
        QnNetworkAnalysis.validate_rmcpot_qn_config(qn_config)

        central_types = list(QnNetworkAnalysis.get_required(qn_config, "central_types"))
        bridge_types = list(QnNetworkAnalysis.get_required(qn_config, "bridge_types"))
        first_min = float(QnNetworkAnalysis.get_required(qn_config, "first_min"))
        first_max = float(QnNetworkAnalysis.get_required(qn_config, "first_max"))
        second_specs = QnNetworkAnalysis.get_second_neighbor_specs(qn_config)

        # RMC_POT の1個の SNC は central type と first-neighbour type を各1種類持つ。
        if len(central_types) != 1:
            raise ConfigError("RMC_POT SNC 互換の Qn 定義では central_types は1種類です。")
        if len(bridge_types) != 1:
            raise ConfigError("RMC_POT SNC 互換の Qn 定義では bridge_types は1種類です。")
        if first_min > first_max:
            raise ConfigError("first_min は first_max 以下である必要があります。")
        for spec in second_specs:
            if float(spec["min"]) > float(spec["max"]):
                raise ConfigError(
                    "second neighbour %s の min は max 以下である必要があります。"
                    % spec["type"]
                )

        central_indices = QnNetworkAnalysis.select_indices(poscar.labels, central_types)
        bridge_indices = QnNetworkAnalysis.select_indices(poscar.labels, bridge_types)

        if len(central_indices) == 0:
            raise ConfigError("central_types に対応する原子がありません。")
        if len(bridge_indices) == 0:
            raise ConfigError("bridge_types に対応する原子がありません。")

        second_indices_by_spec = []
        for spec in second_specs:
            indices = QnNetworkAnalysis.select_indices(poscar.labels, [spec["type"]])
            if len(indices) == 0:
                raise ConfigError(
                    "second neighbour type %s に対応する原子がありません。" % spec["type"]
                )
            second_indices_by_spec.append((spec, indices))

        data = QnNetworkAnalysis()

        # RMC_POT CTYPE_FTYPE_FROM_TO: 距離範囲に入る first neighbour をすべて使う。
        first_map = QnNetworkAnalysis.find_neighbors(
            central_indices,
            bridge_indices,
            poscar,
            first_min,
            first_max,
            True,
        )

        # RMC_POT STYPE_FROM_TO: second type ごとに独立した距離範囲で探索する。
        # bridge_second_map は全 STYPE の和集合を保持する。
        bridge_second_map = {bridge: [] for bridge in bridge_indices}
        for spec, second_indices in second_indices_by_spec:
            partial_map = QnNetworkAnalysis.find_neighbors(
                bridge_indices,
                second_indices,
                poscar,
                float(spec["min"]),
                float(spec["max"]),
                True,
            )
            for bridge, seconds in partial_map.items():
                bridge_second_map[bridge].extend(seconds)
        for bridge in bridge_second_map:
            bridge_second_map[bridge] = sorted(set(bridge_second_map[bridge]))

        data.first_map = first_map
        data.bridge_second_map = bridge_second_map

        # RMC_POT SNC: 「有効な second neighbour を1個以上持つ first neighbour」の本数が Qn。
        for central in central_indices:
            bridges = first_map[central]
            data.bridges_by_central[central] = list(bridges)
            qn_value = 0
            for bridge in bridges:
                second_hits = [
                    second
                    for second in bridge_second_map.get(bridge, [])
                    if second != central
                ]
                data.second_by_bridge_for_central[(central, bridge)] = second_hits
                if second_hits:
                    qn_value += 1
            data.qn_by_central[central] = qn_value

        # 解析・可視化用 network pair。Qn の数え方には使わない。
        for central, bridges in first_map.items():
            for bridge in bridges:
                data.bond_pairs.add(QnNetworkAnalysis.pair_key(central, bridge))

        active_bridges = set()
        for bridges in first_map.values():
            active_bridges.update(bridges)
        for bridge in active_bridges:
            for second in bridge_second_map.get(bridge, []):
                data.bond_pairs.add(QnNetworkAnalysis.pair_key(bridge, second))

        return data

    @staticmethod
    def add_typed_pair(
        pair_types: Dict[Tuple[int, int], int],
        i: int,
        j: int,
        type_id: int,
    ) -> None:
        """FNC=3 用 typed pair を重複なしで登録する。

        同じ原子対が同じ type_id で複数回出る場合は重複として無視する。
        同じ原子対に異なる type_id が割り当てられた場合は、設定または
        Q network 解析の矛盾として止める。

        FNC の type_id は距離範囲 rmin/rmax と一対一に対応するため、
        type_id の大小で自動的に優先順位を決めてはいけない。
        """
        key = QnNetworkAnalysis.pair_key(i, j)

        if key not in pair_types:
            pair_types[key] = type_id
            return

        old_type_id = pair_types[key]
        if old_type_id == type_id:
            return

        raise ConfigError(
            "同一 FNC pair に異なる type_id が割り当てられました: "
            "atom_i=%d atom_j=%d old_type_id=%d new_type_id=%d"
            % (i + 1, j + 1, old_type_id, type_id)
        )

    @staticmethod
    def build_fnc3_pairs(
        poscar: VaspStructure,
        config: Dict[str, Any],
        data,
    ) -> None:
        """QnNetworkAnalysis から FNC=3 用の原子対リストを作る。

        FNC の unit は TOP の independent unit と同じ assignment を使う。
        JSON の ``rmc.fnc.constraints`` では、各 constraint を ``kind`` と
        ``atom_types`` で明示する。

        対応する ``kind`` は次の通りである。

        * ``center_neighbor``: unit 内の center-neighbor pair。
        * ``neighbor_neighbor``: unit 内の neighbor-neighbor pair。

        各 pair は無向 pair として ``data.fnc_pairs`` に保持する。
        実際の ``*.fnc`` 出力では、全原子の行を書き、各 constrained pair は一度だけ出力する。
        """
        fnc3_config = QnNetworkAnalysis.get_fnc_config(config)
        constraints = QnNetworkAnalysis.get_required(fnc3_config, "constraints")

        from rmc_dft.rmc_pot.qn.rmc_constraint_writer import RmcConstraintWriter

        molecules = RmcConstraintWriter.build_independent_unit_molecules(poscar, config)

        pair_types = {}  # type: Dict[Tuple[int, int], int]
        type_names = {}  # type: Dict[int, str]

        for item in constraints:
            type_id = int(QnNetworkAnalysis.get_required(item, "type_id"))
            name = str(QnNetworkAnalysis.get_required(item, "name"))
            kind = str(QnNetworkAnalysis.get_required(item, "kind"))
            atom_types = QnNetworkAnalysis.get_fnc_constraint_atom_types(item)
            rmin = float(QnNetworkAnalysis.get_required(item, "rmin"))
            rmax = float(QnNetworkAnalysis.get_required(item, "rmax"))

            type_names[type_id] = name

            if kind == "center_neighbor":
                QnNetworkAnalysis.add_unit_center_neighbor_pairs_by_atom_types(
                    poscar, molecules, pair_types, type_id, atom_types, rmin, rmax
                )
            elif kind == "neighbor_neighbor":
                QnNetworkAnalysis.add_unit_neighbor_neighbor_pairs_by_atom_types(
                    poscar, molecules, pair_types, type_id, atom_types, rmin, rmax
                )
            else:
                raise ConfigError("未知の FNC constraint kind です: %s" % kind)

        typed_pairs = []
        for pair in sorted(pair_types.keys()):
            typed_pairs.append((pair[0], pair[1], pair_types[pair]))
        data.fnc_pairs = typed_pairs
        data.fnc_type_names = type_names

    @staticmethod
    def get_fnc_constraint_atom_types(
        item: Dict[str, Any],
    ) -> Tuple[str, str]:
        """FNC constraint の対象 atom type pair を取得する。

        JSON では ``atom_types`` に 2 要素の list として書く。
        例: ``["Si", "O"]``、``["O", "O"]``。
        """
        atom_types = QnNetworkAnalysis.get_required(item, "atom_types")

        if not isinstance(atom_types, list):
            raise ConfigError("FNC constraint atom_types は list で指定してください。")
        if len(atom_types) != 2:
            raise ConfigError("FNC constraint atom_types は 2 要素で指定してください。")
        return (str(atom_types[0]), str(atom_types[1]))

    @staticmethod
    def check_fnc_pair_atom_types(
        poscar: VaspStructure,
        atom_i: int,
        atom_j: int,
        atom_types: Tuple[str, str],
    ) -> bool:
        """FNC pair の実 atom type が JSON の ``atom_types`` と一致するか判定する。"""
        label_i = poscar.labels[atom_i]
        label_j = poscar.labels[atom_j]

        expected = sorted([atom_types[0], atom_types[1]])
        actual = sorted([label_i, label_j])

        if expected == actual:
            return True
        return False

    @staticmethod
    def check_fnc_pair_distance(
        poscar: VaspStructure,
        atom_i: int,
        atom_j: int,
        type_id: int,
        rmin: float,
        rmax: float,
    ) -> None:
        """FNC pair の初期距離が rmin/rmax 範囲内にあることを確認する。"""
        distance = QnNetworkAnalysis.pbc_distance(
            atom_i, atom_j, poscar.frac, poscar.cell
        )

        if distance < rmin or distance > rmax:
            raise ConfigError(
                "FNC pair の初期距離が constraint 範囲外です: "
                "atom_i=%d atom_j=%d pair=%s-%s type_id=%d "
                "distance=%.8f rmin=%.8f rmax=%.8f"
                % (
                    atom_i + 1,
                    atom_j + 1,
                    poscar.labels[atom_i],
                    poscar.labels[atom_j],
                    type_id,
                    distance,
                    rmin,
                    rmax,
                )
            )

    @staticmethod
    def add_checked_fnc_pair(
        poscar: VaspStructure,
        pair_types: Dict[Tuple[int, int], int],
        atom_i: int,
        atom_j: int,
        type_id: int,
        atom_types: Tuple[str, str],
        rmin: float,
        rmax: float,
    ) -> None:
        """atom type と距離範囲を確認して FNC pair を登録する。"""
        if not QnNetworkAnalysis.check_fnc_pair_atom_types(
            poscar, atom_i, atom_j, atom_types
        ):
            raise ConfigError(
                "FNC constraint atom_types と実 atom type が一致しません: "
                "atom_i=%d atom_j=%d actual=%s-%s expected=%s-%s type_id=%d"
                % (
                    atom_i + 1,
                    atom_j + 1,
                    poscar.labels[atom_i],
                    poscar.labels[atom_j],
                    atom_types[0],
                    atom_types[1],
                    type_id,
                )
            )

        QnNetworkAnalysis.check_fnc_pair_distance(
            poscar, atom_i, atom_j, type_id, rmin, rmax
        )
        QnNetworkAnalysis.add_typed_pair(pair_types, atom_i, atom_j, type_id)

    @staticmethod
    def add_unit_center_neighbor_pairs_by_atom_types(
        poscar: VaspStructure,
        molecules: Sequence[
            Tuple[str, List[int], List[Tuple[int, int]], List[Tuple[int, int, int]]]
        ],
        pair_types: Dict[Tuple[int, int], int],
        type_id: int,
        atom_types: Tuple[str, str],
        rmin: float,
        rmax: float,
    ) -> None:
        """unit 内の center-neighbor pair を FNC pair として登録する。

        ``build_independent_unit_molecules`` の unit では local index 0 を center、
        local index 1 以降を neighbor とする。
        """
        for molecule_name, atom_indices, bonds, angles in molecules:
            if len(atom_indices) <= 1:
                continue

            center = atom_indices[0]
            for local_index in range(1, len(atom_indices)):
                neighbor = atom_indices[local_index]
                QnNetworkAnalysis.add_checked_fnc_pair(
                    poscar,
                    pair_types,
                    center,
                    neighbor,
                    type_id,
                    atom_types,
                    rmin,
                    rmax,
                )

    @staticmethod
    def add_unit_neighbor_neighbor_pairs_by_atom_types(
        poscar: VaspStructure,
        molecules: Sequence[
            Tuple[str, List[int], List[Tuple[int, int]], List[Tuple[int, int, int]]]
        ],
        pair_types: Dict[Tuple[int, int], int],
        type_id: int,
        atom_types: Tuple[str, str],
        rmin: float,
        rmax: float,
    ) -> None:
        """unit 内の neighbor-neighbor pair を FNC pair として登録する。"""

        for molecule_name, atom_indices, bonds, angles in molecules:
            if len(atom_indices) <= 2:
                continue

            local_indices = list(range(1, len(atom_indices)))
            for local_i, local_j in itertools.combinations(local_indices, 2):
                atom_i = atom_indices[local_i]
                atom_j = atom_indices[local_j]
                QnNetworkAnalysis.add_checked_fnc_pair(
                    poscar,
                    pair_types,
                    atom_i,
                    atom_j,
                    type_id,
                    atom_types,
                    rmin,
                    rmax,
                )

    @staticmethod
    def collect_fnc3_ranges(config: Dict[str, Any]) -> Dict[int, Tuple[float, float]]:
        """FNC=3 type_id ごとの rmin/rmax を JSON 設定から取り出す。"""
        result = {}
        fnc3_config = QnNetworkAnalysis.get_fnc_config(config)
        constraints = QnNetworkAnalysis.get_required(fnc3_config, "constraints")

        for item in constraints:
            type_id = int(QnNetworkAnalysis.get_required(item, "type_id"))
            rmin = float(QnNetworkAnalysis.get_required(item, "rmin"))
            rmax = float(QnNetworkAnalysis.get_required(item, "rmax"))
            result[type_id] = (rmin, rmax)
        return result

    @staticmethod
    def write_fnc(
        path: str,
        atom_count: int,
        typed_pairs: Sequence[Tuple[int, int, int]],
        ranges: Dict[int, Tuple[float, float]],
        title: str,
    ) -> None:
        """RMC_POT FNC=3 用 ``*.fnc`` ファイルを書く。

        ``typed_pairs`` は 0 始まり index で受け取り、ファイル出力時に 1 始まりへ変換する。
        ``ranges`` は FNC type_id -> (rmin, rmax) で、type_id は 1 から連番である必要がある。
        """
        neighbours = defaultdict(list)  # type: Dict[int, List[Tuple[int, int]]]
        for i, j, type_id in typed_pairs:
            # RMC_POT の .fnc では全 atom について neighbour_count 行を書くが、
            # 各 constrained pair 自体は一度だけ書く。
            # Packmol 由来の既存 FNC と同じく、pair_key で正規化した
            # 小さい atom index 側に neighbour を保持する。
            neighbours[i].append((j, type_id))

        max_type = 0
        for type_id in ranges.keys():
            if type_id > max_type:
                max_type = type_id

        with open(path, "w", encoding="utf-8") as f:
            # RMC_POT の .fnc は、先頭に constraint type ごとの rmin/rmax を書き、
            # その後に各 atom の neighbour index と type_id を並べる。
            f.write("Fixed neighbours constraints (.FNC) file for %s\n" % title)
            f.write("\n")
            f.write(" No. of possible rmin-rmax pairs:\n")
            f.write(" %d\n" % max_type)
            rmin_values = []
            rmax_values = []

            for type_id in range(1, max_type + 1):
                if type_id not in ranges:
                    raise ConfigError("FNC type_id が連番ではありません: %d" % type_id)
                rmin_values.append("%.8f" % ranges[type_id][0])
                rmax_values.append("%.8f" % ranges[type_id][1])

            f.write(" ".join(rmin_values) + "\n")
            f.write(" ".join(rmax_values) + "\n")
            f.write("\n")
            f.write(" %d\n" % atom_count)
            f.write("\n")
            for atom in range(atom_count):
                items = neighbours.get(atom, [])
                f.write("%d %d\n" % (atom + 1, len(items)))
                if len(items) > 0:
                    f.write(" ".join(str(item[0] + 1) for item in items) + "\n")
                    f.write(" ".join(str(item[1]) for item in items) + "\n")

    @staticmethod
    def write_qn_csv_file(
        path: str,
        poscar: VaspStructure,
        data,
    ) -> None:
        """Central atom ごとの Qn 判定結果を CSV へ書く。

        ``central_rmc_index`` と ``bridge_rmc_indices`` は RMC_POT へ対応しやすいように
        1 始まり index で出力する。
        """
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "central_rmc_index",
                    "central_type",
                    "q_n",
                    "first_bridge_count",
                    "bridging_bridge_count",
                    "bridge_rmc_indices",
                    "bridge_types",
                ]
            )
            for central in sorted(data.qn_by_central.keys()):
                bridges = data.bridges_by_central[central]
                bridge_index_text = " ".join(str(bridge + 1) for bridge in bridges)
                bridge_type_text = " ".join(poscar.labels[bridge] for bridge in bridges)
                writer.writerow(
                    [
                        central + 1,
                        poscar.labels[central],
                        data.qn_by_central[central],
                        len(bridges),
                        data.qn_by_central[central],
                        bridge_index_text,
                        bridge_type_text,
                    ]
                )

    @staticmethod
    def write_network_bonds_csv_file(
        path: str,
        poscar: VaspStructure,
        data,
    ) -> None:
        """検出された network bond を CSV へ書く。

        この CSV は FNC=4 topology の [ bonds ] に入る bond の確認用である。
        """
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["i_rmc_index", "j_rmc_index", "i_type", "j_type"])
            for i, j in sorted(data.bond_pairs):
                writer.writerow([i + 1, j + 1, poscar.labels[i], poscar.labels[j]])

    @staticmethod
    def summarize_distances(
        poscar: VaspStructure,
        typed_pairs: Sequence[Tuple[int, int, int]],
    ) -> Dict[int, List[float]]:
        """FNC type_id ごとに現在構造での実距離を集計する。"""
        result = defaultdict(list)  # type: Dict[int, List[float]]
        for i, j, type_id in typed_pairs:
            result[type_id].append(
                QnNetworkAnalysis.pbc_distance(i, j, poscar.frac, poscar.cell)
            )
        return result

    @staticmethod
    def percentile_line(values: Sequence[float]) -> str:
        """距離分布の代表値を report 用の 1 行文字列に変換する。"""
        if len(values) == 0:
            return "count=0"
        arr = np.array(values, dtype=float)
        ps = np.percentile(arr, [0, 1, 5, 50, 95, 99, 100])
        text = (
            "count=%d min=%.4f p1=%.4f p5=%.4f p50=%.4f p95=%.4f p99=%.4f max=%.4f"
            % (len(arr), ps[0], ps[1], ps[2], ps[3], ps[4], ps[5], ps[6])
        )
        return text

    @staticmethod
    def write_report_file(
        path: str,
        poscar: VaspStructure,
        config: Dict[str, Any],
        data,
    ) -> None:
        """Qn 解析と FNC 生成結果の確認用 report を書く。

        report には Qn 定義、Qn 分布、検出された bond/angle 数、
        FNC=3 type ごとの距離統計を含める。
        """
        qn_counts = Counter(data.qn_by_central.values())
        distances_by_type = QnNetworkAnalysis.summarize_distances(
            poscar, data.fnc_pairs
        )

        with open(path, "w", encoding="utf-8") as f:
            f.write("Qn/RMC_POT constraint generation report\n")
            f.write("====================================\n\n")
            f.write("Atom count: %d\n" % poscar.atom_count)
            f.write(
                "Atom types in structure: %s\n" % " ".join(sorted(set(poscar.labels)))
            )
            f.write("Cell volume: %.8f A^3\n" % abs(float(np.linalg.det(poscar.cell))))
            f.write(
                "Number density: %.10f 1/A^3\n\n"
                % (float(poscar.atom_count) / abs(float(np.linalg.det(poscar.cell))))
            )

            f.write("Qn definition\n")
            f.write("-------------\n")
            qn_config = QnNetworkAnalysis.get_network_config(config)
            f.write("central_types: %s\n" % " ".join(qn_config["central_types"]))
            f.write("bridge_types : %s\n" % " ".join(qn_config["bridge_types"]))
            second_specs = QnNetworkAnalysis.get_second_neighbor_specs(qn_config)
            f.write(
                "first range  : %.4f - %.4f A\n"
                % (float(qn_config["first_min"]), float(qn_config["first_max"]))
            )
            for spec in second_specs:
                f.write(
                    "second      : %s %.4f - %.4f A\n"
                    % (spec["type"], float(spec["min"]), float(spec["max"]))
                )
            f.write("\n")

            f.write("Qn counts\n")
            f.write("---------\n")
            for qn_value in sorted(qn_counts.keys()):
                f.write("Q%d: %d\n" % (qn_value, qn_counts[qn_value]))
            f.write("\n")

            f.write("Network counts\n")
            f.write("--------------\n")
            f.write(
                "central-bridge entries: %d\n"
                % sum(len(v) for v in data.first_map.values())
            )
            f.write(
                "bridge-second entries : %d\n"
                % sum(len(v) for v in data.bridge_second_map.values())
            )
            f.write("unique bond pairs      : %d\n" % len(data.bond_pairs))
            f.write("FNC=3 pairs            : %d\n" % len(data.fnc_pairs))
            fnc_error = getattr(data, "fnc_error", None)
            if fnc_error is not None:
                f.write("FNC=3 error            : %s\n" % fnc_error)
            f.write("network topology angles: %d\n\n" % len(data.angles))

            f.write("FNC=3 distance statistics\n")
            f.write("-------------------------\n")
            for type_id in sorted(distances_by_type.keys()):
                type_name = data.fnc_type_names.get(type_id, "type_%d" % type_id)
                f.write(
                    "type %d %s: %s\n"
                    % (
                        type_id,
                        type_name,
                        QnNetworkAnalysis.percentile_line(distances_by_type[type_id]),
                    )
                )

            f.write("\nNotes\n")
            f.write("-----\n")
            f.write("FNC=3 is a distance-window constraint.\n")
            f.write(
                "RMC_POT topology uses bonded potential terms. The .top file is not an FNC file.\n"
            )
            f.write(
                "The standard generated topology is the independent-unit topology.\n"
            )
            f.write(
                "The current-network topology is generated only by an explicit call.\n"
            )
            f.write(
                "For production RMC_POT runs, check the .hst output and acceptance ratio, then tune sigmas/force constants.\n"
            )
