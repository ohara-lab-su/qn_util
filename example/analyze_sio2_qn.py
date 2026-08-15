#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VASP 構造から Qn 分布の個数と比率を表示する。

Qn のパラメータは手動設定または RMC_POT free-format *.dat の [ SNC ] から取得する。
どちらの場合も Qn の数え方は RMC_POT SNC 定義の単一実装を使用する。
CSV、TOP、FNC などのファイルは生成せず、解析結果だけを標準出力する。
"""

from pathlib import Path

from rmc_dft import RmcPot, RmcPotLog

# from rmc_dft.rmc_pot.qn.qn_network_analysis import QnNetworkAnalysis
# from rmc_dft.rmc_pot.qn.vasp_structure import VaspStructure
from qn.qn_network_analysis import QnNetworkAnalysis
from qn.vasp_structure import VaspStructure

from x_logger import XLogger


# =============================================================================
# ユーザー設定
# =============================================================================

# 解析するVASP構造ファイル。
#
# 次の順番で存在を確認する。
#
#     1. POSCAR
#     2. POSCAR.vasp
#     3. aaa.vasp
#     4. aaa.cfg
#
# aaa.cfgだけが存在する場合は、RmcPot.cfg2poscar()を直接使用して
# aaa.vaspへ変換してから解析する。
CONTCAR_PATH = Path("CONTCAR")
POSCAR_PATH = Path("POSCAR")

FALLBACK_POSCAR_PATHS = [
    Path("POSCAR.vasp"),
    Path("aaa.vasp"),
]
CFG_PATH = Path("aaa.cfg")

# aaa.cfgをaaa.vaspへ変換するときに使用する元素順。
ELEMENT_INFO = [
    "Si",
    "O",
]

# Qn パラメータの取得方法。
#
#     "manual"
#         下の QN_CONFIG に書いた値を使う。
#
#     "rmc_dat"
#         実際に RMC_POT で使う free-format *.dat の [ SNC ] から
#         CTYPE_FTYPE_FROM_TO / STYPE_FROM_TO を読み、Qn 距離値を手入力しない。
#
# どちらのモードでも Qn の数え方は完全に同じで、
# QnNetworkAnalysis の RMC_POT SNC 互換ロジックだけを使う。
PARAMETER_MODE = "manual"

# PARAMETER_MODE="rmc_dat" のときに探索する RMC_POT parameter file。
RMC_DAT_PATHS = [
    Path("aaa.dat"),
    Path("aaa.freers"),
    Path("aaa.freer"),
    Path("aaa.free"),
]

# PARAMETER_MODE="manual" のときだけ使う。
# RMC_POT free-format の [ SNC ] で言えば、
#   CTYPE_FTYPE_FROM_TO = central first first_min first_max
#   STYPE_FROM_TO       = second second_min second_max
# に対応する。
#
# Qn の数え方を変更するオプションは存在しない。
QN_CONFIG = {
    "structure": {
        # POSCARに元素名行がない場合に使用する元素順。
        "atom_types": [
            "Si",
            "O",
        ],
    },
    "network": {
        # RMC_POT SNC の central type。1種類だけ指定する。
        "central_types": [
            "Si",
        ],

        # RMC_POT SNC の first-neighbour type。1種類だけ指定する。
        "bridge_types": [
            "O",
        ],

        # central-first の距離範囲 [A]。
        "first_min": 1.2,
        "first_max": 1.9,

        # RMC_POT の STYPE_FROM_TO に対応。複数 second type も指定可能。
        "second_neighbors": [
            {
                "type": "Si",
                "min": 1.2,
                "max": 1.9,
            },
        ],
    },
}


def convert_cfg_to_poscar() -> Path:
    """aaa.cfgをaaa.vaspへ変換する。

    ``conv_cfg2poscar.py aaa Si,O`` と同じ処理を、このスクリプト内で
    直接実行する。

    Returns:
        生成されたaaa.vaspのパス。

    Raises:
        RuntimeError:
            RmcPot.cfg2poscar()の実行後もaaa.vaspが生成されなかった場合。
    """
    cfg_name = CFG_PATH.stem

    logger = XLogger(
        log_level="debug",
        logger_name="analyze_sio2_qn",
    )
    rmc_log = RmcPotLog(
        cfg_name=cfg_name,
        logger=logger,
    )
    rmc = RmcPot(
        rmc_log=rmc_log,
        cfg_name=cfg_name,
        logger=logger,
    )
    rmc.cfg2poscar(
        cfg_name=cfg_name,
        element_info=ELEMENT_INFO,
    )

    generated_poscar_path = Path("%s.vasp" % cfg_name)

    if not generated_poscar_path.is_file():
        raise RuntimeError(
            "%s の変換後に %s が生成されませんでした。"
            % (
                CFG_PATH,
                generated_poscar_path,
            )
        )

    return generated_poscar_path


def resolve_poscar_path() -> Path:
    """解析対象のVASP構造ファイルを決定する。

    次の順番で存在を確認する。

        1. POSCAR
        2. POSCAR.vasp
        3. aaa.vasp
        4. aaa.cfg

    aaa.cfgだけが存在する場合は、aaa.vaspへ変換してそのファイルを返す。

    Returns:
        解析対象となるVASP構造ファイルのパス。

    Raises:
        FileNotFoundError:
            解析可能な入力ファイルが存在しない場合。
    """
    if CONTCAR_PATH.is_file():
        return CONTCAR_PATH

    if POSCAR_PATH.is_file():
        return POSCAR_PATH

    for fallback_poscar_path in FALLBACK_POSCAR_PATHS:
        if fallback_poscar_path.is_file():
            return fallback_poscar_path

    if CFG_PATH.is_file():
        return convert_cfg_to_poscar()

    raise FileNotFoundError(
        "POSCAR, POSCAR.vasp, aaa.vasp, aaa.cfg のいずれも見つかりません。"
    )


def resolve_qn_config():
    """Qn parameter source を解決する。

    parameter source が manual でも rmc_dat でも、返された設定は同じ
    RMC_POT SNC 互換 Qn counting routine に渡される。
    """
    mode = PARAMETER_MODE.strip().lower()

    if mode == "manual":
        return QN_CONFIG, "embedded QN_CONFIG"

    if mode == "rmc_dat":
        dat_path = None
        for candidate in RMC_DAT_PATHS:
            if candidate.is_file():
                dat_path = candidate
                break
        if dat_path is None:
            raise FileNotFoundError(
                "RMC_POT parameter file が見つかりません: %s"
                % ", ".join(str(path) for path in RMC_DAT_PATHS)
            )

        config = QnNetworkAnalysis.network_config_from_rmc_dat(
            str(dat_path),
            atom_types=ELEMENT_INFO,
        )
        return config, "RMC_POT SNC from %s" % dat_path

    raise ValueError(
        "PARAMETER_MODE は 'manual' または 'rmc_dat' です: %s" % PARAMETER_MODE
    )


def main() -> None:
    """VASP 構造を解析し、Qnごとの個数と比率を表示する。"""
    poscar_path = resolve_poscar_path()
    qn_config, qn_config_source = resolve_qn_config()

    atom_types = QnNetworkAnalysis.get_structure_atom_types(qn_config)
    structure = VaspStructure.read_poscar(
        str(poscar_path),
        atom_types,
    )
    analysis = QnNetworkAnalysis(
        structure=structure,
        config=qn_config,
    )

    qn_counts = analysis.get_qn_counts()
    total_count = sum(qn_counts.values())

    print("Qn distribution")
    print("VASP structure:", poscar_path)
    print("Qn config:", qn_config_source)
    print("Total central atoms:", total_count)

    if total_count == 0:
        print("Qn network was not found.")
        return

    max_qn = max(qn_counts.keys())

    print()
    print("Qn    Count    Ratio")
    for qn in range(max_qn + 1):
        count = qn_counts.get(qn, 0)
        ratio = count / total_count
        print("Q%-2d %8d %9.4f %%" % (qn, count, ratio * 100.0))


if __name__ == "__main__":
    main()
