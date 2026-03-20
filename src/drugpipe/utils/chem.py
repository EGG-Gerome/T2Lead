"""RDKit chemistry utilities shared across pipeline stages."""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 各流水线阶段共用的 RDKit 化学工具。

from __future__ import annotations

import math
from typing import Any, Dict, Optional

import numpy as np

from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.QED import qed
from rdkit.Chem import rdFingerprintGenerator


# ---------------------------------------------------------------------------
# Singleton filter catalog (PAINS + Brenk structural alerts)
# 单例结构警示目录（PAINS + Brenk）
# ---------------------------------------------------------------------------

# EN: _build_filter_catalog core behavior and intent.
# 中文：_build_filter_catalog 的核心行为与设计意图。
def _build_filter_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


_FILTER_CATALOG: Optional[FilterCatalog] = None


# EN: get_filter_catalog core behavior and intent.
# 中文：get_filter_catalog 的核心行为与设计意图。
def get_filter_catalog() -> FilterCatalog:
    """Return the shared PAINS/Brenk filter catalog. / 返回共用的 PAINS/Brenk 过滤目录。"""
    global _FILTER_CATALOG
    if _FILTER_CATALOG is None:
        _FILTER_CATALOG = _build_filter_catalog()
    return _FILTER_CATALOG


# ---------------------------------------------------------------------------
# SMILES → Mol
# SMILES 转 Mol 对象
# ---------------------------------------------------------------------------

# EN: safe_mol core behavior and intent.
# 中文：safe_mol 的核心行为与设计意图。
def safe_mol(smiles: str) -> Optional[Chem.Mol]:
    """Parse SMILES, returning *None* on any failure."""
    # 解析 SMILES，失败时返回 None。
    if not smiles or not isinstance(smiles, str):
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Fingerprints
# 分子指纹
# ---------------------------------------------------------------------------

_FP_GENERATORS: Dict[tuple, Any] = {}


# EN: _get_morgan_gen core behavior and intent.
# 中文：_get_morgan_gen 的核心行为与设计意图。
def _get_morgan_gen(radius: int, nbits: int) -> Any:
    key = (radius, nbits)
    if key not in _FP_GENERATORS:
        _FP_GENERATORS[key] = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=nbits,
        )
    return _FP_GENERATORS[key]


# EN: morgan_fp_bits core behavior and intent.
# 中文：morgan_fp_bits 的核心行为与设计意图。
def morgan_fp_bits(mol: Chem.Mol, radius: int = 2, nbits: int = 2048) -> np.ndarray:
    """Morgan (ECFP) bit vector. / Morgan（ECFP）位向量。"""
    gen = _get_morgan_gen(radius, nbits)
    fp = gen.GetFingerprintAsNumPy(mol)
    return fp.astype(np.uint8)


# EN: smiles_to_fp core behavior and intent.
# 中文：smiles_to_fp 的核心行为与设计意图。
def smiles_to_fp(smiles: str, radius: int = 2, nbits: int = 2048) -> Optional[np.ndarray]:
    mol = safe_mol(smiles)
    if mol is None:
        return None
    return morgan_fp_bits(mol, radius, nbits)


# EN: batch_fps core behavior and intent.
# 中文：batch_fps 的核心行为与设计意图。
def batch_fps(smiles_list: list[str], radius: int = 2, nbits: int = 2048) -> np.ndarray:
    """Return (N, nbits) float32 fingerprint matrix."""
    # 返回 (N, nbits) 的 float32 指纹矩阵。
    X = np.zeros((len(smiles_list), nbits), dtype=np.float32)
    for i, smi in enumerate(smiles_list):
        mol = safe_mol(smi)
        if mol is not None:
            X[i] = morgan_fp_bits(mol, radius, nbits).astype(np.float32)
    return X


# ---------------------------------------------------------------------------
# IC50 ↔ pIC50 conversion
# IC50 与 pIC50 互算
# ---------------------------------------------------------------------------

# EN: ic50_to_pic50 core behavior and intent.
# 中文：ic50_to_pic50 的核心行为与设计意图。
def ic50_to_pic50(ic50_nM: float) -> Optional[float]:
    """pIC50 = 9 − log10(IC50 in nM)."""
    # pIC50 = 9 − log10(IC50 [nM])。
    try:
        v = float(ic50_nM)
        if v <= 0 or math.isnan(v) or math.isinf(v):
            return None
        return 9.0 - math.log10(v)
    except Exception:
        return None


# EN: pic50_to_ic50 core behavior and intent.
# 中文：pic50_to_ic50 的核心行为与设计意图。
def pic50_to_ic50(pic50: float) -> float:
    """IC50 (nM) = 10^(9 − pIC50)."""
    # IC50 (nM) = 10^(9 − pIC50)。
    return 10.0 ** (9.0 - float(pic50))


# ---------------------------------------------------------------------------
# Molecular descriptors
# 分子描述符
# ---------------------------------------------------------------------------

# EN: calc_descriptors core behavior and intent.
# 中文：calc_descriptors 的核心行为与设计意图。
def calc_descriptors(mol: Chem.Mol) -> Dict[str, Any]:
    """Compute MW, cLogP, TPSA, HBD, HBA, RotB, Rings, HeavyAtoms."""
    # 计算分子量、cLogP、TPSA、氢键供体/受体、可旋转键、环数、重原子数。
    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
        "Rings": Lipinski.RingCount(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
    }


# EN: calc_qed core behavior and intent.
# 中文：calc_qed 的核心行为与设计意图。
def calc_qed(mol: Chem.Mol) -> float:
    """Quantitative Estimate of Drug-likeness. / 药物相似性定量估计（QED）。"""
    return float(qed(mol))


# EN: has_structural_alert core behavior and intent.
# 中文：has_structural_alert 的核心行为与设计意图。
def has_structural_alert(mol: Chem.Mol) -> bool:
    """Return True if the molecule matches any PAINS / Brenk alert."""
    # 若分子命中任一 PAINS/Brenk 警示则返回 True。
    return get_filter_catalog().HasMatch(mol)
