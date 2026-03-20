"""Virtual screening: score all molecules using trained models."""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 虚拟筛选：用训练好的模型对所有分子打分。

from __future__ import annotations

import logging
import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from drugpipe.target_to_hit.featurizer import MorganFeaturizer
from drugpipe.target_to_hit.models import ModelTrainer
from drugpipe.utils.chem import (
    calc_descriptors,
    calc_qed,
    has_structural_alert,
    pic50_to_ic50,
    safe_mol,
)

logger = logging.getLogger(__name__)


# EN: VirtualScreener core behavior and intent.
# 中文：VirtualScreener 的核心行为与设计意图。
class VirtualScreener:
    """Score all molecules with the trained pIC50 model and compute properties."""
    # 用训练好的 pIC50 模型对所有分子打分并计算性质。

    # EN: __init__ core behavior and intent.
    # 中文：__init__ 的核心行为与设计意图。
    def __init__(self, cfg: Dict[str, Any], out_dir: Path, fp_cache_dir: Optional[Path] = None):
        self.out_dir = out_dir
        self.scored_csv = out_dir / "scored_candidates.csv"
        self.featurizer = MorganFeaturizer(cfg, cache_dir=fp_cache_dir)

    # ------------------------------------------------------------------
    # EN: run core behavior and intent.
    # 中文：run 的核心行为与设计意图。
    def run(
        self,
        df_mol: pd.DataFrame,
        trainer: ModelTrainer,
    ) -> pd.DataFrame:
        """
        Predict pIC50 for every molecule in *df_mol*, compute QED /
        descriptors / structural alerts, and save to CSV.

        If ``scored_candidates.csv`` already exists with the correct row
        count, loads it directly instead of recomputing (~15 min saved).
        """
        if self.scored_csv.exists():
            cached = pd.read_csv(self.scored_csv)
            if len(cached) == len(df_mol) and "pred_pIC50_ens" in cached.columns:
                logger.info(
                    "Loaded cached scored candidates: %s (%d rows)",
                    self.scored_csv, len(cached),
                )
                return cached

        smiles = df_mol["canonical_smiles"].tolist()
        logger.info("Featurizing %d molecules ...", len(smiles))
        X = self.featurizer.transform(smiles)

        logger.info("Predicting pIC50 (ensemble) ...")
        pred_ens = trainer.predict(X)
        pred_rf = trainer.predict_rf(X)

        out = df_mol.copy()
        out["pred_pIC50_rf"] = pred_rf
        out["pred_pIC50_ens"] = pred_ens
        out["pred_IC50_nM_ens"] = [pic50_to_ic50(p) for p in pred_ens]

        pred_mlp = trainer.predict_mlp(X)
        if pred_mlp is not None:
            out["pred_pIC50_mlp"] = pred_mlp

        logger.info("Computing molecular properties ...")
        self._add_properties(out)

        out.to_csv(self.scored_csv, index=False)
        logger.info("Scored candidates saved: %s (%d rows)", self.scored_csv, len(out))
        return out

    # ------------------------------------------------------------------
    @staticmethod
    # EN: _add_properties core behavior and intent.
    # 中文：_add_properties 的核心行为与设计意图。
    def _add_properties(df: pd.DataFrame) -> None:
        smiles_list = df["canonical_smiles"].tolist()
        n_workers = min(os.cpu_count() or 1, 16)

        logger.info("Computing molecular properties (%d molecules, %d workers) ...",
                     len(smiles_list), n_workers)

        if n_workers > 1 and len(smiles_list) > 500:
            chunk_size = max(200, len(smiles_list) // (n_workers * 4))
            with Pool(n_workers) as pool:
                results = pool.map(_compute_props_single, smiles_list, chunksize=chunk_size)
        else:
            results = [_compute_props_single(smi) for smi in smiles_list]

        _DESC_KEYS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "HeavyAtoms"]
        qeds, alerts = [], []
        desc_cols: Dict[str, list] = {k: [] for k in _DESC_KEYS}

        for q, a, desc in results:
            qeds.append(q)
            alerts.append(a)
            for k in _DESC_KEYS:
                desc_cols[k].append(desc.get(k, np.nan))

        df["QED"] = qeds
        df["HasAlert"] = alerts
        for k, v in desc_cols.items():
            df[k] = v


# EN: _compute_props_single core behavior and intent.
# 中文：_compute_props_single 的核心行为与设计意图。
def _compute_props_single(smi: str) -> Tuple[float, bool, Dict[str, Any]]:
    """Compute QED, structural alerts, and descriptors for one SMILES.

    Top-level function so it can be pickled by multiprocessing.
    """
    mol = safe_mol(smi)
    if mol is None:
        return (np.nan, True, {})
    return (calc_qed(mol), has_structural_alert(mol), calc_descriptors(mol))
