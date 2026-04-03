"""Virtual screening: score all molecules using trained models."""
# Virtual screening: score all molecules using trained models.
# 说明模块职责、上下游关系与维护注意事项。

# 虚拟筛选：用训练好的模型对所有分子打分。

from __future__ import annotations

import hashlib
import logging
import os
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


def _props_cache_key(smiles_list: List[str]) -> str:
    """Deterministic hash for molecular-property cache."""
    h = hashlib.sha256()
    h.update(f"props_v1,len={len(smiles_list)}\n".encode())
    for smi in smiles_list:
        h.update(smi.encode())
        h.update(b"\n")
    return h.hexdigest()[:16]


class VirtualScreener:
    """Score all molecules with the trained pIC50 model and compute properties."""
    # 用训练好的 pIC50 模型对所有分子打分并计算性质。

    def __init__(
        self,
        cfg: Dict[str, Any],
        out_dir: Path,
        fp_cache_dir: Optional[Path] = None,
        props_cache_dir: Optional[Path] = None,
    ):
        self.out_dir = out_dir
        self.scored_csv = out_dir / "scored_candidates.csv"
        self.featurizer = MorganFeaturizer(cfg, cache_dir=fp_cache_dir)
        self.props_cache_dir = props_cache_dir

    # ------------------------------------------------------------------
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
    def _add_properties(self, df: pd.DataFrame) -> None:
        smiles_list = df["canonical_smiles"].tolist()
        cached = self._load_props_cache(smiles_list)
        if cached is not None:
            for k, v in cached.items():
                df[k] = v
            logger.info("Loaded cached molecular properties (%d molecules).", len(smiles_list))
            return

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

        self._save_props_cache(
            smiles_list,
            {
                "QED": np.asarray(qeds, dtype=np.float32),
                "HasAlert": np.asarray(alerts, dtype=np.bool_),
                **{k: np.asarray(v, dtype=np.float32) for k, v in desc_cols.items()},
            },
        )

    def _props_cache_path(self, smiles_list: List[str]) -> Optional[Path]:
        if self.props_cache_dir is None:
            return None
        key = _props_cache_key(smiles_list)
        return self.props_cache_dir / f"molprops_{key}.npz"

    def _load_props_cache(self, smiles_list: List[str]) -> Optional[Dict[str, np.ndarray]]:
        path = self._props_cache_path(smiles_list)
        if path is None or not path.exists():
            return None
        cols = ["QED", "HasAlert", "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "HeavyAtoms"]
        try:
            data = np.load(path, allow_pickle=False)
            out: Dict[str, np.ndarray] = {}
            for c in cols:
                if c not in data:
                    logger.warning("Property cache missing column %s, recomputing.", c)
                    return None
                arr = np.array(data[c])
                if arr.shape[0] != len(smiles_list):
                    logger.warning("Property cache length mismatch, recomputing.")
                    return None
                out[c] = arr
            return out
        except Exception as exc:
            logger.warning("Failed to load property cache %s: %s", path.name, exc)
            return None

    def _save_props_cache(self, smiles_list: List[str], cols: Dict[str, np.ndarray]) -> None:
        path = self._props_cache_path(smiles_list)
        if path is None:
            return
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, **cols)
            size_mb = path.stat().st_size / (1024 * 1024)
            logger.info("Saved molecular-property cache: %s (%.1f MB)", path.name, size_mb)
        except Exception as exc:
            logger.warning("Failed to save property cache: %s", exc)


def _compute_props_single(smi: str) -> Tuple[float, bool, Dict[str, Any]]:
    """Compute QED, structural alerts, and descriptors for one SMILES.

    Top-level function so it can be pickled by multiprocessing.
    """
    mol = safe_mol(smi)
    if mol is None:
        return (np.nan, True, {})
    return (calc_qed(mol), has_structural_alert(mol), calc_descriptors(mol))
