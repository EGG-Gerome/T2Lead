"""Virtual screening: score all molecules using trained models."""
# 虚拟筛选：用训练好的模型对所有分子打分。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

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


class VirtualScreener:
    """Score all molecules with the trained pIC50 model and compute properties."""
    # 用训练好的 pIC50 模型对所有分子打分并计算性质。

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        self.out_dir = out_dir
        self.scored_csv = out_dir / "scored_candidates.csv"
        self.featurizer = MorganFeaturizer(cfg)

    # ------------------------------------------------------------------
    def run(
        self,
        df_mol: pd.DataFrame,
        trainer: ModelTrainer,
    ) -> pd.DataFrame:
        """
        Predict pIC50 for every molecule in *df_mol*, compute QED /
        descriptors / structural alerts, and save to CSV.
        """
        # 对 df_mol 中每个分子预测 pIC50，计算 QED/描述符/结构警示并保存为 CSV。
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
    def _add_properties(df: pd.DataFrame) -> None:
        qeds, alerts, pass_rules_list = [], [], []
        desc_cols: Dict[str, list] = {
            k: [] for k in ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "HeavyAtoms"]
        }

        for smi in df["canonical_smiles"]:
            mol = safe_mol(smi)
            if mol is None:
                qeds.append(np.nan)
                alerts.append(True)
                pass_rules_list.append(False)
                for k in desc_cols:
                    desc_cols[k].append(np.nan)
                continue

            qeds.append(calc_qed(mol))
            alerts.append(has_structural_alert(mol))
            desc = calc_descriptors(mol)
            for k in desc_cols:
                desc_cols[k].append(float(desc[k]))
            pass_rules_list.append(True)  # actual rule check deferred to filter

        df["QED"] = qeds
        df["HasAlert"] = alerts
        for k, v in desc_cols.items():
            df[k] = v
