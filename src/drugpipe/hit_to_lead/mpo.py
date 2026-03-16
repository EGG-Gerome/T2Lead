"""Multi-Parameter Optimization (MPO) scoring for lead prioritization."""
# 多参数优化（MPO）打分，用于先导化合物排序。

from __future__ import annotations

import logging
from typing import Any, Dict

import numpy as np
import pandas as pd

from drugpipe.utils.chem import (
    calc_descriptors,
    calc_qed,
    has_structural_alert,
    safe_mol,
)

logger = logging.getLogger(__name__)


class MPOScorer:
    """
    Score candidate molecules using a weighted combination of:
      - predicted potency  (pIC50)
      - drug-likeness      (QED)
      - ADMET rule pass    (binary → 0/1)
      - structural novelty (Tanimoto distance to nearest known hit)
    """
    # 用加权组合对候选分子打分：预测效力(pIC50)、类药性(QED)、ADMET 通过(0/1)、结构新颖性。

    def __init__(self, cfg: Dict[str, Any]):
        h2l = cfg.get("hit_to_lead", {})
        mpo = h2l.get("mpo", {})
        self.w_potency = float(mpo.get("w_potency", 0.40))
        self.w_qed = float(mpo.get("w_qed", 0.25))
        self.w_admet = float(mpo.get("w_admet", 0.20))
        self.w_novelty = float(mpo.get("w_novelty", 0.15))

        flt = cfg.get("target_to_hit", {}).get("filter", {})
        self.rules = {
            "mw_min": float(flt.get("mw_min", 150)),
            "mw_max": float(flt.get("mw_max", 500)),
            "logp_min": float(flt.get("logp_min", -1)),
            "logp_max": float(flt.get("logp_max", 5)),
            "tpsa_max": float(flt.get("tpsa_max", 140)),
            "hbd_max": int(flt.get("hbd_max", 5)),
            "hba_max": int(flt.get("hba_max", 10)),
            "rotb_max": int(flt.get("rotb_max", 10)),
            "rings_max": int(flt.get("rings_max", 6)),
        }

    # ------------------------------------------------------------------
    def score(
        self,
        df: pd.DataFrame,
        model_predict_fn=None,
        featurizer_fn=None,
    ) -> pd.DataFrame:
        """
        Compute MPO score for every row in *df*.

        If *model_predict_fn* and *featurizer_fn* are given, predicted pIC50
        will be computed for new analogs that lack it.  Otherwise only rows
        with existing ``pred_pIC50_ens`` are scored.
        """
        # 若提供 model_predict_fn 与 featurizer_fn，会为新类似物补算 pIC50；否则仅对已有 pred_pIC50_ens 的行打分。
        df = df.copy()

        needs_pred = df["pred_pIC50_ens"].isna() if "pred_pIC50_ens" in df.columns else pd.Series([True] * len(df))
        if needs_pred.any() and model_predict_fn is not None and featurizer_fn is not None:
            idx = df.index[needs_pred]
            smi_list = df.loc[idx, "canonical_smiles"].tolist()
            X = featurizer_fn(smi_list)
            preds = model_predict_fn(X)
            df.loc[idx, "pred_pIC50_ens"] = preds

        self._ensure_properties(df)

        potency_norm = self._normalize(df["pred_pIC50_ens"].fillna(0))
        qed_norm = df["QED"].fillna(0).clip(0, 1)
        admet_norm = df["admet_pass"].astype(float)

        novelty_norm = pd.Series(np.ones(len(df)), index=df.index)
        if "source_smiles" in df.columns:
            novelty_norm = (df["canonical_smiles"] != df["source_smiles"]).astype(float)

        df["mpo_score"] = (
            self.w_potency * potency_norm
            + self.w_qed * qed_norm
            + self.w_admet * admet_norm
            + self.w_novelty * novelty_norm
        )

        return df

    # ------------------------------------------------------------------
    def _ensure_properties(self, df: pd.DataFrame) -> None:
        """Compute QED / descriptors / ADMET pass for rows missing them."""
        # 对缺失性质的行计算 QED、描述符与 ADMET 通过情况。
        if "QED" not in df.columns:
            df["QED"] = np.nan
        if "admet_pass" not in df.columns:
            df["admet_pass"] = False
        if "HasAlert" not in df.columns:
            df["HasAlert"] = True

        need_calc = df["QED"].isna()
        for idx in df.index[need_calc]:
            smi = df.at[idx, "canonical_smiles"]
            mol = safe_mol(smi)
            if mol is None:
                df.at[idx, "QED"] = 0.0
                df.at[idx, "HasAlert"] = True
                df.at[idx, "admet_pass"] = False
                continue
            df.at[idx, "QED"] = calc_qed(mol)
            df.at[idx, "HasAlert"] = has_structural_alert(mol)
            desc = calc_descriptors(mol)
            df.at[idx, "admet_pass"] = self._rules_pass(desc) and not df.at[idx, "HasAlert"]

        already = ~need_calc
        if already.any():
            _desc_cols = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "HeavyAtoms"]
            has_desc = all(c in df.columns for c in _desc_cols)
            for idx in df.index[already]:
                if has_desc and pd.notna(df.at[idx, "MW"]):
                    desc = {c: df.at[idx, c] for c in _desc_cols}
                else:
                    mol = safe_mol(df.at[idx, "canonical_smiles"])
                    if mol is None:
                        df.at[idx, "admet_pass"] = False
                        continue
                    desc = calc_descriptors(mol)
                df.at[idx, "admet_pass"] = self._rules_pass(desc) and not df.at[idx, "HasAlert"]

    def _rules_pass(self, desc: Dict[str, Any]) -> bool:
        """Check if descriptor dict satisfies ADMET rules. / 检查描述符是否满足 ADMET 规则。"""
        r = self.rules
        return (
            r["mw_min"] <= desc.get("MW", 0) <= r["mw_max"]
            and r["logp_min"] <= desc.get("cLogP", 0) <= r["logp_max"]
            and desc.get("TPSA", 999) <= r["tpsa_max"]
            and desc.get("HBD", 99) <= r["hbd_max"]
            and desc.get("HBA", 99) <= r["hba_max"]
            and desc.get("RotB", 99) <= r["rotb_max"]
            and desc.get("Rings", 99) <= r["rings_max"]
        )

    @staticmethod
    def _normalize(s: pd.Series) -> pd.Series:
        """Min-max normalize to [0, 1]."""
        # 最小-最大归一化到 [0, 1]。
        mn, mx = s.min(), s.max()
        if mx - mn < 1e-9:
            return pd.Series(np.ones(len(s)), index=s.index)
        return (s - mn) / (mx - mn)
