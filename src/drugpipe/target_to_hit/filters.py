"""ADMET / drug-likeness rule filters and final candidate selection."""
# ADMET / 类药性规则过滤与最终候选筛选。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict

import pandas as pd

logger = logging.getLogger(__name__)


class ADMETFilter:
    """
    Apply Lipinski-like rules, QED threshold, structural-alert exclusion,
    and predicted potency cutoff to select hit candidates.
    """
    # 应用类 Lipinski 规则、QED 阈值、结构警示排除与预测效力 cutoff 筛选 hit 候选。

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        flt = cfg.get("target_to_hit", {}).get("filter", {})
        self.pred_pic50_min = float(flt.get("pred_pIC50_min", 6.0))
        self.qed_min = float(flt.get("qed_min", 0.6))
        self.top_n = int(flt.get("top_n", 200))

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

        self.out_dir = out_dir
        self.final_csv = out_dir / "final_hit_candidates.csv"

    # ------------------------------------------------------------------
    def run(self, df_scored: pd.DataFrame) -> pd.DataFrame:
        """Apply all filters and return top-N hit candidates."""
        # 应用全部过滤条件并返回 top-N hit 候选。
        df = df_scored.copy()
        n0 = len(df)

        df = df[df["pred_pIC50_ens"] >= self.pred_pic50_min]
        logger.info("After pIC50 >= %.1f: %d / %d", self.pred_pic50_min, len(df), n0)

        df = df[df["QED"] >= self.qed_min]
        logger.info("After QED >= %.2f: %d", self.qed_min, len(df))

        df = df[~df["HasAlert"]]
        logger.info("After removing structural alerts: %d", len(df))

        mask = self._rules_mask(df)
        df = df[mask]
        logger.info("After ADMET rules: %d", len(df))

        df = df.sort_values(["pred_pIC50_ens", "QED"], ascending=[False, False])
        df = df.head(self.top_n)

        keep_cols = [
            "molecule_chembl_id", "canonical_smiles",
            "pred_pIC50_ens", "pred_IC50_nM_ens", "QED",
            "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "HeavyAtoms",
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]
        df = df[keep_cols].reset_index(drop=True)

        df.to_csv(self.final_csv, index=False)
        logger.info("Hit candidates saved: %s (%d rows)", self.final_csv, len(df))
        return df

    # ------------------------------------------------------------------
    def _rules_mask(self, df: pd.DataFrame) -> pd.Series:
        r = self.rules
        return (
            df["MW"].between(r["mw_min"], r["mw_max"])
            & df["cLogP"].between(r["logp_min"], r["logp_max"])
            & (df["TPSA"] <= r["tpsa_max"])
            & (df["HBD"] <= r["hbd_max"])
            & (df["HBA"] <= r["hba_max"])
            & (df["RotB"] <= r["rotb_max"])
            & (df["Rings"] <= r["rings_max"])
        )
