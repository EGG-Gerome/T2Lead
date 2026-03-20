"""Build a training dataset for a specific target from crawled ChEMBL data."""
# EN: Build a training dataset for a specific target from crawled ChEMBL data.
# 中文：说明模块职责、上下游关系与维护注意事项。

# 从爬取的 ChEMBL 数据为指定靶点构建训练数据集。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from drugpipe.utils.chem import (
    calc_descriptors,
    calc_qed,
    has_structural_alert,
    ic50_to_pic50,
    safe_mol,
)

logger = logging.getLogger(__name__)


class DatasetBuilder:
    """
    Join molecules + activities for a given target, compute pIC50,
    descriptors, QED, and structural alerts.
    """
    # 合并指定靶点的分子与活性，计算 pIC50、描述符、QED 与结构警示。

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        t2h = cfg.get("target_to_hit", {})
        self.min_train = int(t2h.get("dataset", {}).get("min_train_samples", 200))
        self.out_dir = out_dir
        self.dataset_csv = out_dir / "dataset_target_ic50.csv"

    # ------------------------------------------------------------------
    def build(
        self,
        molecules_csv: Path,
        activities_csv: Path,
        target_chembl_id: Optional[str] = None,
    ) -> pd.DataFrame:
        """Build target-specific IC50 dataset with pIC50 and properties. / 构建靶点专用 IC50 数据集（含 pIC50 与性质）。"""
        df_mol = self._load_molecules(molecules_csv)
        df_act = self._load_activities(activities_csv)

        target_id = target_chembl_id or self._auto_select_target(df_act)
        logger.info("Building dataset for target %s", target_id)

        df_t = df_act[df_act["target_chembl_id"] == target_id].copy()
        if df_t.empty:
            raise RuntimeError(f"No IC50 data found for target {target_id}")

        df_t["pIC50"] = df_t["standard_value"].apply(ic50_to_pic50)
        df_t = df_t.dropna(subset=["pIC50"])

        agg = (
            df_t.groupby("molecule_chembl_id", as_index=False)
            .agg(
                ic50_nM_median=("standard_value", "median"),
                pIC50_median=("pIC50", "median"),
                n_measurements=("pIC50", "count"),
            )
        )

        df = agg.merge(df_mol, on="molecule_chembl_id", how="inner")
        df = df.dropna(subset=["canonical_smiles"]).drop_duplicates("molecule_chembl_id")

        df = self._add_mol_properties(df)

        if len(df) < self.min_train:
            raise RuntimeError(
                f"Only {len(df)} usable samples for target {target_id} "
                f"(minimum {self.min_train}). Try a different target or relax constraints."
            )

        df.to_csv(self.dataset_csv, index=False)
        logger.info("Dataset saved: %s (%d molecules)", self.dataset_csv, len(df))
        return df

    # ------------------------------------------------------------------
    @staticmethod
    def _load_molecules(path: Path) -> pd.DataFrame:
        df = pd.read_csv(path)
        df = df.dropna(subset=["molecule_chembl_id", "canonical_smiles"])
        return df.drop_duplicates("molecule_chembl_id")

    @staticmethod
    def _load_activities(path: Path) -> pd.DataFrame:
        df = pd.read_csv(path)
        return df.dropna(subset=["molecule_chembl_id", "target_chembl_id", "standard_value"])

    @staticmethod
    def _auto_select_target(df_act: pd.DataFrame) -> str:
        """Pick the target with most IC50 records. / 选择 IC50 记录数最多的靶点。"""
        vc = df_act["target_chembl_id"].value_counts()
        if vc.empty:
            raise RuntimeError("Activity data is empty — cannot auto-select target.")
        target = str(vc.index[0])
        logger.info("Auto-selected target %s (%d records)", target, vc.iloc[0])
        return target

    @staticmethod
    def _add_mol_properties(df: pd.DataFrame) -> pd.DataFrame:
        """Compute QED, descriptors, and structural alerts for valid SMILES."""
        # 对有效 SMILES 计算 QED、描述符与结构警示。
        valid_mask = []
        qeds, alerts = [], []
        desc_cols: Dict[str, list] = {k: [] for k in [
            "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "HeavyAtoms",
        ]}

        for smi in df["canonical_smiles"]:
            mol = safe_mol(smi)
            if mol is None:
                valid_mask.append(False)
                qeds.append(np.nan)
                alerts.append(True)
                for k in desc_cols:
                    desc_cols[k].append(np.nan)
                continue

            valid_mask.append(True)
            qeds.append(calc_qed(mol))
            alerts.append(has_structural_alert(mol))
            desc = calc_descriptors(mol)
            for k in desc_cols:
                desc_cols[k].append(float(desc[k]))

        df = df.copy()
        df["QED"] = qeds
        df["HasAlert"] = alerts
        for k, v in desc_cols.items():
            df[k] = v

        return df[valid_mask].reset_index(drop=True)
