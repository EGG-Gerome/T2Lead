"""Murcko scaffold extraction and structure-activity relationship (SAR) analysis."""
# Murcko 骨架提取与构效关系（SAR）分析。

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any, Dict, List, Tuple

import pandas as pd

from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import (
    GetScaffoldForMol,
    MakeScaffoldGeneric,
)

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)


class ScaffoldAnalyzer:
    """
    Extract Murcko scaffolds, group hits by scaffold, and produce a
    per-scaffold SAR summary.
    """
    # 提取 Murcko 骨架，按骨架分组 hit，并生成每骨架 SAR 汇总。

    def __init__(self, cfg: Dict[str, Any]):
        h2l = cfg.get("hit_to_lead", {})
        self.min_cluster = int(h2l.get("scaffold", {}).get("min_cluster_size", 3))

    # ------------------------------------------------------------------
    def analyze(self, df_hits: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parameters
        ----------
        df_hits : DataFrame with at least ``canonical_smiles`` and ``pred_pIC50_ens``.

        Returns
        -------
        df_hits : original frame augmented with ``scaffold_smi`` and ``generic_scaffold_smi``.
        df_scaffold_summary : one row per scaffold with SAR statistics.
        """
        # 返回：带 scaffold_smi、generic_scaffold_smi 的 df_hits；每骨架一行的 SAR 汇总表。
        scaffolds: List[str] = []
        generic: List[str] = []

        for smi in df_hits["canonical_smiles"]:
            mol = safe_mol(smi)
            if mol is None:
                scaffolds.append("")
                generic.append("")
                continue
            try:
                core = GetScaffoldForMol(mol)
                scaffolds.append(Chem.MolToSmiles(core))
                gen = MakeScaffoldGeneric(core)
                generic.append(Chem.MolToSmiles(gen))
            except Exception:
                scaffolds.append("")
                generic.append("")

        df = df_hits.copy()
        df["scaffold_smi"] = scaffolds
        df["generic_scaffold_smi"] = generic

        summary = self._summarize(df)

        logger.info(
            "Scaffold analysis: %d unique scaffolds, %d with >= %d members",
            len(summary),
            (summary["count"] >= self.min_cluster).sum(),
            self.min_cluster,
        )
        return df, summary

    # ------------------------------------------------------------------
    def _summarize(self, df: pd.DataFrame) -> pd.DataFrame:
        groups = df.groupby("generic_scaffold_smi")
        rows = []
        for scaf, grp in groups:
            if not scaf:
                continue
            pic50 = grp["pred_pIC50_ens"]
            rows.append({
                "generic_scaffold_smi": scaf,
                "count": len(grp),
                "pIC50_mean": pic50.mean(),
                "pIC50_std": pic50.std(),
                "pIC50_max": pic50.max(),
                "example_smiles": grp["canonical_smiles"].iloc[0],
            })
        summary = pd.DataFrame(rows)
        if not summary.empty:
            summary = summary.sort_values("count", ascending=False).reset_index(drop=True)
        return summary
