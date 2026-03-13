"""Butina clustering for chemical diversity analysis of hit compounds."""
# Butina 聚类，用于 hit 化合物的化学多样性分析。

from __future__ import annotations

import logging
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)


class DiversityClusterer:
    """
    Cluster hit molecules using Butina algorithm on Tanimoto distance,
    then pick a representative (centroid) per cluster.
    """
    # 基于 Tanimoto 距离用 Butina 算法对 hit 分子聚类，每类取一代表（质心）。

    def __init__(self, cfg: Dict[str, Any]):
        h2l = cfg.get("hit_to_lead", {})
        self.cutoff = float(h2l.get("clustering", {}).get("tanimoto_cutoff", 0.4))

    # ------------------------------------------------------------------
    def cluster(self, df_hits: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parameters
        ----------
        df_hits : DataFrame with ``canonical_smiles`` column.

        Returns
        -------
        df_hits   : augmented with ``cluster_id`` and ``is_centroid``.
        df_summary: one row per cluster with size + centroid SMILES.
        """
        # 返回：带 cluster_id、is_centroid 的 df_hits；每类一行（含 size、centroid SMILES）的汇总。
        smiles_list = df_hits["canonical_smiles"].tolist()
        fps = self._compute_fps(smiles_list)

        valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
        valid_fps = [fps[i] for i in valid_idx]

        if len(valid_fps) < 2:
            logger.warning("Too few valid molecules (%d) for clustering.", len(valid_fps))
            df = df_hits.copy()
            df["cluster_id"] = 0
            df["is_centroid"] = True
            return df, pd.DataFrame()

        dists = self._tanimoto_distance_matrix(valid_fps)
        clusters = Butina.ClusterData(dists, len(valid_fps), self.cutoff, isDistData=True)

        cluster_map = {}
        centroid_set = set()
        for cid, members in enumerate(clusters):
            centroid_set.add(valid_idx[members[0]])
            for m in members:
                cluster_map[valid_idx[m]] = cid

        df = df_hits.copy()
        df["cluster_id"] = [cluster_map.get(i, -1) for i in range(len(df))]
        df["is_centroid"] = [i in centroid_set for i in range(len(df))]

        summary_rows = []
        for cid, members in enumerate(clusters):
            centroid_i = valid_idx[members[0]]
            summary_rows.append({
                "cluster_id": cid,
                "size": len(members),
                "centroid_smiles": smiles_list[centroid_i],
            })
        summary = pd.DataFrame(summary_rows).sort_values("size", ascending=False).reset_index(drop=True)

        logger.info(
            "Butina clustering (cutoff=%.2f): %d clusters from %d molecules, "
            "largest cluster size=%d",
            self.cutoff, len(clusters), len(valid_fps),
            summary["size"].max() if not summary.empty else 0,
        )
        return df, summary

    # ------------------------------------------------------------------
    @staticmethod
    def _compute_fps(smiles_list: List[str]) -> List:
        fps = []
        for smi in smiles_list:
            mol = safe_mol(smi)
            if mol is None:
                fps.append(None)
            else:
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
        return fps

    @staticmethod
    def _tanimoto_distance_matrix(fps: list) -> List[float]:
        """Flat lower-triangular distance list required by Butina."""
        # Butina 所需的扁平下三角距离列表。
        n = len(fps)
        dists = []
        for i in range(1, n):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend([1.0 - s for s in sims])
        return dists
