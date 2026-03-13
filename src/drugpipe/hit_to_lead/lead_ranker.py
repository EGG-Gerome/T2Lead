"""Final lead ranking: combine all candidates, filter, rank by MPO, and output."""
# 最终先导排序：合并所有候选、过滤、按 MPO 排序并输出。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from drugpipe.hit_to_lead.analog_gen import AnalogGenerator
from drugpipe.hit_to_lead.clustering import DiversityClusterer
from drugpipe.hit_to_lead.mpo import MPOScorer
from drugpipe.hit_to_lead.scaffold import ScaffoldAnalyzer

logger = logging.getLogger(__name__)


class LeadRanker:
    """
    Orchestrate the full Hit-to-Lead stage:
      1. Scaffold analysis
      2. Diversity clustering
      3. Analog generation (CReM)
      4. MPO scoring
      5. Final ranking and output
    """
    # 编排完整 Hit-to-Lead 阶段：骨架分析 → 多样性聚类 → 类似物生成(CReM) → MPO 打分 → 最终排序输出。

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        self.cfg = cfg
        self.out_dir = out_dir
        h2l = cfg.get("hit_to_lead", {})
        self.max_hits = int(h2l.get("max_hits_input", 100))
        self.top_n = int(h2l.get("output", {}).get("top_n_leads", 50))

        self.scaffold_analyzer = ScaffoldAnalyzer(cfg)
        self.clusterer = DiversityClusterer(cfg)
        self.analog_gen = AnalogGenerator(cfg)
        self.mpo_scorer = MPOScorer(cfg)

        self.scaffold_csv = out_dir / "h2l_scaffold_summary.csv"
        self.cluster_csv = out_dir / "h2l_cluster_summary.csv"
        self.leads_csv = out_dir / "final_lead_candidates.csv"

    # ------------------------------------------------------------------
    def run(
        self,
        df_hits: pd.DataFrame,
        model_predict_fn=None,
        featurizer_fn=None,
    ) -> pd.DataFrame:
        """
        Execute the full H2L pipeline.

        Parameters
        ----------
        df_hits : Hit compounds from Stage 2 (must have ``canonical_smiles``,
                  ``pred_pIC50_ens``).
        model_predict_fn : callable(X_fp) → pIC50 array, for scoring new analogs.
        featurizer_fn : callable(smiles_list) → X_fp array.

        Returns
        -------
        DataFrame of ranked lead candidates.
        """
        # 返回排序后的先导候选 DataFrame。
        if len(df_hits) > self.max_hits:
            logger.info("Capping hits from %d to %d", len(df_hits), self.max_hits)
            df_hits = df_hits.head(self.max_hits).copy()

        # 1. Scaffold analysis
        logger.info("=== H2L Step 1: Scaffold Analysis ===")
        df_hits, df_scaffolds = self.scaffold_analyzer.analyze(df_hits)
        if not df_scaffolds.empty:
            df_scaffolds.to_csv(self.scaffold_csv, index=False)

        # 2. Diversity clustering
        logger.info("=== H2L Step 2: Diversity Clustering ===")
        df_hits, df_clusters = self.clusterer.cluster(df_hits)
        if not df_clusters.empty:
            df_clusters.to_csv(self.cluster_csv, index=False)

        # 3. Analog generation
        logger.info("=== H2L Step 3: Analog Generation ===")
        df_pool = self.analog_gen.generate(df_hits)
        logger.info("Candidate pool size: %d", len(df_pool))

        # 4. MPO scoring
        logger.info("=== H2L Step 4: MPO Scoring ===")
        df_scored = self.mpo_scorer.score(
            df_pool,
            model_predict_fn=model_predict_fn,
            featurizer_fn=featurizer_fn,
        )

        # 5. Rank and output
        logger.info("=== H2L Step 5: Lead Ranking ===")
        df_leads = self._rank_and_output(df_scored)
        return df_leads

    # ------------------------------------------------------------------
    def _rank_and_output(self, df: pd.DataFrame) -> pd.DataFrame:
        """Sort by MPO, take top-N, write CSV. / 按 MPO 排序、取 top-N、写入 CSV。"""
        df = df.dropna(subset=["mpo_score"])
        df = df.sort_values("mpo_score", ascending=False)
        df = df.head(self.top_n).reset_index(drop=True)

        keep_cols = [
            "canonical_smiles", "source_smiles", "origin",
            "pred_pIC50_ens", "QED", "mpo_score",
            "admet_pass", "HasAlert",
            "scaffold_smi", "cluster_id",
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]
        df_out = df[keep_cols].copy()

        df_out.to_csv(self.leads_csv, index=False)
        logger.info("Lead candidates saved: %s (%d rows)", self.leads_csv, len(df_out))
        return df_out
