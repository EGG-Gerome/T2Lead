"""Combine and rank targets from OpenTargets + OriGene, output ChEMBL IDs."""
# EN: Combine and rank targets from OpenTargets + OriGene, output ChEMBL IDs.
# 中文：说明模块职责、上下游关系与维护注意事项。

# 合并并排序 OpenTargets 与 OriGene 的靶点，输出 ChEMBL ID。

from __future__ import annotations

import logging
from typing import Any, Dict, List

import pandas as pd

from drugpipe.target_discovery.opentargets import OpenTargetsClient
from drugpipe.target_discovery.origene_client import OriGeneClient

logger = logging.getLogger(__name__)


class TargetRanker:
    """
    Orchestrate Stage 1: merge data-driven (OpenTargets) and AI-driven
    (OriGene) target lists, rank, and emit the top-N ChEMBL target IDs.
    """
    # 编排阶段一：合并数据驱动（OpenTargets）与 AI 驱动（OriGene）靶点列表，排序并输出 top-N ChEMBL 靶点 ID。

    def __init__(self, cfg: Dict[str, Any]):
        self.cfg = cfg
        td_cfg = cfg.get("target_discovery", {})
        self.disease = td_cfg.get("disease", "")
        self.top_n = int(td_cfg.get("top_n_targets", 5))
        self.use_origene = td_cfg.get("origene", {}).get("enabled", False)
        self.use_opentargets = td_cfg.get("opentargets", {}).get("enabled", True)

    # ------------------------------------------------------------------
    def run(self) -> List[Dict[str, Any]]:
        """
        Execute target discovery and return ranked target list.

        Each entry: ``{chembl_id, symbol, name, score, source}``.
        """
        # 执行靶点发现并返回排序后的靶点列表，每项含 chembl_id、symbol、name、score、source。
        if not self.disease:
            raise ValueError(
                "target_discovery.disease is empty. "
                "Set it in the config or pass TARGET_CHEMBL_ID directly to Stage 2."
            )

        all_targets: List[Dict[str, Any]] = []

        if self.use_opentargets:
            logger.info("Querying OpenTargets for '%s' ...", self.disease)
            ot = OpenTargetsClient(self.cfg)
            ot_targets = ot.discover(self.disease, top_n=self.top_n * 3)
            for t in ot_targets:
                for cid in t.get("chembl_ids", []):
                    all_targets.append({
                        "chembl_id": cid,
                        "symbol": t.get("symbol", ""),
                        "name": t.get("name", ""),
                        "score": t.get("score", 0.0),
                        "source": "opentargets",
                    })

        if self.use_origene:
            logger.info("Querying OriGene for '%s' ...", self.disease)
            og = OriGeneClient(self.cfg)
            og_targets = og.discover(self.disease)
            for t in og_targets:
                for cid in t.get("chembl_ids", []):
                    all_targets.append({
                        "chembl_id": cid,
                        "symbol": t.get("symbol", ""),
                        "name": t.get("name", ""),
                        "score": t.get("score", 0.0),
                        "source": "origene",
                    })

        if not all_targets:
            logger.warning("No targets found for disease '%s'", self.disease)
            return []

        df = pd.DataFrame(all_targets)

        # Deduplicate on chembl_id, keeping the highest score
        # 按 chembl_id 去重，保留最高分
        df = df.sort_values("score", ascending=False).drop_duplicates("chembl_id")

        # Boost targets that appear in both sources
        # 对同时出现在两个来源的靶点加分
        source_counts = (
            pd.DataFrame(all_targets)
            .groupby("chembl_id")["source"]
            .nunique()
            .rename("n_sources")
        )
        df = df.join(source_counts, on="chembl_id")
        df["rank_score"] = df["score"] + 0.1 * (df["n_sources"] - 1)

        df = df.sort_values("rank_score", ascending=False).head(self.top_n)

        ranked = df.to_dict(orient="records")
        for i, t in enumerate(ranked, 1):
            logger.info(
                "  Target #%d: %s (%s) score=%.3f [%s]",
                i, t["chembl_id"], t["symbol"], t["rank_score"], t["source"],
            )

        return ranked
