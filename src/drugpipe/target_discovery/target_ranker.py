"""Rank targets from OpenTargets and output ChEMBL IDs."""

from __future__ import annotations

import logging
from typing import Any, Dict, List

import pandas as pd

from drugpipe.target_discovery.opentargets import OpenTargetsClient

logger = logging.getLogger(__name__)


class TargetRanker:
    """
    Orchestrate Stage 1: query OpenTargets for disease-associated targets,
    rank, and emit the top-N ChEMBL target IDs.
    """

    def __init__(self, cfg: Dict[str, Any]):
        self.cfg = cfg
        td_cfg = cfg.get("target_discovery", {})
        self.disease = td_cfg.get("disease", "")
        self.top_n = int(td_cfg.get("top_n_targets", 5))
        self.use_opentargets = td_cfg.get("opentargets", {}).get("enabled", True)

    # ------------------------------------------------------------------
    def run(self) -> List[Dict[str, Any]]:
        """
        Execute target discovery and return ranked target list.

        Each entry: ``{chembl_id, symbol, name, score, source}``.
        """
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

        if not all_targets:
            logger.warning("No targets found for disease '%s'", self.disease)
            return []

        df = pd.DataFrame(all_targets)
        df = df.sort_values("score", ascending=False).drop_duplicates("chembl_id")
        df["rank_score"] = df["score"]
        df = df.sort_values("rank_score", ascending=False).head(self.top_n)

        ranked = df.to_dict(orient="records")
        for i, t in enumerate(ranked, 1):
            logger.info(
                "  Target #%d: %s (%s) score=%.3f [%s]",
                i, t["chembl_id"], t["symbol"], t["rank_score"], t["source"],
            )

        return ranked
