"""Benchmark approved drugs against the same scoring stack as optimized leads."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import pandas as pd

from drugpipe.benchmark.drug_finder import ApprovedDrugFinder
from drugpipe.benchmark.scorer import score_benchmark_molecules

logger = logging.getLogger(__name__)


def run_benchmark(
    cfg: Dict[str, Any],
    target_chembl_id: Optional[str],
    protein_info: Dict[str, Any],
    out_dir: Path,
    model_predict_fn: Optional[Callable[..., Any]] = None,
    featurizer_fn: Optional[Callable[..., Any]] = None,
) -> Optional[pd.DataFrame]:
    """Find reference drugs for *target_chembl_id* and write ``benchmark_drugs.csv``."""
    bm = cfg.get("benchmark", {}) or {}
    if not bm.get("enabled", True):
        logger.info("Benchmark module disabled in config.")
        return None
    if not target_chembl_id:
        logger.info("Benchmark: no target_chembl_id — skipping.")
        return None

    try:
        finder = ApprovedDrugFinder(cfg)
        drugs: List[Dict[str, Any]] = finder.find(str(target_chembl_id))
    except Exception as exc:
        logger.warning("Benchmark drug lookup failed: %s", exc, exc_info=True)
        return None

    if not drugs:
        return None

    try:
        return score_benchmark_molecules(
            cfg,
            drugs,
            protein_info,
            out_dir,
            model_predict_fn=model_predict_fn,
            featurizer_fn=featurizer_fn,
        )
    except Exception as exc:
        logger.warning("Benchmark scoring failed: %s", exc, exc_info=True)
        return None
