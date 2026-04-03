"""ChEMBL REST API crawler with checkpoint-based resume and rate limiting."""
# ChEMBL REST API crawler with checkpoint-based resume and rate limiting.
# 说明模块职责、上下游关系与维护注意事项。

# ChEMBL REST API 爬虫，支持断点续跑与请求限速。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from drugpipe.utils.http import HTTPClient
from drugpipe.utils.io import append_csv, load_state, save_state

logger = logging.getLogger(__name__)


def resolve_chembl_target_by_gene(
    gene_symbol: str,
    base_url: str = "https://www.ebi.ac.uk/chembl/api/data",
) -> Optional[str]:
    """Look up a single-protein ChEMBL target ID by gene symbol.
    通过基因符号查询 ChEMBL 单蛋白靶点 ID（人源优先）。

    Returns the best-matching ``CHEMBLXXX`` target ID, or *None* if no
    human single-protein target is found for *gene_symbol*.
    """
    from drugpipe.utils.http import HTTPClient

    http = HTTPClient(timeout=30, retries=4, polite_sleep=0.1)
    url = f"{base_url}/target/search.json"
    try:
        data = http.get_json(url, {"q": gene_symbol, "limit": 20})
    except Exception as exc:
        logger.warning("ChEMBL target search failed for %s: %s", gene_symbol, exc)
        return None

    targets = data.get("targets") or data.get("target") or []
    gene_upper = gene_symbol.upper()

    for t in targets:
        if t.get("target_type") != "SINGLE PROTEIN":
            continue
        org = (t.get("organism") or "").lower()
        if "homo sapiens" not in org:
            continue
        components = t.get("target_components") or []
        for comp in components:
            synonyms = comp.get("target_component_synonyms") or []
            for syn in synonyms:
                if (syn.get("component_synonym") or "").upper() == gene_upper:
                    tid = t.get("target_chembl_id")
                    logger.info(
                        "Gene %s → ChEMBL target %s (%s)",
                        gene_symbol, tid, t.get("pref_name", ""),
                    )
                    return tid

    if targets:
        for t in targets:
            if t.get("target_type") != "SINGLE PROTEIN":
                continue
            pref = (t.get("pref_name") or "").upper()
            if gene_upper in pref:
                tid = t.get("target_chembl_id")
                logger.info(
                    "Gene %s → ChEMBL target %s (name match: %s)",
                    gene_symbol, tid, t.get("pref_name", ""),
                )
                return tid

    logger.warning("No ChEMBL SINGLE PROTEIN target found for gene %s", gene_symbol)
    return None


class ChEMBLCrawler:
    """
    Download molecules (SMILES) and IC50 activity data from ChEMBL.

    Supports:
      - Configurable page size and maximum record counts.
      - Checkpoint/resume via a JSON state file.
      - Polite sleep between requests.
    """
    # 从 ChEMBL 下载分子（SMILES）与 IC50 活性数据；支持可配置分页与上限、JSON 断点续跑、请求间延时。

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        ch = cfg.get("target_to_hit", {}).get("chembl", {})
        self.base_url = ch.get("base_url", "https://www.ebi.ac.uk/chembl/api/data")
        self.page_limit = int(ch.get("page_limit", 1000))
        self.max_molecules = int(ch.get("max_molecules", 0))
        self.max_activities = int(ch.get("max_activities", 0))
        self.require_nM = bool(ch.get("require_nM", True))

        self.http = HTTPClient(
            timeout=60,
            retries=6,
            polite_sleep=float(ch.get("polite_sleep", 0.05)),
        )

        self.out_dir = out_dir
        self.molecules_csv = out_dir / "molecules_chemblid_smiles.csv"
        self.activities_csv = out_dir / "activities_ic50.csv"
        self.state_path = out_dir / "crawl_state.json"

    # ==================================================================
    # Molecules
    # 分子爬取
    # ==================================================================
    def crawl_molecules(self) -> Path:
        """Crawl molecule records (chembl_id + canonical_smiles)."""
        # 爬取分子记录（chembl_id + canonical_smiles）。
        state = load_state(self.state_path)
        if self.molecules_csv.exists() and state.get("molecule_done"):
            logger.info("Molecule crawl marked done in state, skipping.")
            return self.molecules_csv
        if self.molecules_csv.exists() and not state.get("molecule_done"):
            logger.info(
                "Resuming molecule crawl from existing CSV (offset=%d).",
                int(state.get("molecule_offset", 0)),
            )

        url = f"{self.base_url}/molecule.json"
        offset = int(state.get("molecule_offset", 0))
        total = 0

        logger.info("Starting molecule crawl (offset=%d, limit=%d)", offset, self.page_limit)

        while True:
            if self.max_molecules and total >= self.max_molecules:
                logger.info("Reached max_molecules=%d, stopping.", self.max_molecules)
                break

            data = self.http.get_json(url, {"limit": self.page_limit, "offset": offset})
            items = data.get("molecules") or data.get("molecule") or []
            if not items:
                state["molecule_done"] = True
                break

            rows = self._extract_molecules(items)
            if rows:
                append_csv(self.molecules_csv, pd.DataFrame(rows).drop_duplicates())
                total += len(rows)

            offset += self.page_limit
            state["molecule_offset"] = offset
            save_state(state, self.state_path)
            logger.info("  molecules offset=%d total_written=%d", offset, total)

        save_state(state, self.state_path)
        return self.molecules_csv

    # ==================================================================
    # IC50 Activities
    # IC50 活性爬取
    # ==================================================================
    def crawl_activities(self) -> Path:
        """Crawl IC50 activity records."""
        # 爬取 IC50 活性记录。
        state = load_state(self.state_path)
        if self.activities_csv.exists() and state.get("activity_done"):
            logger.info("Activity crawl marked done in state, skipping.")
            return self.activities_csv
        if self.activities_csv.exists() and not state.get("activity_done"):
            logger.info(
                "Resuming activity crawl from existing CSV (offset=%d).",
                int(state.get("activity_offset", 0)),
            )

        url = f"{self.base_url}/activity.json"
        offset = int(state.get("activity_offset", 0))
        total = 0

        logger.info("Starting activity crawl (offset=%d, require_nM=%s)", offset, self.require_nM)

        while True:
            if self.max_activities and total >= self.max_activities:
                logger.info("Reached max_activities=%d, stopping.", self.max_activities)
                break

            params: Dict[str, Any] = {
                "limit": self.page_limit,
                "offset": offset,
                "standard_type": "IC50",
            }
            if self.require_nM:
                params["standard_units"] = "nM"

            data = self.http.get_json(url, params)
            items = data.get("activities") or data.get("activity") or []
            if not items:
                state["activity_done"] = True
                break

            rows = self._extract_activities(items)
            if rows:
                df = pd.DataFrame(rows).drop_duplicates(subset=["activity_id"])
                append_csv(self.activities_csv, df)
                total += len(df)

            offset += self.page_limit
            state["activity_offset"] = offset
            save_state(state, self.state_path)
            logger.info("  activities offset=%d total_written=%d", offset, total)

        save_state(state, self.state_path)
        return self.activities_csv

    # ==================================================================
    # Internals
    # 内部方法
    # ==================================================================
    @staticmethod
    def _extract_molecules(items: List[Dict]) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        for it in items:
            cid = it.get("molecule_chembl_id")
            structs = it.get("molecule_structures") or {}
            smi = structs.get("canonical_smiles")
            if cid and smi:
                rows.append({"molecule_chembl_id": cid, "canonical_smiles": smi})
        return rows

    def _extract_activities(self, items: List[Dict]) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        for it in items:
            std_val = it.get("standard_value")
            std_rel = (it.get("standard_relation") or "").strip()
            std_units = (it.get("standard_units") or "").strip()

            try:
                v = float(std_val) if std_val is not None else None
            except (TypeError, ValueError):
                v = None
            if v is None or v <= 0:
                continue
            if self.require_nM and std_units != "nM":
                continue
            if std_rel and std_rel not in ("=", "<", "<=", ">", ">="):
                continue

            rows.append({
                "activity_id": it.get("activity_id"),
                "molecule_chembl_id": it.get("molecule_chembl_id"),
                "target_chembl_id": it.get("target_chembl_id"),
                "assay_chembl_id": it.get("assay_chembl_id"),
                "standard_type": it.get("standard_type"),
                "standard_relation": std_rel,
                "standard_value": v,
                "standard_units": std_units,
                "pchembl_value": it.get("pchembl_value"),
                "confidence_score": it.get("confidence_score"),
                "data_validity_comment": it.get("data_validity_comment"),
            })
        return rows
