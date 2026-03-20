"""Wrapper for calling a locally deployed OriGene agent.

OriGene (https://github.com/GENTEL-lab/OriGene) is a multi-agent system that
functions as a virtual disease biologist.  This module provides a thin HTTP
interface to a running OriGene instance (which in turn relies on OrigeneMCP).

If the user has not deployed OriGene, the rest of the pipeline still works —
this stage is optional.
"""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 调用本地部署的 OriGene 代理的封装。OriGene 为多智能体虚拟疾病生物学家，本模块提供对其 HTTP 接口的薄封装；未部署时流水线其余阶段仍可运行。

from __future__ import annotations

import logging
import re
from typing import Any, Dict, List

from drugpipe.utils.http import HTTPClient

logger = logging.getLogger(__name__)


# EN: OriGeneClient core behavior and intent.
# 中文：OriGeneClient 的核心行为与设计意图。
class OriGeneClient:
    """Send a target-discovery query to a running OriGene service."""
    # 向运行中的 OriGene 服务发送靶点发现查询。

    # EN: __init__ core behavior and intent.
    # 中文：__init__ 的核心行为与设计意图。
    def __init__(self, cfg: Dict[str, Any]):
        og_cfg = cfg.get("target_discovery", {}).get("origene", {})
        self.server_url = og_cfg.get("server_url", "http://127.0.0.1:8788")
        self.mode = og_cfg.get("mode", "quick")
        self.http = HTTPClient(timeout=300, retries=2, backoff_base=5.0)

    # ------------------------------------------------------------------
    # EN: query core behavior and intent.
    # 中文：query 的核心行为与设计意图。
    def query(self, disease: str) -> Dict[str, Any]:
        """
        Ask OriGene for therapeutic target recommendations.

        Returns the raw response dict from the OriGene API (typically contains
        a ``report`` text field and potentially structured ``targets``).
        """
        # 向 OriGene 请求治疗靶点推荐，返回 API 原始响应（通常含 report 文本及可能的 targets 结构）。
        prompt = (
            f"What are the most promising therapeutic targets for {disease}? "
            "Please provide ChEMBL target IDs and gene symbols where possible."
        )
        url = f"{self.server_url.rstrip('/')}/api/research"
        payload = {"query": prompt, "mode": self.mode}
        try:
            result = self.http.post_json(url, payload)
            return result
        except RuntimeError as exc:
            logger.warning("OriGene query failed: %s", exc)
            return {}

    # ------------------------------------------------------------------
    # EN: discover core behavior and intent.
    # 中文：discover 的核心行为与设计意图。
    def discover(self, disease: str) -> List[Dict[str, Any]]:
        """
        High-level: query OriGene and extract target information.

        Attempts to parse ChEMBL IDs (``CHEMBL\\d+``) and gene symbols from
        the report text.  Returns a list of dicts with ``chembl_ids``,
        ``symbol``, and ``source``.
        """
        # 高层接口：查询 OriGene 并从报告中解析 ChEMBL ID 与基因符号，返回含 chembl_ids、symbol、source 的列表。
        raw = self.query(disease)
        if not raw:
            return []

        report = raw.get("report", "") or raw.get("result", "") or str(raw)
        return self._parse_targets(report)

    # ------------------------------------------------------------------
    @staticmethod
    # EN: _parse_targets core behavior and intent.
    # 中文：_parse_targets 的核心行为与设计意图。
    def _parse_targets(text: str) -> List[Dict[str, Any]]:
        """Best-effort extraction of target info from free-text report."""
        # 从自由文本报告中尽力提取靶点信息。
        chembl_ids = list(dict.fromkeys(re.findall(r"CHEMBL\d+", text)))

        gene_pattern = re.compile(r"\b([A-Z][A-Z0-9]{1,10})\b")
        # Common false-positive words to exclude / 常见误匹配词排除
        _stopwords = {
            "THE", "AND", "FOR", "WITH", "FROM", "THAT", "THIS", "ARE",
            "NOT", "BUT", "HAS", "WAS", "BEEN", "HAVE", "WILL", "CAN",
            "MAY", "ITS", "ALL", "ONE", "TWO", "NEW", "USE", "DNA", "RNA",
            "MCP", "API", "LLM", "SMILES", "ADMET", "IC50",
        }
        raw_symbols = list(dict.fromkeys(gene_pattern.findall(text)))
        symbols = [s for s in raw_symbols if s not in _stopwords and len(s) >= 2]

        targets: List[Dict[str, Any]] = []
        for cid in chembl_ids:
            targets.append({
                "chembl_ids": [cid],
                "symbol": "",
                "name": "",
                "score": 0.0,
                "source": "origene",
            })

        for sym in symbols[:20]:
            if not any(sym in str(t.get("chembl_ids", [])) for t in targets):
                targets.append({
                    "chembl_ids": [],
                    "symbol": sym,
                    "name": "",
                    "score": 0.0,
                    "source": "origene",
                })

        return targets
