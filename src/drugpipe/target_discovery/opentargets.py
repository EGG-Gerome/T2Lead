"""Query the OpenTargets Platform GraphQL API for disease-target associations."""
# EN: Query the OpenTargets Platform GraphQL API for disease-target associations.
# 中文：说明模块职责、上下游关系与维护注意事项。

# 查询 OpenTargets Platform GraphQL API 获取疾病-靶点关联。

from __future__ import annotations

import logging
from typing import Any, Dict, List

from drugpipe.utils.http import HTTPClient

logger = logging.getLogger(__name__)

# GraphQL query: search for a disease, then retrieve associated targets.
# GraphQL 查询：按关键词搜索疾病，再获取关联靶点。
_DISEASE_SEARCH_QUERY = """
query DiseaseSearch($keyword: String!, $size: Int!) {
  search(queryString: $keyword, entityNames: ["disease"], page: {index: 0, size: $size}) {
    hits {
      id
      name
      entity
    }
  }
}
"""

_ASSOCIATIONS_QUERY = """
query DiseaseTargets($diseaseId: String!, $size: Int!) {
  disease(efoId: $diseaseId) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""

# Map Ensembl gene IDs to ChEMBL target IDs through the target detail query.
# 通过靶点详情查询将 Ensembl 基因 ID 映射为 ChEMBL 靶点 ID。
_TARGET_CHEMBL_QUERY = """
query TargetChembl($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    dbXrefs {
      id
      source
    }
  }
}
"""


class OpenTargetsClient:
    """Lightweight wrapper for OpenTargets Platform v4 GraphQL API."""
    # OpenTargets Platform v4 GraphQL API 的轻量封装。

    def __init__(self, cfg: Dict[str, Any]):
        ot_cfg = cfg.get("target_discovery", {}).get("opentargets", {})
        self.api_url = ot_cfg.get("api_url", "https://api.platform.opentargets.org/api/v4/graphql")
        self.min_score = float(ot_cfg.get("min_score", 0.1))
        self.http = HTTPClient(timeout=30, retries=4)

    # ------------------------------------------------------------------
    def search_disease(self, keyword: str, size: int = 10) -> List[Dict[str, Any]]:
        """Free-text search for disease entities, returning list of ``{id, name}``."""
        # 疾病实体自由文本搜索，返回 {id, name} 列表。
        payload = {
            "query": _DISEASE_SEARCH_QUERY,
            "variables": {"keyword": keyword, "size": size},
        }
        data = self.http.post_json(self.api_url, payload)
        hits = data.get("data", {}).get("search", {}).get("hits", [])
        return [{"id": h["id"], "name": h["name"]} for h in hits if h.get("entity") == "disease"]

    # ------------------------------------------------------------------
    def get_associated_targets(
        self,
        disease_id: str,
        size: int = 50,
    ) -> List[Dict[str, Any]]:
        """
        Return targets associated with *disease_id* (an EFO ID like ``EFO_0001378``).

        Each entry: ``{ensembl_id, symbol, name, score, datatype_scores}``.
        """
        # 返回与 disease_id（EFO ID）关联的靶点，每项含 ensembl_id、symbol、name、score、datatype_scores。
        payload = {
            "query": _ASSOCIATIONS_QUERY,
            "variables": {"diseaseId": disease_id, "size": size},
        }
        data = self.http.post_json(self.api_url, payload)
        disease_data = data.get("data", {}).get("disease")
        if disease_data is None:
            logger.warning("No disease found for ID %s", disease_id)
            return []

        rows = disease_data.get("associatedTargets", {}).get("rows", [])
        results = []
        for row in rows:
            score = row.get("score", 0)
            if score < self.min_score:
                continue
            tgt = row.get("target", {})
            results.append({
                "ensembl_id": tgt.get("id", ""),
                "symbol": tgt.get("approvedSymbol", ""),
                "name": tgt.get("approvedName", ""),
                "score": score,
                "datatype_scores": {
                    d["id"]: d["score"] for d in row.get("datatypeScores", [])
                },
            })
        return results

    # ------------------------------------------------------------------
    def ensembl_to_chembl(self, ensembl_id: str) -> List[str]:
        """Map an Ensembl gene ID to ChEMBL target IDs."""
        # 将 Ensembl 基因 ID 映射为 ChEMBL 靶点 ID 列表。
        payload = {
            "query": _TARGET_CHEMBL_QUERY,
            "variables": {"ensemblId": ensembl_id},
        }
        try:
            data = self.http.post_json(self.api_url, payload)
        except RuntimeError:
            return []
        tgt = data.get("data", {}).get("target")
        if tgt is None:
            return []

        # OpenTargets schema now exposes ChEMBL IDs via dbXrefs.
        xrefs = tgt.get("dbXrefs") or []
        chembl_ids = []
        for x in xrefs:
            xid = (x or {}).get("id")
            src = ((x or {}).get("source") or "").lower()
            if not xid:
                continue
            if src == "chembl" or str(xid).startswith("CHEMBL"):
                chembl_ids.append(str(xid))

        # Preserve order while removing duplicates.
        return list(dict.fromkeys(chembl_ids))

    # ------------------------------------------------------------------
    def discover(self, disease_query: str, top_n: int = 5) -> List[Dict[str, Any]]:
        """
        High-level helper: disease name → ranked targets with ChEMBL IDs.

        Returns a list sorted by association score, each entry including
        ``chembl_ids`` resolved from Ensembl.
        """
        # 高层接口：疾病名 → 带 ChEMBL ID 的排序靶点列表，按关联得分排序。
        diseases = self.search_disease(disease_query)
        if not diseases:
            logger.warning("No disease matches for '%s'", disease_query)
            return []

        disease_id = diseases[0]["id"]
        logger.info("Resolved disease '%s' → %s (%s)", disease_query, disease_id, diseases[0]["name"])

        targets = self.get_associated_targets(disease_id, size=max(top_n * 3, 50))

        for tgt in targets:
            tgt["chembl_ids"] = self.ensembl_to_chembl(tgt["ensembl_id"])

        # Keep only targets that have at least one ChEMBL ID
        targets = [t for t in targets if t["chembl_ids"]]
        targets.sort(key=lambda t: t["score"], reverse=True)

        return targets[:top_n]
