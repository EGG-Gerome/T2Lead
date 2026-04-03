"""Discover approved / late-clinical drugs for a ChEMBL target via ChEMBL API."""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from drugpipe.utils.http import HTTPClient

logger = logging.getLogger(__name__)

_CHEMBL_COMPOUND_URL = "https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"


class ApprovedDrugFinder:
    """Fetch benchmark molecules: mechanism-of-action first, then activity fallback."""

    def __init__(self, cfg: Dict[str, Any]):
        bm = cfg.get("benchmark", {}) or {}
        t2h = cfg.get("target_to_hit", {}) or {}
        chembl = t2h.get("chembl", {}) or {}
        raw_base = bm.get("chembl_base_url") or chembl.get("base_url", "https://www.ebi.ac.uk/chembl/api/data")
        self.base_url = str(raw_base or "https://www.ebi.ac.uk/chembl/api/data").rstrip("/")
        self.max_drugs = int(bm.get("max_drugs", 5))
        self.min_phase = float(bm.get("min_phase", 3))
        self.http = HTTPClient(timeout=90, retries=4, polite_sleep=0.2)

    def find(self, target_chembl_id: str) -> List[Dict[str, Any]]:
        """Return up to ``max_drugs`` dicts with smiles, chembl_id, phase, pref_name, best_pchembl."""
        if not target_chembl_id or not str(target_chembl_id).strip():
            logger.warning("Benchmark: empty target_chembl_id — skipping drug lookup.")
            return []

        tid = str(target_chembl_id).strip()
        mol_ids = self._mechanism_molecule_ids(tid)
        if not mol_ids:
            logger.info(
                "Benchmark: no mechanism records for %s — trying activity fallback.",
                tid,
            )
            mol_ids = self._activity_molecule_ids(tid)

        if not mol_ids:
            logger.warning(
                "Benchmark: no ChEMBL drugs/activities found for target %s. "
                "The target may be undrugged in ChEMBL or use a different target id.",
                tid,
            )
            return []

        candidates: List[Dict[str, Any]] = []
        for mid in mol_ids:
            info = self._fetch_molecule(mid)
            if info is None:
                continue
            phase = info.get("max_phase")
            try:
                ph = float(phase) if phase is not None else 0.0
            except (TypeError, ValueError):
                ph = 0.0
            if ph < self.min_phase:
                continue
            smi = info.get("canonical_smiles") or ""
            if not smi:
                continue
            cid = info.get("molecule_chembl_id") or mid
            best_pchembl = self._best_pchembl_for_pair(tid, cid)
            candidates.append(
                {
                    "chembl_id": cid,
                    "canonical_smiles": smi,
                    "pref_name": info.get("pref_name") or "",
                    "max_phase": ph,
                    "first_approval": info.get("first_approval"),
                    "chembl_url": _CHEMBL_COMPOUND_URL.format(chembl_id=cid),
                    "best_pchembl": best_pchembl,
                    "mechanism_of_action": info.get("_moa") or "",
                }
            )

        if not candidates:
            logger.warning(
                "Benchmark: found %d ChEMBL molecule id(s) for %s but none with max_phase >= %.1f.",
                len(mol_ids),
                tid,
                self.min_phase,
            )
            return []

        candidates.sort(key=lambda x: (-x["max_phase"], -(x["best_pchembl"] or 0.0)))
        out = candidates[: self.max_drugs]
        logger.info(
            "Benchmark: selected %d reference drug(s) for target %s (min_phase >= %.1f).",
            len(out),
            tid,
            self.min_phase,
        )
        return out

    def _mechanism_molecule_ids(self, target_chembl_id: str) -> List[str]:
        seen_order: List[str] = []
        seen_set: set[str] = set()
        offset = 0
        limit = 500
        url = f"{self.base_url}/mechanism.json"
        while True:
            try:
                data = self.http.get_json(
                    url,
                    {"target_chembl_id": target_chembl_id, "limit": limit, "offset": offset},
                )
            except Exception as exc:
                logger.debug("Mechanism query failed at offset %d: %s", offset, exc)
                break
            rows = data.get("mechanisms") or []
            for row in rows:
                mid = row.get("molecule_chembl_id") or row.get("drug_chembl_id")
                if mid and mid not in seen_set:
                    seen_set.add(str(mid))
                    seen_order.append(str(mid))
            meta = data.get("page_meta", {})
            if not meta.get("next"):
                break
            offset += limit
            if offset > 5000:
                break
        return seen_order

    def _activity_molecule_ids(self, target_chembl_id: str) -> List[str]:
        """Collect molecules with reported pchembl on target, then caller filters by phase."""
        best_pchembl: Dict[str, float] = {}
        offset = 0
        limit = 1000
        url = f"{self.base_url}/activity.json"
        while True:
            try:
                data = self.http.get_json(
                    url,
                    {
                        "target_chembl_id": target_chembl_id,
                        "pchembl_value__isnull": False,
                        "limit": limit,
                        "offset": offset,
                    },
                )
            except Exception as exc:
                logger.debug("Activity query failed at offset %d: %s", offset, exc)
                break
            rows = data.get("activities") or []
            for row in rows:
                mid = row.get("molecule_chembl_id")
                pc = row.get("pchembl_value")
                if not mid or pc is None:
                    continue
                try:
                    v = float(pc)
                except (TypeError, ValueError):
                    continue
                k = str(mid)
                best_pchembl[k] = max(best_pchembl.get(k, 0.0), v)
            meta = data.get("page_meta", {})
            if not meta.get("next"):
                break
            offset += limit
            if offset > 10000:
                break
        ordered = sorted(best_pchembl.keys(), key=lambda m: -best_pchembl[m])
        return ordered[:200]

    def _fetch_molecule(self, molecule_chembl_id: str) -> Optional[Dict[str, Any]]:
        url = f"{self.base_url}/molecule.json"
        try:
            data = self.http.get_json(url, {"molecule_chembl_id": molecule_chembl_id, "limit": 1})
        except Exception:
            return None
        mols = data.get("molecules") or []
        if not mols:
            return None
        mol = mols[0]
        structs = mol.get("molecule_structures")
        smi = ""
        if isinstance(structs, list) and structs:
            smi = structs[0].get("canonical_smiles") or ""
        elif isinstance(structs, dict):
            smi = structs.get("canonical_smiles") or ""
        return {
            "molecule_chembl_id": mol.get("molecule_chembl_id") or molecule_chembl_id,
            "pref_name": mol.get("pref_name") or "",
            "max_phase": mol.get("max_phase"),
            "first_approval": mol.get("first_approval"),
            "canonical_smiles": smi,
            "_moa": "",
        }

    def _best_pchembl_for_pair(self, target_chembl_id: str, molecule_chembl_id: str) -> Optional[float]:
        url = f"{self.base_url}/activity.json"
        try:
            data = self.http.get_json(
                url,
                {
                    "target_chembl_id": target_chembl_id,
                    "molecule_chembl_id": molecule_chembl_id,
                    "pchembl_value__isnull": False,
                    "limit": 1,
                    "order_by": "-pchembl_value",
                },
            )
        except Exception:
            return None
        acts = data.get("activities") or []
        if not acts:
            return None
        pc = acts[0].get("pchembl_value")
        try:
            return float(pc) if pc is not None else None
        except (TypeError, ValueError):
            return None
