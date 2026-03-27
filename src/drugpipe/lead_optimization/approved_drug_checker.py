"""Check whether output molecules are known approved / clinical-stage drugs.

Queries the ChEMBL REST API (flexmatch on canonical SMILES) and annotates
a DataFrame with clinical-phase metadata so users can immediately tell
if an optimised lead is already a marketed drug.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional
from urllib.parse import quote

import numpy as np
import pandas as pd

from drugpipe.utils.http import HTTPClient

logger = logging.getLogger(__name__)

_CHEMBL_COMPOUND_URL = (
    "https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"
)


class ApprovedDrugChecker:
    """Annotate molecules with ChEMBL approval / clinical-phase metadata."""

    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        adc = lo.get("approved_drug_check", {})
        self.enabled = bool(adc.get("enabled", True))
        self.base_url = adc.get(
            "chembl_base_url",
            "https://www.ebi.ac.uk/chembl/api/data",
        )
        self.http = HTTPClient(timeout=60, retries=4, polite_sleep=0.3)

    # ------------------------------------------------------------------
    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add approval columns to *df* (in-place) and return it.

        New columns: is_approved, max_phase, pref_name,
                     chembl_id, chembl_url, first_approval.
        """
        if not self.enabled:
            logger.info("Approved-drug check disabled, skipping.")
            return df

        n = len(df)
        logger.info(
            "=== Lead Optimization Step 6: Approved Drug Check (%d molecules) ===", n
        )

        is_approved = []
        max_phase = []
        pref_name = []
        chembl_id = []
        chembl_url = []
        first_approval = []

        for idx, smi in enumerate(df["canonical_smiles"]):
            info = self._query_chembl(smi)
            if info is not None:
                raw = info.get("max_phase")
                try:
                    phase = float(raw) if raw is not None else 0.0
                except (TypeError, ValueError):
                    phase = 0.0
                is_approved.append(phase >= 4.0)
                max_phase.append(phase)
                pref_name.append(info.get("pref_name") or "")
                cid = info.get("molecule_chembl_id") or ""
                chembl_id.append(cid)
                chembl_url.append(
                    _CHEMBL_COMPOUND_URL.format(chembl_id=cid) if cid else ""
                )
                yr = info.get("first_approval")
                first_approval.append(int(yr) if yr else np.nan)
            else:
                is_approved.append(False)
                max_phase.append(np.nan)
                pref_name.append("")
                chembl_id.append("")
                chembl_url.append("")
                first_approval.append(np.nan)

            logger.debug(
                "  [%d/%d] %s → %s",
                idx + 1, n, smi[:60],
                f"phase={max_phase[-1]}, name={pref_name[-1]}" if pref_name[-1] else "not found",
            )

        df["is_approved"] = is_approved
        df["max_phase"] = max_phase
        df["pref_name"] = pref_name
        df["chembl_id"] = chembl_id
        df["chembl_url"] = chembl_url
        df["first_approval"] = first_approval

        n_approved = sum(is_approved)
        logger.info(
            "Approved drug check: %d/%d molecules matched approved drugs (max_phase=4)",
            n_approved, n,
        )
        return df

    # ------------------------------------------------------------------
    def _query_chembl(self, smiles: str) -> Optional[Dict[str, Any]]:
        """Query ChEMBL molecule endpoint by SMILES flexmatch.

        Returns the first matching molecule dict, or *None* on miss / error.
        """
        url = f"{self.base_url}/molecule.json"
        params = {
            "molecule_structures__canonical_smiles__flexmatch": smiles,
            "limit": 1,
        }
        try:
            data = self.http.get_json(url, params)
        except Exception:
            logger.warning("ChEMBL query failed for SMILES: %s", smiles[:80], exc_info=True)
            return None

        molecules = data.get("molecules") or []
        if not molecules:
            return None
        return molecules[0]
