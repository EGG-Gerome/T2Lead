"""Enhanced ADMET profiling beyond basic Lipinski / QED.

All predictions use RDKit — no extra dependencies required.

Computed properties:
  - Synthetic Accessibility (SA) score  (1=easy … 10=hard)
  - hERG cardiotoxicity liability       (SMARTS toxicophores)
  - CYP3A4 / CYP2D6 inhibition risk    (SMARTS pharmacophores)
  - Veber oral bioavailability rules    (RotBonds ≤ 10, TPSA ≤ 140)
  - Composite ADMET risk score          (0 = low risk … 1 = high risk)
"""
# EN: Enhanced ADMET profiling beyond basic Lipinski / QED.
# 中文：说明模块职责、上下游关系与维护注意事项。

# 增强 ADMET 评估：合成可及性、hERG 毒性、CYP 抑制风险、Veber 规则、综合风险评分。

from __future__ import annotations

import logging
from typing import Any, Dict, List

import importlib.util
import os
import sys

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, RDConfig

from drugpipe.utils.chem import safe_mol, calc_descriptors

# SA Score lives in RDKit's Contrib directory and isn't directly importable.
_SA_SCORER = None
_sa_scorer_path = os.path.join(RDConfig.RDContribDir, "SA_Score", "sascorer.py")
if os.path.exists(_sa_scorer_path):
    spec = importlib.util.spec_from_file_location("sascorer", _sa_scorer_path)
    _SA_SCORER = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(_SA_SCORER)
    except Exception:
        _SA_SCORER = None

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------
# SMARTS patterns for hERG liability (literature-derived toxicophores)
# hERG 毒性 SMARTS（文献衍生毒性基团）
# -----------------------------------------------------------------------
_HERG_SMARTS = [
    Chem.MolFromSmarts("[NX3](C)(C)CCc1ccccc1"),       # basic nitrogen + phenylethyl
    Chem.MolFromSmarts("[NX3]CCOc1ccccc1"),             # aminoethoxy-phenyl
    Chem.MolFromSmarts("c1ccc2c(c1)CCN(C2)"),          # tetrahydroisoquinoline
    Chem.MolFromSmarts("[NX3](C)(C)CCc1ccc(F)cc1"),    # fluoro-phenylethylamine
    Chem.MolFromSmarts("c1ccc(cc1)C(O)CCN"),           # phenyl-propanolamine
]
_HERG_SMARTS = [p for p in _HERG_SMARTS if p is not None]

# -----------------------------------------------------------------------
# SMARTS patterns for CYP3A4 / CYP2D6 inhibition risk
# CYP 抑制 SMARTS（常见药效团）
# -----------------------------------------------------------------------
_CYP_SMARTS = [
    Chem.MolFromSmarts("c1cc[nH]c1"),                  # imidazole (CYP3A4)
    Chem.MolFromSmarts("c1ccncc1"),                    # pyridine  (CYP2D6)
    Chem.MolFromSmarts("[NX3](C)(C)c1ccccc1"),         # N,N-dialkylaniline
    Chem.MolFromSmarts("c1ccc2[nH]ccc2c1"),            # indole
    Chem.MolFromSmarts("c1cnc2ccccc2n1"),              # quinazoline
]
_CYP_SMARTS = [p for p in _CYP_SMARTS if p is not None]


class DeepADMET:
    """Enhanced ADMET profiling for lead candidates."""

    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        ad = lo.get("admet_deep", {})
        self.enabled = bool(ad.get("enabled", True))
        self.sa_max = float(ad.get("sa_score_max", 6.0))
        self.herg_filter = bool(ad.get("herg_filter", True))

    # ------------------------------------------------------------------
    def profile(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add ADMET columns and a composite ``admet_risk`` score."""
        df = df.copy()

        if not self.enabled:
            logger.info("Enhanced ADMET profiling disabled.")
            df["admet_risk"] = 0.0
            return df

        df["sa_score"] = np.nan
        df["herg_flag"] = False
        df["cyp_risk"] = False
        df["veber_pass"] = True
        df["admet_risk"] = np.nan

        for idx in df.index:
            smi = df.at[idx, "canonical_smiles"]
            mol = safe_mol(smi)
            if mol is None:
                df.at[idx, "admet_risk"] = 1.0
                continue

            sa = self._sa_score(mol)
            herg = self._herg_liability(mol) if self.herg_filter else False
            cyp = self._cyp_inhibition(mol)
            veber = self._veber_rules(mol)

            df.at[idx, "sa_score"] = sa
            df.at[idx, "herg_flag"] = herg
            df.at[idx, "cyp_risk"] = cyp
            df.at[idx, "veber_pass"] = veber
            df.at[idx, "admet_risk"] = self._composite_risk(sa, herg, cyp, veber)

        n_risky = (df["admet_risk"] > 0.5).sum()
        logger.info(
            "Enhanced ADMET complete: %d / %d flagged as high-risk (>0.5).",
            n_risky, len(df),
        )
        return df

    # ------------------------------------------------------------------
    @staticmethod
    def _sa_score(mol: Chem.Mol) -> float:
        if _SA_SCORER is not None:
            try:
                return float(_SA_SCORER.calculateScore(mol))
            except Exception:
                pass
        # Fallback: rough SA estimate from heavy-atom count + ring complexity
        ha = mol.GetNumHeavyAtoms()
        rings = Lipinski.RingCount(mol)
        return min(1.0 + ha * 0.05 + rings * 0.3, 10.0)

    @staticmethod
    def _herg_liability(mol: Chem.Mol) -> bool:
        return any(mol.HasSubstructMatch(pat) for pat in _HERG_SMARTS)

    @staticmethod
    def _cyp_inhibition(mol: Chem.Mol) -> bool:
        matches = sum(1 for pat in _CYP_SMARTS if mol.HasSubstructMatch(pat))
        return matches >= 2  # flag when ≥2 pharmacophores present

    @staticmethod
    def _veber_rules(mol: Chem.Mol) -> bool:
        """Veber rules for oral bioavailability."""
        rotb = Lipinski.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        return rotb <= 10 and tpsa <= 140

    def _composite_risk(
        self,
        sa: float,
        herg: bool,
        cyp: bool,
        veber: bool,
    ) -> float:
        """0 = low risk, 1 = high risk.

        Weighted combination:
          SA contribution:   normalised SA/10, weight 0.30
          hERG:              1 if flagged, weight 0.30
          CYP:               1 if flagged, weight 0.15
          Veber violation:   1 if fails,  weight 0.25
        """
        risk = 0.0
        risk += 0.30 * min(sa / 10.0, 1.0)
        risk += 0.30 * float(herg)
        risk += 0.15 * float(cyp)
        risk += 0.25 * float(not veber)
        return round(risk, 4)
