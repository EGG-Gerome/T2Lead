"""Stage 4 orchestrator: protein prep → docking → ADMET → MD → ranking.

Coordinates all lead-optimization substeps and produces
``optimized_leads.csv`` with composite scores.
"""
# Stage 4 orchestrator: protein prep → docking → ADMET → MD → ranking.
# 说明模块职责、上下游关系与维护注意事项。

# 阶段四编排器：蛋白准备 → 对接 → ADMET → MD → 综合排序 → optimized_leads.csv。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from drugpipe.lead_optimization.admet_deep import DeepADMET
from drugpipe.lead_optimization.approved_drug_checker import ApprovedDrugChecker
from drugpipe.lead_optimization.docking import VinaDocking
from drugpipe.lead_optimization.md_simulation import MDSimulator
from drugpipe.lead_optimization.protein_prep import ProteinPreparator

logger = logging.getLogger(__name__)


class LeadOptimizer:
    """Orchestrate the full Lead Optimization stage.

    Substeps
    --------
    1. Protein preparation  (fetch PDB, fix, detect binding site)
    2. Molecular docking    (AutoDock Vina)
    3. Enhanced ADMET       (SA score, hERG, CYP, Veber)
    4. MD simulation        (OpenMM MM-GBSA, top-N only)
    5. Composite scoring & ranking
    """

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        self.cfg = cfg
        self.out_dir = out_dir
        lo = cfg.get("lead_optimization", {})
        self.enabled = bool(lo.get("enabled", True))
        self.top_n = int(lo.get("output", {}).get("top_n_optimized", 10))

        sc = lo.get("scoring", {})
        self.w_docking = float(sc.get("w_docking", 0.35))
        self.w_admet = float(sc.get("w_admet", 0.20))
        self.w_md = float(sc.get("w_md_energy", 0.30))
        self.w_stability = float(sc.get("w_stability", 0.15))

        self.protein_prep = ProteinPreparator(cfg)
        self.docker = VinaDocking(cfg)
        self.admet = DeepADMET(cfg)
        self.md_sim = MDSimulator(cfg)
        self.drug_checker = ApprovedDrugChecker(cfg)

        self.output_csv = out_dir / "optimized_leads.csv"

    # ------------------------------------------------------------------
    def run(self, df_leads: pd.DataFrame) -> pd.DataFrame:
        """Execute the full lead-optimization pipeline."""
        if not self.enabled:
            logger.info("Lead optimization stage disabled.")
            return df_leads

        logger.info("=== Lead Optimization Step 1: Protein Preparation ===")
        protein_info = self.protein_prep.prepare(self.out_dir)
        if protein_info is None:
            logger.warning(
                "Protein preparation skipped (no pdb_id). "
                "Docking and MD will be disabled."
            )
            protein_info = {}

        logger.info("=== Lead Optimization Step 2: Molecular Docking ===")
        if protein_info:
            df_leads = self.docker.dock_leads(df_leads, protein_info, self.out_dir)
        else:
            logger.info("Skipping docking — no protein structure.")
            df_leads = df_leads.copy()
            df_leads["docking_score"] = np.nan

        logger.info("=== Lead Optimization Step 3: Enhanced ADMET ===")
        df_leads = self.admet.profile(df_leads)

        logger.info("=== Lead Optimization Step 4: MD Simulation ===")
        if protein_info:
            df_leads = self.md_sim.run(df_leads, protein_info, self.out_dir)
        else:
            logger.info("Skipping MD — no protein structure.")
            df_leads["md_binding_energy"] = np.nan
            df_leads["md_rmsd_mean"] = np.nan

        logger.info("=== Lead Optimization Step 5: Composite Scoring ===")
        df_leads = self._composite_score(df_leads)
        df_out = self._rank_and_output(df_leads)
        return df_out

    # ------------------------------------------------------------------
    def _composite_score(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute a weighted composite optimization score.

        Higher is better.  Individual components are normalised to [0,1].
        When MD data is missing, its weight is redistributed proportionally
        to the components that have real data, so NaN components don't
        inflate or deflate scores with arbitrary defaults.
        """
        df = df.copy()

        dock_s = df.get("docking_score")
        dock_available = dock_s is not None and dock_s.notna().any()
        dock_norm = self._normalise_lower_better(dock_s) if dock_available else None

        admet_norm = 1.0 - df.get("admet_risk", pd.Series(0.5, index=df.index)).fillna(0.5)

        md_s = df.get("md_binding_energy")
        md_available = md_s is not None and md_s.notna().any()
        md_norm = self._normalise_lower_better(md_s) if md_available else None

        rmsd_s = df.get("md_rmsd_mean")
        stab_available = rmsd_s is not None and rmsd_s.notna().any()
        stab_norm = self._normalise_lower_better(rmsd_s) if stab_available else None

        components: list[tuple[float, pd.Series]] = []
        components.append((self.w_admet, admet_norm))
        if dock_norm is not None:
            components.append((self.w_docking, dock_norm))
        else:
            logger.warning("Docking data missing — excluded from composite score.")
        if md_norm is not None:
            components.append((self.w_md, md_norm))
        else:
            logger.warning("MD binding energy missing — excluded from composite score.")
        if stab_norm is not None:
            components.append((self.w_stability, stab_norm))
        else:
            logger.warning("MD stability (RMSD) missing — excluded from composite score.")

        w_total = sum(w for w, _ in components)
        score = pd.Series(0.0, index=df.index)
        for w, s in components:
            score += w * s
        df["opt_score"] = score / max(w_total, 1e-9)

        return df

    @staticmethod
    def _normalise_lower_better(s: Optional[pd.Series]) -> pd.Series:
        """Min-max normalise where *lower* raw values → *higher* normalised
        scores (closer to 1.0)."""
        if s is None or s.isna().all():
            return pd.Series(0.5, index=s.index if s is not None else pd.RangeIndex(0))
        filled = s.fillna(s.max() if s.max() == s.max() else 0)
        mn, mx = filled.min(), filled.max()
        if mx - mn < 1e-9:
            return pd.Series(0.5, index=s.index)
        return 1.0 - (filled - mn) / (mx - mn)

    # ------------------------------------------------------------------
    def _rank_and_output(self, df: pd.DataFrame) -> pd.DataFrame:
        """Sort by composite score, keep top-N, annotate approval status, write CSV."""
        df = df.sort_values("opt_score", ascending=False)
        df = df.head(self.top_n).reset_index(drop=True)

        df = self.drug_checker.annotate(df)

        keep_cols = [
            "canonical_smiles", "origin",
            "pred_pIC50_ens", "QED", "mpo_score",
            "docking_score", "sa_score", "herg_flag", "cyp_risk",
            "veber_pass", "admet_risk",
            "md_binding_energy", "md_rmsd_mean",
            "opt_score",
            "is_approved", "max_phase", "pref_name",
            "chembl_id", "chembl_url", "first_approval",
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]
        df_out = df[keep_cols].copy()

        df_out.to_csv(self.output_csv, index=False)
        logger.info("Optimized leads saved: %s (%d rows)", self.output_csv, len(df_out))
        return df_out
