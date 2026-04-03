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
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd

from drugpipe.lead_optimization.admet_deep import DeepADMET
from drugpipe.lead_optimization.approved_drug_checker import ApprovedDrugChecker
from drugpipe.lead_optimization.docking import VinaDocking
from drugpipe.lead_optimization.md_simulation import ExplicitSolventRefiner, MDSimulator
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
    4b. Explicit solvent MD (optional; triggered when top scores are too close)
    5. Composite scoring & ranking
    """

    def __init__(self, cfg: Dict[str, Any], out_dir: Path):
        self.cfg = cfg
        self.out_dir = out_dir
        lo = cfg.get("lead_optimization", {})
        self.enabled = bool(lo.get("enabled", True))
        self.top_n = int(lo.get("output", {}).get("top_n_optimized", 10))

        sc = lo.get("scoring", {})
        self.w_docking = float(sc.get("w_docking", 0.25))
        self.w_md = float(sc.get("w_md_energy", 0.40))
        self.w_stability = float(sc.get("w_stability", 0.35))
        self.w_cyp_soft = float(sc.get("w_cyp_soft", 0.08))

        rel = lo.get("md_reliability", {})
        self.md_reliable_delta_g_max = float(rel.get("delta_g_max_kcal_per_mol", -1.0))
        self.md_reliable_rmsd_max = float(rel.get("rmsd_max_angstrom", 15.0))
        self.md_fast_fallback = bool(rel.get("use_fast_score_fallback", True))

        self.protein_prep = ProteinPreparator(cfg)
        self.docker = VinaDocking(cfg)
        self.admet = DeepADMET(cfg)
        self.md_sim = MDSimulator(cfg)
        self.explicit_refiner = ExplicitSolventRefiner(cfg)
        self.drug_checker = ApprovedDrugChecker(cfg)

        self.output_csv = out_dir / "optimized_leads.csv"

    # ------------------------------------------------------------------
    def run(self, df_leads: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Execute the full lead-optimization pipeline.

        Returns
        -------
        (optimized_leads, protein_info)
            *protein_info* is empty when preparation was skipped or stage disabled.
        """
        if not self.enabled:
            logger.info("Lead optimization stage disabled.")
            return df_leads, {}

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
        df_leads = self.admet.hard_filter(df_leads)
        if df_leads.empty:
            logger.warning("All lead candidates removed by ADMET hard filters.")
            return df_leads, protein_info

        logger.info("=== Lead Optimization Step 4: MD Simulation ===")
        if protein_info:
            df_leads = self.md_sim.run(df_leads, protein_info, self.out_dir)
        else:
            logger.info("Skipping MD — no protein structure.")
            df_leads["md_binding_energy"] = np.nan
            df_leads["md_rmsd_mean"] = np.nan

        logger.info("=== Lead Optimization Step 5: Composite Scoring ===")
        df_leads = self._composite_score(df_leads)

        # Step 4b — optional explicit solvent refinement
        if protein_info and self.explicit_refiner.should_trigger(df_leads):
            logger.info("=== Lead Optimization Step 4b: Explicit Solvent Refinement ===")
            df_leads = self.explicit_refiner.run(df_leads, protein_info, self.out_dir)
            df_leads = self._rescore_with_explicit(df_leads)

        df_out = self._rank_and_output(df_leads)
        return df_out, protein_info

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

        # Soft CYP signal only (hERG / Veber / SA are hard-filtered upstream).
        if self.w_cyp_soft > 1e-9 and "cyp_risk" in df.columns:
            cyp = df["cyp_risk"].fillna(False).astype(bool)
            cyp_norm = (~cyp).astype(float) + 0.55 * cyp.astype(float)
        else:
            cyp_norm = None

        md_s = df.get("md_binding_energy")
        md_available = md_s is not None and md_s.notna().any()
        md_norm = self._md_binding_abs_norm(md_s) if md_available else None

        rmsd_s = df.get("md_rmsd_mean")
        stab_available = rmsd_s is not None and rmsd_s.notna().any()
        stab_norm = self._rmsd_abs_norm(rmsd_s) if stab_available else None

        components: list[tuple[float, pd.Series]] = []
        if cyp_norm is not None:
            components.append((self.w_cyp_soft, cyp_norm))
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
    def _clip01_series(s: pd.Series) -> pd.Series:
        return s.clip(lower=0.0, upper=1.0)

    def _fast_score(self, df: pd.DataFrame) -> pd.Series:
        """Fast non-MD composite score used as fallback when MD is unreliable."""
        idx = df.index
        num_sum = pd.Series(0.0, index=idx, dtype=float)
        den_sum = pd.Series(0.0, index=idx, dtype=float)

        def num_col(name: str) -> pd.Series:
            if name in df.columns:
                return pd.to_numeric(df[name], errors="coerce")
            return pd.Series(np.nan, index=idx, dtype=float)

        def add_component(weight: float, values: pd.Series) -> None:
            mask = values.notna()
            if not mask.any():
                return
            num_sum.loc[mask] += float(weight) * values.loc[mask]
            den_sum.loc[mask] += float(weight)

        pic = num_col("pred_pIC50_ens")
        add_component(0.30, self._clip01_series((pic - 5.0) / 6.0))

        qed = num_col("QED")
        add_component(0.10, self._clip01_series(qed))

        mpo = num_col("mpo_score")
        add_component(0.20, self._clip01_series(mpo))

        dock = num_col("docking_score")
        add_component(0.20, self._clip01_series(((-dock) - 4.0) / 6.0))

        sa = num_col("sa_score")
        add_component(0.10, self._clip01_series((6.0 - sa) / 5.0))

        admet = num_col("admet_risk")
        add_component(0.10, self._clip01_series(1.0 - admet))

        out = num_sum / den_sum.replace(0.0, np.nan)
        return out

    def _compute_md_reliable(self, df: pd.DataFrame) -> pd.Series:
        """Flag rows whose MD signal is considered reliable for ranking."""
        idx = df.index
        md = (
            pd.to_numeric(df["md_binding_energy"], errors="coerce")
            if "md_binding_energy" in df.columns
            else pd.Series(np.nan, index=idx, dtype=float)
        )
        rmsd = (
            pd.to_numeric(df["md_rmsd_mean"], errors="coerce")
            if "md_rmsd_mean" in df.columns
            else pd.Series(np.nan, index=idx, dtype=float)
        )
        return (
            md.notna()
            & rmsd.notna()
            & (md <= self.md_reliable_delta_g_max)
            & (rmsd <= self.md_reliable_rmsd_max)
        )

    def _attach_ranking_signals(self, df: pd.DataFrame) -> pd.DataFrame:
        """Attach fast_score / md_reliable / rank_score columns."""
        out = df.copy()
        out["fast_score"] = self._fast_score(out)
        out["md_reliable"] = self._compute_md_reliable(out).astype(bool)
        out["rank_score"] = out["opt_score"]
        if self.md_fast_fallback:
            use_fb = (~out["md_reliable"]) & out["fast_score"].notna()
            out.loc[use_fb, "rank_score"] = out.loc[use_fb, "fast_score"]
        return out

    @staticmethod
    def _dock_score_abs_norm(s: pd.Series) -> pd.Series:
        """kcal/mol: more negative is better. Maps roughly [-10, -4] → [1, 0]."""
        v = pd.to_numeric(s, errors="coerce").astype(float)
        return LeadOptimizer._clip01_series((-v - 4.0) / 6.0)

    @staticmethod
    def _md_binding_abs_norm(s: pd.Series) -> pd.Series:
        """kcal/mol: more negative is better. Maps roughly [-80, 0] → [1, 0]; +ΔG → 0."""
        v = pd.to_numeric(s, errors="coerce").astype(float)
        return LeadOptimizer._clip01_series((-v) / 80.0)

    @staticmethod
    def _rmsd_abs_norm(s: pd.Series) -> pd.Series:
        """Å: lower is better. Maps roughly [0, 15] → [1, 0]."""
        v = pd.to_numeric(s, errors="coerce").astype(float)
        return LeadOptimizer._clip01_series((15.0 - v) / 15.0)

    @staticmethod
    def _normalise_lower_better(s: Optional[pd.Series]) -> pd.Series:
        """Min-max normalise where *lower* raw values → *higher* normalised
        scores (closer to 1.0).  Used for explicit-solvent drift tie-break only."""
        if s is None or s.isna().all():
            return pd.Series(0.5, index=s.index if s is not None else pd.RangeIndex(0))
        filled = s.fillna(s.max() if s.max() == s.max() else 0)
        mn, mx = filled.min(), filled.max()
        if mx - mn < 1e-9:
            return pd.Series(0.5, index=s.index)
        return 1.0 - (filled - mn) / (mx - mn)

    def _rescore_with_explicit(self, df: pd.DataFrame) -> pd.DataFrame:
        """Re-rank candidates that have explicit-solvent data.

        For molecules with explicit MD results, replace the implicit energy /
        RMSD in the composite score and add drift as a tie-breaker.  Molecules
        without explicit data keep their original scores.
        """
        has_explicit = df["explicit_binding_energy"].notna()
        if not has_explicit.any():
            return df

        df = df.copy()
        df.loc[has_explicit, "md_binding_energy"] = df.loc[has_explicit, "explicit_binding_energy"]
        df.loc[has_explicit, "md_rmsd_mean"] = df.loc[has_explicit, "explicit_rmsd_mean"]

        df = self._composite_score(df)

        drift = df.get("explicit_rmsd_drift")
        if drift is not None and drift.notna().any():
            drift_penalty = self._normalise_lower_better(drift)
            df["opt_score"] = df["opt_score"] * 0.9 + drift_penalty * 0.1

        logger.info(
            "Re-scored %d molecules with explicit-solvent data.",
            int(has_explicit.sum()),
        )
        return df

    # ------------------------------------------------------------------
    def _rank_and_output(self, df: pd.DataFrame) -> pd.DataFrame:
        """Sort by composite score, keep top-N, annotate approval status, write CSV."""
        df = self._attach_ranking_signals(df)
        if "md_reliable" in df.columns:
            n_rel = int(df["md_reliable"].sum())
            n_tot = int(len(df))
            n_fb = int(((~df["md_reliable"]) & df["fast_score"].notna()).sum())
            logger.info(
                "MD reliability: %d/%d rows reliable; fallback ranking applied on %d rows.",
                n_rel, n_tot, n_fb,
            )

        df = df.sort_values(["rank_score", "opt_score"], ascending=False)
        df = df.head(self.top_n).reset_index(drop=True)

        df = self.drug_checker.annotate(df)

        keep_cols = [
            "canonical_smiles", "origin",
            "pred_pIC50_ens", "QED", "mpo_score",
            "docking_score", "sa_score", "herg_flag", "cyp_risk",
            "veber_pass", "admet_risk",
            "md_binding_energy", "md_binding_energy_std", "md_rmsd_mean",
            "explicit_binding_energy", "explicit_rmsd_mean", "explicit_rmsd_drift",
            "fast_score", "md_reliable", "rank_score",
            "opt_score",
            "is_approved", "max_phase", "pref_name",
            "chembl_id", "chembl_url", "first_approval",
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]
        df_out = df[keep_cols].copy()

        df_out.to_csv(self.output_csv, index=False)
        logger.info("Optimized leads saved: %s (%d rows)", self.output_csv, len(df_out))
        return df_out
