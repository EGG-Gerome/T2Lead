"""Score benchmark (approved) drugs with the same scoring stack as Stage 4."""
from __future__ import annotations

import copy
import logging
import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from drugpipe.hit_to_lead.mpo import MPOScorer
from drugpipe.lead_optimization.admet_deep import DeepADMET
from drugpipe.lead_optimization.docking import VinaDocking
from drugpipe.lead_optimization.lead_optimizer import LeadOptimizer
from drugpipe.lead_optimization.md_simulation import MDSimulator

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem

    _RDKIT_OK = True
except Exception:
    Chem = None  # type: ignore[assignment]
    _RDKIT_OK = False


def _recover_protein_info_from_stage4(cfg: Dict[str, Any], out_dir: Path) -> Dict[str, Any]:
    """Best-effort recovery of protein info from Stage4 artifacts (for post-hoc benchmark runs)."""
    receptor = next(iter(sorted(out_dir.glob("*_receptor.pdbqt"))), None)
    fixed = next(iter(sorted(out_dir.glob("*_fixed.pdb"))), None)
    if fixed is None:
        fixed = next(iter(sorted(out_dir.glob("*.pdb"))), None)
    if receptor is None or fixed is None:
        return {}

    bs = (cfg.get("lead_optimization", {}) or {}).get("binding_site", {}) or {}
    center = list(bs.get("center", [0.0, 0.0, 0.0]))
    box_size = list(bs.get("box_size", [25.0, 25.0, 25.0]))

    logs_dir = out_dir.parent / "logs"
    if logs_dir.is_dir():
        logs = sorted(logs_dir.glob("*_full.log"), key=lambda p: p.stat().st_mtime, reverse=True)
        pat = re.compile(
            r"Auto-detected binding site center=([-\d.]+),([-\d.]+),([-\d.]+)\s+box=([-\d.]+),([-\d.]+),([-\d.]+)"
        )
        for lp in logs:
            try:
                text = lp.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            m = pat.search(text)
            if m:
                center = [float(m.group(1)), float(m.group(2)), float(m.group(3))]
                box_size = [float(m.group(4)), float(m.group(5)), float(m.group(6))]
                break

    return {
        "pdb_path": str(fixed),
        "pdbqt_path": str(receptor),
        "center": center,
        "box_size": box_size,
    }


def _maybe_run_md_for_benchmark(
    cfg: Dict[str, Any],
    df: pd.DataFrame,
    protein_info: Dict[str, Any],
    out_dir: Path,
) -> pd.DataFrame:
    bm = cfg.get("benchmark", {}) or {}
    if not bm.get("run_md", True):
        return df
    if not protein_info:
        logger.warning("Benchmark MD requested but protein info missing; skipping benchmark MD.")
        return df

    cfg_md = copy.deepcopy(cfg)
    md_cfg = (
        cfg_md.setdefault("lead_optimization", {})
        .setdefault("md_simulation", {})
    )
    md_cfg["enabled"] = True
    md_cfg["top_n_for_md"] = int(bm.get("md_top_n", min(2, len(df))))
    # Benchmark rows often reuse small integer indices (0..N-1), so sharing the
    # lead-stage checkpoint file can silently graft lead MD results onto drugs.
    md_cfg["resume_from_checkpoint"] = False
    ens = md_cfg.setdefault("ensemble", {})
    ens["enabled"] = bool(bm.get("md_ensemble_enabled", True))
    ens["n_runs"] = int(bm.get("md_ensemble_n_runs", 1))
    ens["equilibration_ps"] = float(bm.get("md_equilibration_ps", 20))
    ens["production_ps"] = float(bm.get("md_production_ps", 50))
    ens["sample_interval_ps"] = float(bm.get("md_sample_interval_ps", 5))

    logger.info(
        "Benchmark MD enabled: top_n=%d, ensemble_runs=%d, production_ps=%.1f, resume=%s",
        int(md_cfg["top_n_for_md"]),
        int(ens["n_runs"]),
        float(ens["production_ps"]),
        bool(md_cfg["resume_from_checkpoint"]),
    )
    md = MDSimulator(cfg_md)
    df_md = df.copy()
    if "canonical_smiles" in df_md.columns:
        # MD parameterization is fragile for salts/disconnected species.
        # Run on the largest fragment, then map MD outputs back.
        df_md["canonical_smiles"] = df_md["canonical_smiles"].astype(str).map(_largest_fragment_smiles)
    md_out = md.run(df_md, protein_info, out_dir)
    for c in ("md_binding_energy", "md_binding_energy_std", "md_rmsd_mean"):
        if c in md_out.columns:
            df[c] = md_out[c]
    return df


def _cfg_for_benchmark_docking(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Use lighter docking defaults for benchmark to avoid long post-hoc runs."""
    bm = cfg.get("benchmark", {}) or {}
    cfg_dock = copy.deepcopy(cfg)
    dk = cfg_dock.setdefault("lead_optimization", {}).setdefault("docking", {})
    dk["enabled"] = True
    default_exh = int(dk.get("exhaustiveness", 32))
    dk["exhaustiveness"] = int(bm.get("docking_exhaustiveness", min(default_exh, 8)))
    dk["n_poses"] = int(bm.get("docking_n_poses", 1))
    if "docking_n_cpus" in bm:
        dk["n_cpus"] = int(bm.get("docking_n_cpus", 0))
    return cfg_dock


def _largest_fragment_smiles(smiles: str) -> str:
    """Convert salts/multi-component SMILES to largest fragment for docking."""
    if not smiles or not _RDKIT_OK:
        return smiles
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) <= 1:
            return smiles
        main = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        return Chem.MolToSmiles(main, isomericSmiles=True)
    except Exception:
        return smiles


def score_benchmark_molecules(
    cfg: Dict[str, Any],
    drugs: List[Dict[str, Any]],
    protein_info: Dict[str, Any],
    out_dir: Path,
    model_predict_fn: Optional[Callable[..., Any]] = None,
    featurizer_fn: Optional[Callable[..., Any]] = None,
) -> Optional[pd.DataFrame]:
    """Score reference drugs with Dock/ADMET/MPO and optional MD.

    Defaults (see ``benchmark`` in config):
    - run_docking=True
    - run_md=True (full MD stack, same family as Stage 4 leads)
    """
    if not drugs:
        return None

    bm = cfg.get("benchmark", {}) or {}
    run_docking = bool(bm.get("run_docking", True))

    rows = []
    for d in drugs:
        rows.append(
            {
                "canonical_smiles": d["canonical_smiles"],
                "source_smiles": "__BENCHMARK__",
                "origin": "approved_benchmark",
                "pref_name": d.get("pref_name", ""),
                "chembl_id": d.get("chembl_id", ""),
                "chembl_url": d.get("chembl_url", ""),
                "max_phase": d.get("max_phase", np.nan),
                "first_approval": d.get("first_approval", np.nan),
                "reference_pchembl": d.get("best_pchembl", np.nan),
            }
        )
    df = pd.DataFrame(rows)

    if "pred_pIC50_ens" not in df.columns:
        df["pred_pIC50_ens"] = np.nan
    ref = df.get("reference_pchembl")
    if ref is not None:
        mask = df["pred_pIC50_ens"].isna() & ref.notna()
        df.loc[mask, "pred_pIC50_ens"] = ref[mask]

    mpo = MPOScorer(cfg)
    df = mpo.score(df, model_predict_fn=model_predict_fn, featurizer_fn=featurizer_fn)

    if ref is not None and df["pred_pIC50_ens"].isna().any():
        mask2 = df["pred_pIC50_ens"].isna() & ref.notna()
        df.loc[mask2, "pred_pIC50_ens"] = ref[mask2]

    admet = DeepADMET(cfg)
    df = admet.profile(df)

    if run_docking and not protein_info:
        protein_info = _recover_protein_info_from_stage4(cfg, out_dir)

    if run_docking and protein_info:
        cfg_dock = _cfg_for_benchmark_docking(cfg)
        docker = VinaDocking(cfg_dock)
        bench_pose_dir = out_dir / "benchmark_docking_poses"
        bench_pose_dir.mkdir(parents=True, exist_ok=True)
        # Keep original display SMILES, but dock with largest fragment when salts are present.
        df_dock = df.copy()
        if "canonical_smiles" in df_dock.columns:
            df_dock["canonical_smiles"] = df_dock["canonical_smiles"].astype(str).map(_largest_fragment_smiles)
        docked = docker.dock_leads(df_dock, protein_info, bench_pose_dir)
        df["docking_score"] = docked.get("docking_score")
        if "docking_pose_file" in docked.columns:
            df["docking_pose_file"] = docked["docking_pose_file"]
    else:
        df = df.copy()
        df["docking_score"] = np.nan

    if "md_binding_energy" not in df.columns:
        df["md_binding_energy"] = np.nan
    if "md_binding_energy_std" not in df.columns:
        df["md_binding_energy_std"] = np.nan
    if "md_rmsd_mean" not in df.columns:
        df["md_rmsd_mean"] = np.nan

    df = _maybe_run_md_for_benchmark(cfg, df, protein_info, out_dir)

    optimizer = LeadOptimizer(cfg, out_dir)
    df = optimizer._composite_score(df)
    df = optimizer._attach_ranking_signals(df)

    df["is_approved"] = True
    df["herg_flag"] = df.get("herg_flag", False)
    df["cyp_risk"] = df.get("cyp_risk", False)
    df["veber_pass"] = df.get("veber_pass", True)

    keep_cols = [
        "canonical_smiles",
        "origin",
        "pred_pIC50_ens",
        "QED",
        "mpo_score",
        "docking_score",
        "sa_score",
        "herg_flag",
        "cyp_risk",
        "veber_pass",
        "admet_risk",
        "md_binding_energy",
        "md_binding_energy_std",
        "md_rmsd_mean",
        "fast_score",
        "md_reliable",
        "rank_score",
        "opt_score",
        "is_approved",
        "max_phase",
        "pref_name",
        "chembl_id",
        "chembl_url",
        "first_approval",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df_out = df[keep_cols].copy()
    sort_col = "rank_score" if "rank_score" in df_out.columns else "opt_score"
    df_out = df_out.sort_values(sort_col, ascending=False).reset_index(drop=True)

    out_csv = out_dir / "benchmark_drugs.csv"
    df_out.to_csv(out_csv, index=False)
    logger.info("Benchmark drugs saved: %s (%d rows)", out_csv, len(df_out))
    return df_out
