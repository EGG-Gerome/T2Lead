"""Composite scoring (Stage 4) after ADMET hard gate."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from drugpipe.config import load_config
from drugpipe.lead_optimization.lead_optimizer import LeadOptimizer


def test_composite_score_weights(tmp_path: Path):
    cfg = load_config(overrides={"lead_optimization": {"enabled": True}})
    opt = LeadOptimizer(cfg, tmp_path)
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CCO", "CCC"],
            "docking_score": [-8.0, -7.0],
            "md_binding_energy": [-50.0, -40.0],
            "md_rmsd_mean": [1.0, 2.0],
            "cyp_risk": [False, True],
        }
    )
    out = opt._composite_score(df)
    assert "opt_score" in out.columns
    assert out["opt_score"].notna().all()
    # Lower docking (better) → higher score
    assert out.loc[0, "opt_score"] > out.loc[1, "opt_score"]


def test_composite_normalizes_when_md_missing(tmp_path: Path):
    cfg = load_config(overrides={"lead_optimization": {"enabled": True}})
    opt = LeadOptimizer(cfg, tmp_path)
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CC", "CCC"],
            "docking_score": [-9.0, -8.0],
            "md_binding_energy": [np.nan, np.nan],
            "md_rmsd_mean": [np.nan, np.nan],
            "cyp_risk": [False, False],
        }
    )
    out = opt._composite_score(df)
    assert out["opt_score"].between(0.0, 1.0 + 1e-9).all()
    assert out.loc[0, "opt_score"] > out.loc[1, "opt_score"]
