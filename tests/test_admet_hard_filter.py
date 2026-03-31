"""Tests for Stage 4 ADMET hard filtering."""
from __future__ import annotations

import pandas as pd

from drugpipe.config import load_config
from drugpipe.lead_optimization.admet_deep import DeepADMET


def _base_cfg():
    return load_config(overrides={"lead_optimization": {"enabled": True}})


def test_hard_filter_drops_herg():
    cfg = _base_cfg()
    ad = DeepADMET(cfg)
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CC", "CCC", "CCCC"],
            "herg_flag": [True, False, False],
            "veber_pass": [True, True, True],
            "cyp_risk": [False, False, True],
            "sa_score": [3.0, 3.0, 3.0],
        }
    )
    out = ad.hard_filter(df)
    assert len(out) == 2
    assert not out["herg_flag"].any()


def test_hard_filter_drops_veber_and_high_sa():
    cfg = _base_cfg()
    ad = DeepADMET(cfg)
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CC", "CCC"],
            "herg_flag": [False, False],
            "veber_pass": [False, True],
            "cyp_risk": [False, False],
            "sa_score": [3.0, 8.0],
        }
    )
    out = ad.hard_filter(df)
    assert len(out) == 0


def test_hard_filter_cyp_survives_and_sets_soft_risk():
    cfg = _base_cfg()
    ad = DeepADMET(cfg)
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CC"],
            "herg_flag": [False],
            "veber_pass": [True],
            "cyp_risk": [True],
            "sa_score": [4.0],
        }
    )
    out = ad.hard_filter(df)
    assert len(out) == 1
    assert out.loc[0, "admet_risk"] == 0.15


def test_hard_filter_disabled_noop():
    cfg = _base_cfg()
    cfg["lead_optimization"]["admet_deep"]["hard_filter"] = {"enabled": False}
    ad = DeepADMET(cfg)
    df = pd.DataFrame(
        {
            "canonical_smiles": ["CC"],
            "herg_flag": [True],
            "veber_pass": [True],
            "cyp_risk": [False],
            "sa_score": [3.0],
        }
    )
    out = ad.hard_filter(df)
    assert len(out) == 1
