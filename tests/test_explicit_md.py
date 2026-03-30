"""Tests for the explicit-solvent refinement trigger logic."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from drugpipe.lead_optimization.md_simulation import ExplicitSolventRefiner


def _make_cfg(**overrides) -> dict:
    """Build a minimal config for ExplicitSolventRefiner."""
    emd = {
        "enabled": True,
        "top_n": 3,
        "production_ns": 0.01,
        "trigger": {"delta_opt_score_threshold": 0.05},
    }
    emd.update(overrides)
    return {"lead_optimization": {"explicit_md": emd}, "pipeline": {"device": "cpu"}}


class TestTriggerLogic:
    def test_trigger_when_scores_close(self):
        cfg = _make_cfg()
        refiner = ExplicitSolventRefiner(cfg)
        df = pd.DataFrame({
            "opt_score": [0.80, 0.81, 0.82, 0.79, 0.805],
            "canonical_smiles": ["C"] * 5,
        })
        assert refiner.should_trigger(df) is True

    def test_no_trigger_when_scores_spread(self):
        cfg = _make_cfg()
        refiner = ExplicitSolventRefiner(cfg)
        df = pd.DataFrame({
            "opt_score": [0.90, 0.50, 0.30, 0.20, 0.10],
            "canonical_smiles": ["C"] * 5,
        })
        assert refiner.should_trigger(df) is False

    def test_no_trigger_when_disabled(self):
        cfg = _make_cfg(enabled=False)
        refiner = ExplicitSolventRefiner(cfg)
        df = pd.DataFrame({
            "opt_score": [0.80, 0.81, 0.82],
            "canonical_smiles": ["C"] * 3,
        })
        assert refiner.should_trigger(df) is False

    def test_no_trigger_without_opt_score(self):
        cfg = _make_cfg()
        refiner = ExplicitSolventRefiner(cfg)
        df = pd.DataFrame({"canonical_smiles": ["C"] * 3})
        assert refiner.should_trigger(df) is False

    def test_no_trigger_single_candidate(self):
        cfg = _make_cfg()
        refiner = ExplicitSolventRefiner(cfg)
        df = pd.DataFrame({
            "opt_score": [0.90],
            "canonical_smiles": ["C"],
        })
        assert refiner.should_trigger(df) is False

    def test_threshold_boundary(self):
        cfg = _make_cfg()
        refiner = ExplicitSolventRefiner(cfg)
        df = pd.DataFrame({
            "opt_score": [0.70, 0.85, 0.80],
            "canonical_smiles": ["C"] * 3,
        })
        assert refiner.should_trigger(df) is False

    def test_config_defaults(self):
        cfg = {"lead_optimization": {}, "pipeline": {"device": "cpu"}}
        refiner = ExplicitSolventRefiner(cfg)
        assert refiner.enabled is False
        assert refiner.top_n == 5
        assert refiner.production_ns == 10.0
