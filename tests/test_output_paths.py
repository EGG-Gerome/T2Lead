"""Stage output path layout."""
from __future__ import annotations

from pathlib import Path

import pytest

from drugpipe.config import load_config
from drugpipe.paths import STAGE2, STAGE4, disease_slug, run_root_for_config, stage_paths


def test_disease_slug():
    assert disease_slug("Breast Cancer") == "breast_cancer"


def test_stage_paths_nested(tmp_path: Path):
    cfg = load_config(
        overrides={
            "pipeline": {
                "out_dir": str(tmp_path / "root"),
                "output_layout": {"use_stage_subdirs": True},
            },
            "target_discovery": {"disease": "test disease"},
        }
    )
    run_root = run_root_for_config(cfg)
    assert "test_disease" in str(run_root)
    layout = stage_paths(run_root, cfg)
    assert layout[STAGE2].name == "stage2_hits"
    assert layout[STAGE4].name == "stage4_optimization"
    assert (layout[STAGE2]).is_dir()


def test_stage_paths_flat(tmp_path: Path):
    cfg = load_config(
        overrides={
            "pipeline": {
                "out_dir": str(tmp_path / "root"),
                "output_layout": {"use_stage_subdirs": False},
            },
        }
    )
    run_root = tmp_path / "root"
    run_root.mkdir(parents=True)
    layout = stage_paths(run_root, cfg)
    assert layout[STAGE2] == run_root


def test_variant_run_id_env_keeps_underscore(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("DP_PIPELINE__OUTPUT_LAYOUT__VARIANT_RUN_ID", "20260413_191441")
    cfg = load_config(
        overrides={
            "pipeline": {"out_dir": str(tmp_path / "root")},
            "target_discovery": {"disease": "breast cancer"},
            "variant_analysis": {"enabled": True, "sample_id": "s1"},
        },
    )
    run_root = run_root_for_config(cfg)
    assert "20260413_191441" in str(run_root)
