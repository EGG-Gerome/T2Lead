#!/usr/bin/env python3
"""Regenerate dashboard.html from existing CSV data — NO benchmark re-run.

Usage:
    python scripts/regen_dashboard_html_only.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from drugpipe.config import load_config
from drugpipe.paths import stage_paths, STAGE4
from drugpipe.report import generate_report

OUT_ROOT = Path("/root/autodl-fs/T2Lead")

TARGET_MAP = {
    "lung cancer": "CHEMBL203",
    "breast cancer": "CHEMBL4005",
}

for disease in ("lung cancer", "breast cancer"):
    slug = disease.replace(" ", "_")
    run_root = OUT_ROOT / slug
    cfg = load_config(
        config_path=str(ROOT / "configs" / "default_config.yaml"),
        overrides={
            "pipeline": {"out_dir": str(OUT_ROOT)},
            "target_discovery": {"disease": disease},
        },
    )
    layout = stage_paths(run_root, cfg)
    s4 = layout[STAGE4]
    print(f"\n--- {disease} ---")
    print(f"  Stage4: {s4}")
    db = generate_report(cfg, run_root, layout, stage4_dir=s4)
    print(f"  Dashboard: {db}")

print("\nDone.")
