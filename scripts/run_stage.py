#!/usr/bin/env python3
"""Run a single pipeline stage independently."""
# Module overview and key intent for maintainers.
# 模块总览与关键设计意图，便于后续维护。

# 独立运行流水线某一阶段。

import argparse
import copy
import logging
import os
import sys
from pathlib import Path

conda_lib = Path(sys.executable).resolve().parents[1] / "lib"
if conda_lib.exists():
    os.environ["LD_LIBRARY_PATH"] = str(conda_lib) + ":" + os.environ.get("LD_LIBRARY_PATH", "")

# 将 src 加入路径以便导入 drugpipe
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from drugpipe.config import load_config
from drugpipe.pipeline import run_hit_to_lead, run_target_discovery, run_target_to_hit
from drugpipe.paths import STAGE2, chembl_library_root, run_root_for_config, stage_paths


# main 的核心行为与设计意图。
def main() -> None:
    parser = argparse.ArgumentParser(description="Run a single DrugPipe stage.")
    parser.add_argument(
        "stage",
        choices=["target_discovery", "target_to_hit", "hit_to_lead"],
        help="Which stage to run.",
    )
    parser.add_argument("-c", "--config", default=None, help="Custom YAML config.")
    parser.add_argument("--disease", default=None, help="Disease name (Stage 1).")
    parser.add_argument("--target", default=None, help="ChEMBL target ID (Stage 2).")
    parser.add_argument("--activities-csv", default=None,
                        help="User-supplied IC50 CSV for novel targets.")
    parser.add_argument("--screening-library", default=None,
                        help="User-supplied SMILES library CSV for screening.")
    parser.add_argument("--docking-only", action="store_true",
                        help="Skip ML, use docking only (novel targets w/o IC50 data).")
    parser.add_argument("--hits-csv", default=None, help="Hits CSV path (Stage 3).")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    overrides = {}
    if args.disease:
        overrides.setdefault("target_discovery", {})["disease"] = args.disease
    if args.target:
        overrides.setdefault("target_to_hit", {})["target_chembl_id"] = args.target
    if args.activities_csv:
        overrides.setdefault("target_to_hit", {})["external_activities_csv"] = args.activities_csv
    if args.screening_library:
        overrides.setdefault("target_to_hit", {})["screening_library_csv"] = args.screening_library
    if args.docking_only:
        overrides.setdefault("target_to_hit", {})["docking_only"] = True

    cfg = load_config(config_path=args.config, overrides=overrides)

    if args.stage == "target_discovery":
        targets = run_target_discovery(cfg)
        for t in targets:
            print(f"  {t['chembl_id']}  {t.get('symbol', '')}  score={t.get('rank_score', t.get('score', 0)):.3f}")

    elif args.stage == "target_to_hit":
        run_root = run_root_for_config(cfg)
        layout = stage_paths(run_root, cfg)
        s2_local = layout[STAGE2]
        library_root = chembl_library_root(cfg, s2_local)
        cfg_s2 = copy.deepcopy(cfg)
        cfg_s2.setdefault("pipeline", {})["out_dir"] = str(s2_local)
        df_hits = run_target_to_hit(
            cfg_s2,
            target_chembl_id=args.target,
            crawl_out_dir=library_root,
        )
        print(f"\nHit candidates: {len(df_hits)} rows")
        print(df_hits.head(10).to_string(index=False))

    elif args.stage == "hit_to_lead":
        import pandas as pd
        run_root = run_root_for_config(cfg)
        default_hits = stage_paths(run_root, cfg)[STAGE2] / "final_hit_candidates.csv"
        hits_path = args.hits_csv or str(default_hits)
        df_hits = pd.read_csv(hits_path)
        print(f"Loaded {len(df_hits)} hits from {hits_path}")
        df_leads = run_hit_to_lead(cfg, df_hits)
        print(f"\nLead candidates: {len(df_leads)} rows")
        print(df_leads.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
