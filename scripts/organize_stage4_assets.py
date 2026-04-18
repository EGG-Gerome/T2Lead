#!/usr/bin/env python3
"""Reorganize legacy Stage 4 PNG/SVG files into visual_assets/.

This script is safe to run repeatedly:
- missing legacy files are ignored
- existing destination files are replaced
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def _move_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    shutil.move(str(src), str(dst))
    return True


def organize_stage4_visual_assets(stage4_dir: Path) -> int:
    mapping = {
        "lead_structures_2d": "visual_assets/leads/structures",
        "benchmark_structures_2d": "visual_assets/benchmark/structures",
        "lead_structures_2d_grid.png": "visual_assets/leads/grids/lead_structures_2d_grid.png",
        "lead_structures_2d_grid.svg": "visual_assets/leads/grids/lead_structures_2d_grid.svg",
        "lead_structures_2d_grid_clean.png": "visual_assets/leads/grids/lead_structures_2d_grid_clean.png",
        "lead_structures_2d_grid_clean.svg": "visual_assets/leads/grids/lead_structures_2d_grid_clean.svg",
        "lead_structures_2d_grid_metrics.png": "visual_assets/leads/grids/lead_structures_2d_grid_metrics.png",
        "lead_structures_2d_grid_metrics.svg": "visual_assets/leads/grids/lead_structures_2d_grid_metrics.svg",
        "benchmark_structures_2d_grid_clean.png": "visual_assets/benchmark/grids/benchmark_structures_2d_grid_clean.png",
        "benchmark_structures_2d_grid_clean.svg": "visual_assets/benchmark/grids/benchmark_structures_2d_grid_clean.svg",
        "benchmark_structures_2d_grid_metrics.png": "visual_assets/benchmark/grids/benchmark_structures_2d_grid_metrics.png",
        "benchmark_structures_2d_grid_metrics.svg": "visual_assets/benchmark/grids/benchmark_structures_2d_grid_metrics.svg",
    }
    moved = 0
    for src_name, dst_rel in mapping.items():
        if _move_if_exists(stage4_dir / src_name, stage4_dir / dst_rel):
            moved += 1
    return moved


def _stage4_from_args(out_dir: str, disease: str) -> Path:
    slug = disease.strip().lower().replace(" ", "_")
    return Path(out_dir) / slug / "stage4_optimization"


def main() -> None:
    parser = argparse.ArgumentParser(description="Move legacy Stage 4 visual files into visual_assets/")
    parser.add_argument("--stage4-dir", default="", help="Direct stage4_optimization path")
    parser.add_argument("--out-dir", default="/root/autodl-fs/T2Lead", help="Pipeline output root")
    parser.add_argument("--disease", default="", help="Disease label, e.g. 'breast cancer'")
    args = parser.parse_args()

    if args.stage4_dir:
        stage4_dir = Path(args.stage4_dir)
    elif args.disease:
        stage4_dir = _stage4_from_args(args.out_dir, args.disease)
    else:
        raise SystemExit("Provide either --stage4-dir or --disease (with optional --out-dir).")

    if not stage4_dir.exists():
        raise FileNotFoundError(f"Stage 4 directory not found: {stage4_dir}")

    moved = organize_stage4_visual_assets(stage4_dir)
    print(f"Stage4 dir: {stage4_dir}")
    print(f"Moved {moved} item(s) into {stage4_dir / 'visual_assets'}")


if __name__ == "__main__":
    main()

