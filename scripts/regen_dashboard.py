"""Regenerate benchmark + dashboard WITHOUT re-running Stage 1~4.

Usage (inside t2lead conda env):
  PYTHONUNBUFFERED=1 python scripts/regen_dashboard.py --disease "lung cancer"
  PYTHONUNBUFFERED=1 python scripts/regen_dashboard.py --disease "breast cancer"

Benchmark always runs full MD when ``benchmark.run_md`` is true in config (same as Stage 4).
"""
import argparse
import os
import shutil
from pathlib import Path
from copy import deepcopy

os.environ.setdefault("PYTHONUNBUFFERED", "1")

from drugpipe.config import load_config
from drugpipe.paths import stage_paths, STAGE4
from drugpipe.pipeline import _apply_target_specific_stage4_pdb, _resolve_target_chembl_id
from drugpipe.lead_optimization.protein_prep import ProteinPreparator
from drugpipe.benchmark import run_benchmark
from drugpipe.report import generate_report


TARGET_MAP = {
    "lung cancer": "CHEMBL203",
    "breast cancer": "CHEMBL4005",
}


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


def _organize_stage4_visual_assets(stage4_dir: Path) -> int:
    """Move legacy root-level PNG/SVG assets into a structured folder."""
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


def main():
    parser = argparse.ArgumentParser(description="Regenerate benchmark + dashboard only")
    parser.add_argument("--disease", required=True, help="e.g. 'lung cancer' or 'breast cancer'")
    parser.add_argument("--out-dir", default="/root/autodl-fs/T2Lead", help="Pipeline output root")
    parser.add_argument(
        "--no-organize-assets",
        action="store_true",
        help="Do not move legacy PNG/SVG files into stage4_optimization/visual_assets/",
    )
    args = parser.parse_args()

    disease = args.disease.strip()
    slug = disease.lower().replace(" ", "_")

    cfg = load_config(
        config_path="/root/T2Lead/configs/default_config.yaml",
        overrides={
            "pipeline": {"out_dir": args.out_dir},
            "target_discovery": {"disease": disease},
        },
    )

    run_root = Path(args.out_dir) / slug
    if not run_root.exists():
        raise FileNotFoundError(f"Run directory not found: {run_root}")

    layout = stage_paths(run_root, cfg)
    s4 = layout[STAGE4]

    target_id = TARGET_MAP.get(disease.lower())
    if not target_id:
        target_id = _resolve_target_chembl_id(cfg, layout)
    print(f"Disease: {disease}")
    print(f"Target: {target_id}")
    print(f"Stage4 dir: {s4}")

    cfg_s4 = deepcopy(cfg)
    cfg_s4 = _apply_target_specific_stage4_pdb(cfg_s4, layout, target_id)

    protein_info = ProteinPreparator(cfg_s4).prepare(s4) or {}
    print(f"Receptor: {protein_info.get('pdbqt_path', 'N/A')}")

    print("\n--- Running benchmark (find approved drugs + score) ---")
    out = run_benchmark(cfg_s4, target_id, protein_info, s4)
    n = 0 if out is None else len(out)
    print(f"Benchmark drugs found: {n}")
    if out is not None and not out.empty:
        for _, r in out.iterrows():
            name = r.get("pref_name", "?")
            dock = r.get("docking_score", "—")
            md = r.get("md_binding_energy", "—")
            rmsd = r.get("md_rmsd_mean", "—")
            opt = r.get("opt_score", 0)
            print(f"  {name:25s}  Dock={dock}  MD_dG={md}  RMSD={rmsd}  Opt={opt:.3f}")

    print("\n--- Generating dashboard.html ---")
    db = generate_report(cfg_s4, run_root, layout, stage4_dir=s4)
    print(f"Dashboard: {db}")
    if not args.no_organize_assets:
        moved = _organize_stage4_visual_assets(s4)
        print(f"Visual assets organized: moved {moved} item(s) into {s4 / 'visual_assets'}")
    print("Done!")


if __name__ == "__main__":
    main()
