#!/usr/bin/env python3
"""Run a small A/B/C MD validation on selected Lung molecules.

Conditions:
  A) implicit baseline
  B) implicit + ligand positional restraint
  C) explicit solvent
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List

import pandas as pd

from drugpipe.config import load_config
from drugpipe.lead_optimization.md_simulation import ExplicitSolventRefiner, MDSimulator


def _parse_mol_ids(raw: str) -> List[int]:
    out: List[int] = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        out.append(int(token))
    if not out:
        raise ValueError("No molecule IDs parsed.")
    return out


def _build_input_df(stage3_csv: Path, stage4_csv: Path, poses_dir: Path, mol_ids: List[int]) -> pd.DataFrame:
    df3 = pd.read_csv(stage3_csv)
    df4 = pd.read_csv(stage4_csv)
    dock_map = (
        df4.dropna(subset=["canonical_smiles", "docking_score"])
        .drop_duplicates(subset=["canonical_smiles"], keep="first")
        .set_index("canonical_smiles")["docking_score"]
        .to_dict()
    )

    rows = []
    for mid in mol_ids:
        if mid < 0 or mid >= len(df3):
            raise IndexError(f"mol_id {mid} out of range for {stage3_csv} ({len(df3)} rows)")
        smi = str(df3.iloc[mid]["canonical_smiles"])
        pose = poses_dir / f"pose_{mid}.pdbqt"
        if not pose.is_file():
            raise FileNotFoundError(f"Docking pose not found: {pose}")
        rows.append(
            {
                "mol_id": mid,
                "canonical_smiles": smi,
                "docking_score": float(dock_map.get(smi, -7.0)),
                "docking_pose_file": str(pose),
                "opt_score": 1.0,  # for explicit top-N selection
            }
        )

    df = pd.DataFrame(rows).set_index("mol_id", drop=False)
    return df


def _run_implicit(
    df_in: pd.DataFrame,
    protein_info: Dict[str, str],
    out_dir: Path,
    restraint_k: float,
) -> pd.DataFrame:
    overrides = {
        "lead_optimization": {
            "md_simulation": {
                "enabled": True,
                "top_n_for_md": int(len(df_in)),
                "ligand_restraint_k_kcal_per_a2": float(restraint_k),
                "ligand_restraint_heavy_only": True,
                "ensemble": {
                    "enabled": True,
                    "n_runs": 1,
                    "equilibration_ps": 20,
                    "production_ps": 50,
                    "sample_interval_ps": 5,
                },
            }
        }
    }
    cfg = load_config(overrides=overrides)
    md = MDSimulator(cfg)
    return md.run(df_in.copy(), protein_info, out_dir)


def _run_explicit(
    df_in: pd.DataFrame,
    protein_info: Dict[str, str],
    out_dir: Path,
) -> pd.DataFrame:
    overrides = {
        "lead_optimization": {
            "explicit_md": {
                "enabled": True,
                "top_n": int(len(df_in)),
                "production_ns": 0.2,
                "equilibration_nvt_ps": 20,
                "equilibration_npt_ps": 20,
                "padding_nm": 1.0,
                "ionic_strength_M": 0.15,
            }
        }
    }
    cfg = load_config(overrides=overrides)
    refiner = ExplicitSolventRefiner(cfg)
    return refiner.run(df_in.copy(), protein_info, out_dir)


def main() -> None:
    ap = argparse.ArgumentParser(description="Run 2-molecule MD mini validation.")
    ap.add_argument(
        "--stage4-dir",
        default="/autodl-fs/data/T2Lead/lung_cancer/stage4_optimization",
        help="Stage4 optimization directory with receptor/pdbqt/poses.",
    )
    ap.add_argument(
        "--stage3-csv",
        default="/autodl-fs/data/T2Lead/lung_cancer/stage3_leads/final_lead_candidates.csv",
        help="Stage3 lead candidates CSV.",
    )
    ap.add_argument(
        "--stage4-csv",
        default="/autodl-fs/data/T2Lead/lung_cancer/stage4_optimization/optimized_leads.csv",
        help="Stage4 optimized leads CSV (for docking-score lookup).",
    )
    ap.add_argument(
        "--mol-ids",
        default="48,7",
        help="Comma-separated molecule IDs from stage3 rows.",
    )
    ap.add_argument(
        "--restraint-k",
        type=float,
        default=1.0,
        help="Implicit restraint strength in kcal/mol/A^2 for condition B.",
    )
    ap.add_argument(
        "--out-dir",
        default="",
        help="Output directory for mini-validation artifacts. Default: <stage4-dir>/mini_validation/<timestamp>",
    )
    args = ap.parse_args()

    stage4_dir = Path(args.stage4_dir).resolve()
    stage3_csv = Path(args.stage3_csv).resolve()
    stage4_csv = Path(args.stage4_csv).resolve()
    poses_dir = stage4_dir / "docking_poses"
    receptor_pdb = stage4_dir / "1M17_fixed.pdb"
    if not receptor_pdb.is_file():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")

    mol_ids = _parse_mol_ids(args.mol_ids)
    df_in = _build_input_df(stage3_csv, stage4_csv, poses_dir, mol_ids)
    protein_info = {"pdb_path": str(receptor_pdb)}

    stamp = time.strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (stage4_dir / "mini_validation" / stamp)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    conditions = [
        ("implicit_baseline", lambda: _run_implicit(df_in, protein_info, out_dir / "implicit_baseline", 0.0)),
        ("implicit_restrained", lambda: _run_implicit(df_in, protein_info, out_dir / "implicit_restrained", args.restraint_k)),
        ("explicit_solvent", lambda: _run_explicit(df_in, protein_info, out_dir / "explicit_solvent")),
    ]

    for name, runner in conditions:
        t0 = time.perf_counter()
        df_out = runner()
        elapsed = time.perf_counter() - t0

        csv_path = out_dir / f"{name}.csv"
        df_out.to_csv(csv_path, index=False)

        for _, row in df_out.reset_index(drop=True).iterrows():
            summary_rows.append(
                {
                    "condition": name,
                    "mol_id": int(row.get("mol_id")),
                    "canonical_smiles": row.get("canonical_smiles"),
                    "docking_score": row.get("docking_score"),
                    "md_binding_energy": row.get("md_binding_energy"),
                    "md_binding_energy_std": row.get("md_binding_energy_std"),
                    "md_rmsd_mean": row.get("md_rmsd_mean"),
                    "explicit_binding_energy": row.get("explicit_binding_energy"),
                    "explicit_rmsd_mean": row.get("explicit_rmsd_mean"),
                    "explicit_rmsd_drift": row.get("explicit_rmsd_drift"),
                    "condition_elapsed_s": elapsed,
                }
            )

        print(f"[done] {name} elapsed={elapsed:.1f}s csv={csv_path}")

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = out_dir / "summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    meta = {
        "timestamp": stamp,
        "mol_ids": mol_ids,
        "restraint_k_kcal_per_a2": float(args.restraint_k),
        "out_dir": str(out_dir),
        "summary_csv": str(summary_csv),
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[ok] summary={summary_csv}")


if __name__ == "__main__":
    main()
