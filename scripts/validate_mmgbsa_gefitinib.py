#!/usr/bin/env python3
"""Quick MM-GBSA sanity check: GEFITINIB + 1M17 with docked pose (same-frame ΔG).

Expect ``md_binding_energy`` roughly in [-120, 5] kcal/mol for this crude
implicit model (not experimental ΔG). Run inside ``t2lead`` conda env::

  conda activate t2lead
  python scripts/validate_mmgbsa_gefitinib.py

Override paths with env vars: MMDEMO_PDB, MMDEMO_POSE, MMDEMO_DEVICE.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from drugpipe.lead_optimization.md_simulation import MDSimulator  # noqa: E402


def main() -> int:
    pdb = os.environ.get(
        "MMDEMO_PDB",
        "/root/autodl-fs/T2Lead/lung_cancer/stage4_optimization/1M17_fixed.pdb",
    )
    pose = os.environ.get(
        "MMDEMO_POSE",
        "/root/autodl-fs/T2Lead/lung_cancer/stage4_optimization/"
        "benchmark_docking_poses/docking_poses/pose_1.pdbqt",
    )
    device = os.environ.get("MMDEMO_DEVICE", "auto")

    if not Path(pdb).is_file():
        print(f"Missing receptor PDB: {pdb}", file=sys.stderr)
        return 2
    if not Path(pose).is_file():
        print(f"Missing pose PDBQT: {pose}", file=sys.stderr)
        return 2

    cfg = {
        "pipeline": {"device": device},
        "lead_optimization": {
            "md_simulation": {
                "enabled": True,
                "top_n_for_md": 1,
                "equilibration_ps": 0,
                "production_ns": 0.002,
                "gb_model": "OBC2",
                "ensemble": {"enabled": False},
            }
        },
    }

    protein_info = {
        "pdb_path": pdb,
        "pdbqt_path": str(Path(pdb).with_name(Path(pdb).stem.replace("_fixed", "") + "_receptor.pdbqt")),
        "center": [-12.8, -9.2, 12.5],
        "box_size": [30.0, 20.0, 30.0],
    }
    rec_qt = Path(pdb).parent / (Path(pdb).stem.replace("_fixed", "") + "_receptor.pdbqt")
    if rec_qt.is_file():
        protein_info["pdbqt_path"] = str(rec_qt)

    smi = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"  # GEFITINIB
    df = pd.DataFrame(
        [
            {
                "canonical_smiles": smi,
                "docking_score": -7.7,
                "docking_pose_file": pose,
            }
        ]
    )

    out_dir = Path(os.environ.get("MMDEMO_OUT", "/tmp/validate_mmgbsa_gefitinib"))
    out_dir.mkdir(parents=True, exist_ok=True)

    md = MDSimulator(cfg)
    out = md.run(df, protein_info, out_dir)
    dg = out["md_binding_energy"].iloc[0]
    print("md_binding_energy_kcal_mol:", dg)
    if dg != dg:
        print("FAIL: MD returned NaN", file=sys.stderr)
        return 1
    dg_f = float(dg)
    # Wide band: catches +3500-style overflow while allowing crude GBSA variance
    if not (-120.0 < dg_f < 5.0):
        print(f"FAIL: dG={dg_f} outside sanity band (-120, 5)", file=sys.stderr)
        return 1
    print("PASS: same-frame MM-GBSA in expected numerical range.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
