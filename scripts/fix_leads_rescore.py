#!/usr/bin/env python3
"""Fix Stage-4 optimized_leads.csv: overflow md_binding_energy → NaN, rescore with absolute scoring.

Run once after the MM-GBSA / scoring code fix and before dashboard generation:

    conda activate t2lead
    python scripts/fix_leads_rescore.py

Does NOT re-run MD; it cleans overflow values and recalculates opt_score using
the new absolute-range transforms (same scale as benchmark drugs).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from drugpipe.config import load_config
from drugpipe.lead_optimization.lead_optimizer import LeadOptimizer


OUT_ROOT = Path("/root/autodl-fs/T2Lead")
DISEASES = ["lung_cancer", "breast_cancer"]


def _fix_one(disease_slug: str) -> None:
    s4 = OUT_ROOT / disease_slug / "stage4_optimization"
    csv_path = s4 / "optimized_leads.csv"
    if not csv_path.is_file():
        print(f"[{disease_slug}] No optimized_leads.csv found at {csv_path}, skipping.")
        return

    df = pd.read_csv(csv_path)
    n = len(df)
    print(f"\n[{disease_slug}] Loaded {n} leads from {csv_path}")

    overflow_mask = df["md_binding_energy"].notna() & (df["md_binding_energy"] > 0)
    n_overflow = int(overflow_mask.sum())
    if n_overflow:
        print(f"  Clearing {n_overflow} overflow md_binding_energy values (were +ve)")
        df.loc[overflow_mask, "md_binding_energy"] = np.nan
        df.loc[overflow_mask, "md_binding_energy_std"] = np.nan

    extreme_mask = df["md_binding_energy"].notna() & (df["md_binding_energy"].abs() > 5000)
    n_extreme = int(extreme_mask.sum())
    if n_extreme:
        print(f"  Clearing {n_extreme} extreme md_binding_energy values (|dG|>5000)")
        df.loc[extreme_mask, "md_binding_energy"] = np.nan
        df.loc[extreme_mask, "md_binding_energy_std"] = np.nan

    cfg = load_config(
        config_path=str(ROOT / "configs" / "default_config.yaml"),
        overrides={"pipeline": {"out_dir": str(OUT_ROOT)}},
    )
    optimizer = LeadOptimizer(cfg, s4)
    df = optimizer._composite_score(df)

    df = df.sort_values("opt_score", ascending=False).reset_index(drop=True)

    bak = csv_path.with_suffix(".csv.bak_overflow")
    if not bak.exists():
        import shutil
        shutil.copy2(csv_path, bak)
        print(f"  Backup: {bak}")

    df.to_csv(csv_path, index=False)
    print(f"  Saved rescored leads → {csv_path}")

    for _, row in df.head(5).iterrows():
        smi = str(row.get("canonical_smiles", ""))[:40]
        dock = row.get("docking_score", np.nan)
        md = row.get("md_binding_energy", np.nan)
        rmsd = row.get("md_rmsd_mean", np.nan)
        opt = row.get("opt_score", np.nan)
        print(f"    {smi:40s}  dock={dock:>8.3f}  md_dG={'NaN' if pd.isna(md) else f'{md:.1f}':>8s}  rmsd={rmsd:>6.2f}  opt={opt:.4f}")


def main() -> None:
    for slug in DISEASES:
        _fix_one(slug)
    print("\nDone. Leads CSVs rescored with absolute-range transforms.")


if __name__ == "__main__":
    main()
