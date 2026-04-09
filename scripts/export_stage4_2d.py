#!/usr/bin/env python3
"""Export Stage 4 lead/benchmark 2D structure assets into visual_assets/.

Reads:
  - optimized_leads.csv
  - benchmark_drugs.csv (optional)

Writes:
  - visual_assets/leads/grids/*.svg|png
  - visual_assets/leads/structures/*.svg + index.tsv
  - visual_assets/benchmark/grids/*.svg|png
  - visual_assets/benchmark/structures/*.svg + index.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def _slug(text: str) -> str:
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(text or "").strip())
    return x.strip("_") or "item"


def _read_rows(csv_path: Path, kind: str) -> List[Dict[str, object]]:
    if not csv_path.exists():
        return []
    df = pd.read_csv(csv_path)
    out: List[Dict[str, object]] = []
    for i, row in df.reset_index(drop=True).iterrows():
        smi = str(row.get("canonical_smiles", "") or "").strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        AllChem.Compute2DCoords(mol)
        raw_name = row.get("pref_name", "")
        if pd.isna(raw_name):
            name = ""
        else:
            name = str(raw_name or "").strip()
        if (not name) or (name.lower() == "nan"):
            name = f"{kind}_{i+1:02d}"
        out.append(
            {
                "idx": i + 1,
                "name": name,
                "smiles": smi,
                "dock": row.get("docking_score"),
                "sa": row.get("sa_score"),
                "md": row.get("md_binding_energy"),
                "rmsd": row.get("md_rmsd_mean"),
                "mol": mol,
            }
        )
    return out


def _as_float(v: object, default: float = 0.0) -> float:
    try:
        if v is None or v == "":
            return default
        return float(v)
    except Exception:
        return default


def _fmt(v: object, digits: int) -> str:
    try:
        if v is None or v == "" or (isinstance(v, float) and pd.isna(v)):
            return "NA"
        return f"{float(v):.{digits}f}"
    except Exception:
        return "NA"


def _write_grids(rows: List[Dict[str, object]], out_dir: Path, prefix: str) -> None:
    if not rows:
        return
    mols = [r["mol"] for r in rows]
    legends_clean = [f"{prefix} #{int(r['idx'])}: {r['name']}" for r in rows]
    legends_metrics = [
        (
            f"{r['name']}\n"
            f"Dock {_fmt(r['dock'], 2)} | SA {_fmt(r['sa'], 2)}\n"
            f"MD {_fmt(r['md'], 1)} | RMSD {_fmt(r['rmsd'], 1)} A"
        )
        for r in rows
    ]

    out_dir.mkdir(parents=True, exist_ok=True)

    clean_svg = Draw.MolsToGridImage(
        mols,
        molsPerRow=2,
        subImgSize=(860, 560),
        legends=legends_clean,
        useSVG=True,
    )
    (out_dir / f"{prefix.lower()}_structures_2d_grid_clean.svg").write_text(clean_svg, encoding="utf-8")

    metric_svg = Draw.MolsToGridImage(
        mols,
        molsPerRow=2,
        subImgSize=(860, 600),
        legends=legends_metrics,
        useSVG=True,
    )
    (out_dir / f"{prefix.lower()}_structures_2d_grid_metrics.svg").write_text(metric_svg, encoding="utf-8")

    # PNG export is optional (depends on RDKit/Pillow build support).
    try:
        clean_png = Draw.MolsToGridImage(
            mols,
            molsPerRow=2,
            subImgSize=(860, 560),
            legends=legends_clean,
            useSVG=False,
        )
        clean_png.save(str(out_dir / f"{prefix.lower()}_structures_2d_grid_clean.png"))

        metric_png = Draw.MolsToGridImage(
            mols,
            molsPerRow=2,
            subImgSize=(860, 600),
            legends=legends_metrics,
            useSVG=False,
        )
        metric_png.save(str(out_dir / f"{prefix.lower()}_structures_2d_grid_metrics.png"))
    except Exception:
        pass


def _write_structures(rows: List[Dict[str, object]], out_dir: Path, prefix: str) -> None:
    if not rows:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.svg"):
        old.unlink()
    idx_path = out_dir / "index.tsv"
    with idx_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "name", "docking_score", "sa_score", "md_binding_energy", "md_rmsd_mean", "smiles"])
        for r in rows:
            file_base = f"{prefix.lower()}_{int(r['idx']):02d}_{_slug(str(r['name']))}"
            legend = (
                f"{r['name']}\n"
                f"Dock {_fmt(r['dock'], 2)} | SA {_fmt(r['sa'], 2)}\n"
                f"MD {_fmt(r['md'], 1)} | RMSD {_fmt(r['rmsd'], 1)} A"
            )
            drawer = Draw.MolDraw2DSVG(820, 620)
            drawer.DrawMolecule(r["mol"], legend=legend)
            drawer.FinishDrawing()
            (out_dir / f"{file_base}.svg").write_text(drawer.GetDrawingText(), encoding="utf-8")
            w.writerow([r["idx"], r["name"], r["dock"], r["sa"], r["md"], r["rmsd"], r["smiles"]])


def _resolve_stage4(args: argparse.Namespace) -> Path:
    if args.stage4_dir:
        return Path(args.stage4_dir)
    slug = args.disease.strip().lower().replace(" ", "_")
    return Path(args.out_dir) / slug / "stage4_optimization"


def main() -> None:
    parser = argparse.ArgumentParser(description="Export Stage4 2D structure assets into visual_assets/")
    parser.add_argument("--stage4-dir", default="", help="Direct stage4_optimization path")
    parser.add_argument("--out-dir", default="/root/autodl-fs/T2Lead", help="Pipeline output root")
    parser.add_argument("--disease", default="", help='Disease name, e.g. "breast cancer"')
    args = parser.parse_args()

    stage4 = _resolve_stage4(args)
    if not stage4.exists():
        raise FileNotFoundError(f"Stage4 directory not found: {stage4}")

    leads = _read_rows(stage4 / "optimized_leads.csv", "lead")
    bench = _read_rows(stage4 / "benchmark_drugs.csv", "benchmark")

    root = stage4 / "visual_assets"
    _write_grids(leads, root / "leads" / "grids", "lead")
    _write_structures(leads, root / "leads" / "structures", "lead")
    _write_grids(bench, root / "benchmark" / "grids", "benchmark")
    _write_structures(bench, root / "benchmark" / "structures", "benchmark")

    print(f"Stage4 dir: {stage4}")
    print(f"Lead structures exported: {len(leads)}")
    print(f"Benchmark structures exported: {len(bench)}")
    print(f"Output root: {root}")


if __name__ == "__main__":
    main()

