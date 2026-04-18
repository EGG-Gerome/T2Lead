"""Tests for patient-level aggregation across per-variant Stage4 outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from drugpipe.variant_analysis.patient_aggregation import build_patient_recommendations


def _write_opt_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def test_build_patient_recommendations(tmp_path: Path) -> None:
    run_root = tmp_path / "run"
    s4_a = run_root / "stage4_optimization" / "EGFR_L858R"
    s4_b = run_root / "stage4_optimization" / "PIK3CA_H1047R"

    _write_opt_csv(
        s4_a / "optimized_leads.csv",
        [
            {
                "chembl_id": "CHEMBL1",
                "pref_name": "DrugA",
                "canonical_smiles": "CCO",
                "rank_score": 0.90,
                "md_reliable": True,
                "admet_risk": 0.1,
            },
            {
                "chembl_id": "CHEMBL2",
                "pref_name": "DrugB",
                "canonical_smiles": "CCC",
                "rank_score": 0.70,
                "md_reliable": True,
                "admet_risk": 0.2,
            },
        ],
    )
    _write_opt_csv(
        s4_b / "optimized_leads.csv",
        [
            {
                "chembl_id": "CHEMBL1",
                "pref_name": "DrugA",
                "canonical_smiles": "CCO",
                "rank_score": 0.80,
                "md_reliable": True,
                "admet_risk": 0.1,
            },
            {
                "chembl_id": "CHEMBL3",
                "pref_name": "DrugC",
                "canonical_smiles": "CCN",
                "rank_score": 0.65,
                "md_reliable": False,
                "admet_risk": 0.3,
            },
        ],
    )

    cfg = {
        "variant_analysis": {
            "driver_genes": ["EGFR", "PIK3CA"],
            "patient_aggregation": {
                "enabled": True,
                "top_k_per_variant": 2,
                "top_n_patient_recommendations": 20,
            },
        },
    }
    ctx = [
        {
            "variant_key": "EGFR::L858R",
            "gene": "EGFR",
            "mutation": "L858R",
            "impact": "HIGH",
            "method": "experimental_mutated_localopt(1M17)",
            "stage4_dir": str(s4_a),
        },
        {
            "variant_key": "PIK3CA::H1047R",
            "gene": "PIK3CA",
            "mutation": "H1047R",
            "impact": "HIGH",
            "method": "experimental_mutated_localopt(4OVU)",
            "stage4_dir": str(s4_b),
        },
    ]

    out = build_patient_recommendations(cfg, run_root, ctx)
    assert out is not None

    patient_csv = Path(str(out["patient_csv"]))
    notes_md = Path(str(out["notes_md"]))
    dashboard = Path(str(out["dashboard_html"]))

    assert patient_csv.is_file()
    assert notes_md.is_file()
    assert dashboard.is_file()

    df = pd.read_csv(patient_csv)
    assert not df.empty
    # CHEMBL1 appears in two variants and should be retained.
    hit = df[df["chembl_id"] == "CHEMBL1"]
    assert not hit.empty
    assert int(hit.iloc[0]["support_variants"]) >= 2
