"""Patient-level aggregation across per-variant Stage-4 outputs."""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd


def _to_num(v: Any) -> Optional[float]:
    try:
        x = float(v)
        if math.isnan(x) or math.isinf(x):
            return None
        return x
    except (TypeError, ValueError):
        return None


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def _pick_score_column(df: pd.DataFrame) -> Optional[str]:
    for col in ("rank_score", "opt_score", "fast_score"):
        if col in df.columns:
            return col
    return None


def _compound_key(row: pd.Series) -> str:
    cid = str(row.get("chembl_id", "") or "").strip()
    if cid:
        return cid
    smi = str(row.get("canonical_smiles", "") or "").strip()
    if smi:
        return hashlib.sha1(smi.encode("utf-8")).hexdigest()[:16]
    return ""


def _render_dashboard_html(payload: Dict[str, Any]) -> str:
    data = json.dumps(payload, ensure_ascii=False)
    return f"""<!DOCTYPE html>
<html lang="zh">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Patient-level Recommendations</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; background: #0f172a; color: #e2e8f0; }}
    h1 {{ margin-bottom: 8px; }}
    .sub {{ color: #94a3b8; margin-bottom: 18px; }}
    .card {{ background: #111827; border: 1px solid #334155; border-radius: 10px; padding: 14px; margin-bottom: 16px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
    th, td {{ border-bottom: 1px solid #334155; padding: 8px 6px; text-align: left; }}
    th {{ color: #93c5fd; }}
    .chart {{ height: 360px; }}
  </style>
</head>
<body>
  <h1>Patient-level Drug Recommendations</h1>
  <div class="sub" id="sub"></div>
  <div class="card"><div class="chart"><canvas id="c1"></canvas></div></div>
  <div class="card"><table id="tbl"><thead></thead><tbody></tbody></table></div>
<script>
const DATA = {data};
document.getElementById("sub").textContent =
  `run: ${{DATA.run_root}} | variants with Stage4: ${{DATA.n_variant_outputs}} | generated: ${{DATA.generated_at}}`;
const rows = DATA.top_recommendations || [];
const labels = rows.map((r, i) => `#${{i+1}} ${{r.compound_label || r.compound_key}}`);
const vals = rows.map(r => Number(r.patient_score || 0));
new Chart(document.getElementById("c1"), {{
  type: "bar",
  data: {{ labels, datasets: [{{ data: vals, backgroundColor: "rgba(14,165,233,.65)" }}] }},
  options: {{
    indexAxis: "y",
    responsive: true,
    maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{ x: {{ title: {{ display: true, text: "Patient score" }} }} }}
  }}
}});
document.querySelector("#tbl thead").innerHTML =
  "<tr><th>Rank</th><th>Compound</th><th>Patient score</th><th>Support variants</th><th>Genes</th><th>Chembl</th></tr>";
const tb = document.querySelector("#tbl tbody");
rows.forEach((r, i) => {{
  const tr = document.createElement("tr");
  tr.innerHTML = `<td>${{i+1}}</td><td>${{r.compound_label || ""}}</td><td>${{Number(r.patient_score||0).toFixed(4)}}</td><td>${{r.support_variants||0}}</td><td>${{r.genes||""}}</td><td>${{r.chembl_id||""}}</td>`;
  tb.appendChild(tr);
}});
</script>
</body>
</html>
"""


def build_patient_recommendations(
    cfg: Dict[str, Any],
    run_out_dir: Path,
    variant_contexts: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Aggregate per-variant optimized leads into one patient-level ranking."""
    if not variant_contexts:
        return None

    va = cfg.get("variant_analysis", {}) or {}
    agg_cfg = va.get("patient_aggregation", {}) or {}
    if not bool(agg_cfg.get("enabled", True)):
        return None

    top_k = max(1, int(agg_cfg.get("top_k_per_variant", 5)))
    top_n = max(1, int(agg_cfg.get("top_n_patient_recommendations", 50)))
    impact_weights = {
        "HIGH": float(agg_cfg.get("impact_weight_high", 1.0)),
        "MODERATE": float(agg_cfg.get("impact_weight_moderate", 0.75)),
        "LOW": float(agg_cfg.get("impact_weight_low", 0.45)),
    }
    md_reliable_weight = float(agg_cfg.get("md_unreliable_penalty_weight", 0.7))
    admet_risk_penalty = float(agg_cfg.get("admet_risk_penalty_weight", 0.5))

    driver_genes = {str(g).upper() for g in (va.get("driver_genes", []) or [])}
    driver_gene_weights_raw = agg_cfg.get("driver_gene_weights", {}) or {}
    driver_gene_weights = {str(k).upper(): float(v) for k, v in driver_gene_weights_raw.items()}

    gene_counts: Dict[str, int] = {}
    for ctx in variant_contexts:
        gene = str(ctx.get("gene", "")).strip().upper()
        if gene:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
    max_gene_count = max(gene_counts.values()) if gene_counts else 1

    per_variant_rows: List[Dict[str, Any]] = []
    for ctx in variant_contexts:
        stage4_dir = Path(str(ctx.get("stage4_dir", "")))
        opt_csv = stage4_dir / "optimized_leads.csv"
        if not opt_csv.is_file():
            continue

        try:
            df = pd.read_csv(opt_csv)
        except Exception:
            continue
        if df.empty:
            continue

        score_col = _pick_score_column(df)
        if score_col is None:
            continue

        df = df.copy()
        df["_base_score"] = pd.to_numeric(df[score_col], errors="coerce")
        df = df.dropna(subset=["_base_score"]).sort_values("_base_score", ascending=False).head(top_k)
        if df.empty:
            continue

        gene = str(ctx.get("gene", "")).strip()
        gene_u = gene.upper()
        impact = str(ctx.get("impact", "")).upper()
        impact_w = impact_weights.get(impact, impact_weights["MODERATE"])

        driver_w = driver_gene_weights.get(gene_u)
        if driver_w is None:
            driver_w = 1.25 if gene_u in driver_genes else 1.0

        gene_count = gene_counts.get(gene_u, 1)
        gene_freq_w = 1.0 + (math.log1p(gene_count) / max(math.log1p(max_gene_count), 1e-9))

        for _, row in df.iterrows():
            base = _to_num(row.get("_base_score"))
            if base is None:
                continue
            base = _clip01(base)

            md_rel = row.get("md_reliable", True)
            md_w = 1.0 if bool(md_rel) else md_reliable_weight

            admet = _to_num(row.get("admet_risk"))
            admet_w = 1.0 if admet is None else max(0.0, 1.0 - admet_risk_penalty * _clip01(admet))

            contrib = base * driver_w * impact_w * gene_freq_w * md_w * admet_w
            per_variant_rows.append(
                {
                    "variant_key": str(ctx.get("variant_key", "")),
                    "gene": gene,
                    "mutation": str(ctx.get("mutation", "")),
                    "impact": impact,
                    "method": str(ctx.get("method", "")),
                    "stage4_dir": str(stage4_dir),
                    "score_col": score_col,
                    "base_score": base,
                    "driver_weight": driver_w,
                    "impact_weight": impact_w,
                    "gene_frequency_weight": gene_freq_w,
                    "md_weight": md_w,
                    "admet_weight": admet_w,
                    "patient_score_contrib": contrib,
                    "chembl_id": row.get("chembl_id", ""),
                    "pref_name": row.get("pref_name", ""),
                    "canonical_smiles": row.get("canonical_smiles", ""),
                    "rank_score": row.get("rank_score", ""),
                    "opt_score": row.get("opt_score", ""),
                    "md_reliable": row.get("md_reliable", ""),
                    "admet_risk": row.get("admet_risk", ""),
                },
            )

    out_dir = run_out_dir / "patient_aggregation"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not per_variant_rows:
        return None

    df_var = pd.DataFrame(per_variant_rows).sort_values("patient_score_contrib", ascending=False)
    var_csv = out_dir / "variant_top_leads.csv"
    df_var.to_csv(var_csv, index=False)

    agg: Dict[str, Dict[str, Any]] = {}
    for _, row in df_var.iterrows():
        ck = _compound_key(row)
        if not ck:
            continue
        rec = agg.setdefault(
            ck,
            {
                "compound_key": ck,
                "chembl_id": str(row.get("chembl_id", "") or ""),
                "pref_name": str(row.get("pref_name", "") or ""),
                "canonical_smiles": str(row.get("canonical_smiles", "") or ""),
                "patient_score": 0.0,
                "support_variants": 0,
                "genes": set(),
                "mutations": set(),
                "best_rank_score": None,
                "mean_base_score": [],
            },
        )
        rec["patient_score"] += float(row.get("patient_score_contrib", 0.0) or 0.0)
        rec["support_variants"] += 1
        rec["genes"].add(str(row.get("gene", "") or ""))
        rec["mutations"].add(str(row.get("mutation", "") or ""))
        b = _to_num(row.get("base_score"))
        if b is not None:
            rec["mean_base_score"].append(b)
        r = _to_num(row.get("rank_score"))
        if r is not None and (rec["best_rank_score"] is None or r > rec["best_rank_score"]):
            rec["best_rank_score"] = r

    rows: List[Dict[str, Any]] = []
    for rec in agg.values():
        rows.append(
            {
                "compound_key": rec["compound_key"],
                "chembl_id": rec["chembl_id"],
                "pref_name": rec["pref_name"],
                "compound_label": rec["pref_name"] or rec["chembl_id"] or rec["compound_key"],
                "canonical_smiles": rec["canonical_smiles"],
                "patient_score": rec["patient_score"],
                "support_variants": rec["support_variants"],
                "genes": ";".join(sorted(g for g in rec["genes"] if g)),
                "mutations": ";".join(sorted(m for m in rec["mutations"] if m)),
                "best_rank_score": rec["best_rank_score"],
                "mean_base_score": (
                    sum(rec["mean_base_score"]) / len(rec["mean_base_score"])
                    if rec["mean_base_score"] else None
                ),
            },
        )

    df_pat = pd.DataFrame(rows).sort_values("patient_score", ascending=False).head(top_n)
    patient_csv = out_dir / "patient_recommendations.csv"
    df_pat.to_csv(patient_csv, index=False)

    top20 = df_pat.head(20)
    md_lines = [
        "# Patient-level 推荐说明",
        "",
        f"- run_root: `{run_out_dir}`",
        f"- 参与聚合的变异候选行数: `{len(df_var)}`",
        f"- 参与聚合的化合物数: `{len(df_pat)}`",
        "",
        "## 聚合方法（当前默认）",
        "",
        "- 每个变异取 top-K（按 `rank_score/opt_score/fast_score` 优先顺序）。",
        "- 贡献分 = 基础分 × 驱动基因权重 × 影响等级权重 × 基因频率权重 × MD 可靠性权重 × ADMET 风险权重。",
        "- 跨变异对同一化合物累加贡献分，得到 `patient_score`。",
        "",
        "## Top 20",
        "",
    ]
    if top20.empty:
        md_lines.append("无可用候选。")
    else:
        md_lines.append("| Rank | Compound | Patient score | Support variants | Genes |")
        md_lines.append("|---:|---|---:|---:|---|")
        for i, (_, r) in enumerate(top20.iterrows(), start=1):
            md_lines.append(
                f"| {i} | {r.get('compound_label', '')} | {float(r.get('patient_score', 0.0)):.4f} | "
                f"{int(r.get('support_variants', 0))} | {r.get('genes', '')} |",
            )
    notes_md = out_dir / "recommendation_notes_zh.md"
    notes_md.write_text("\n".join(md_lines), encoding="utf-8")

    payload = {
        "run_root": str(run_out_dir),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_variant_outputs": len({str(x.get("variant_key", "")) for x in per_variant_rows}),
        "top_recommendations": df_pat.head(50).to_dict(orient="records"),
    }
    (out_dir / "patient_report_payload.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    dashboard = out_dir / "dashboard.html"
    dashboard.write_text(_render_dashboard_html(payload), encoding="utf-8")

    return {
        "out_dir": str(out_dir),
        "variant_csv": str(var_csv),
        "patient_csv": str(patient_csv),
        "notes_md": str(notes_md),
        "dashboard_html": str(dashboard),
        "n_patient_rows": int(len(df_pat)),
    }

