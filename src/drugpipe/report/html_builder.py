"""Build self-contained dashboard.html from pipeline artifacts."""
from __future__ import annotations

import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from drugpipe.paths import STAGE1, STAGE2, STAGE3, STAGE4
from drugpipe.report.i18n import DASHBOARD_I18N
from drugpipe.report.log_parser import find_latest_full_log, issues_to_json_serializable, parse_full_log


def _csv_rows(path: Path, max_file_bytes: int = 80_000_000) -> int:
    """Row count minus header; skip very large files (use log-derived funnel instead)."""
    if not path.is_file():
        return 0
    try:
        if path.stat().st_size > max_file_bytes:
            return 0
        with open(path, encoding="utf-8", errors="replace") as f:
            return max(0, sum(1 for _ in f) - 1)
    except OSError:
        return 0


def _df_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    if df is None or df.empty:
        return []
    out: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        d: Dict[str, Any] = {}
        for k, v in row.items():
            if isinstance(v, (float, np.floating)) and (math.isnan(v) or math.isinf(v)):
                d[k] = None
            elif v is pd.NA:
                d[k] = None
            elif isinstance(v, (np.integer, np.floating)):
                d[k] = float(v) if isinstance(v, np.floating) else int(v)
            elif isinstance(v, (bool, np.bool_)):
                d[k] = bool(v)
            else:
                d[k] = v
        out.append(d)
    return out


def _funnel_from_log(text: str, stage2: Path) -> Tuple[List[float], bool]:
    lib = 0
    m = re.search(r"scored_candidates\.csv\s*\((\d+)\s+rows?\)", text, re.I)
    if m:
        lib = int(m.group(1))
    if lib <= 0:
        lib = _csv_rows(stage2 / "scored_candidates.csv")

    def grab(pat: str) -> Optional[int]:
        mm = re.search(pat, text)
        return int(mm.group(1)) if mm else None

    a = grab(r"After pIC50[^:]*:\s*(\d+)\s*/")
    b = grab(r"After QED[^:]*:\s*(\d+)")
    c = grab(r"After removing structural alerts:\s*(\d+)")
    d = grab(r"After ADMET rules:\s*(\d+)")
    e = grab(r"Hit candidates saved:.*\((\d+)\s+rows?\)")
    if e is None:
        e = _csv_rows(stage2 / "final_hit_candidates.csv")

    vals = [float(lib or 0), float(a or 0), float(b or 0), float(c or 0), float(d or 0), float(e or 0)]
    ok = lib > 0 or any(x > 0 for x in vals[1:])
    return vals, ok


def _parse_cv_r2(text: str) -> Optional[float]:
    m = re.search(r"cv_r2_mean['\"]?\s*:\s*([\d.]+)", text)
    if m:
        return float(m.group(1))
    m = re.search(r"CV summary:.*R²=([\d.]+)\s*±", text)
    if m:
        return float(m.group(1))
    return None


def _parse_docking_ok(text: str) -> Tuple[Optional[int], Optional[int]]:
    m = re.search(r"Docking complete:\s*(\d+)/(\d+)", text)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def _parse_admet_risk(text: str) -> Optional[int]:
    m = re.search(r"Enhanced ADMET complete:\s*(\d+)\s*/\s*(\d+)\s+flagged", text)
    if m:
        return int(m.group(1))
    return None


def _parse_reinvent(text: str) -> Optional[int]:
    m = re.search(r"REINVENT4 produced\s+(\d+)", text)
    return int(m.group(1)) if m else None


def _parse_selected_target(text: str) -> Tuple[Optional[str], Optional[str]]:
    m = re.search(r"Selected target\s+(CHEMBL\d+)\s+\(([^)]+)\)", text)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def _safe_num(v: Any) -> Optional[float]:
    try:
        if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
            return None
        x = float(v)
        if math.isnan(x) or math.isinf(x):
            return None
        return x
    except (TypeError, ValueError):
        return None


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def _fast_score_row(row: pd.Series) -> Optional[float]:
    """Fast benchmark-like score from visible medicinal chemistry metrics.

    Uses pIC50/QED/MPO/Docking/SA/ADMET (no MD), with row-wise weight
    renormalization when some metrics are missing.
    """
    # Absolute transforms (not per-batch min-max) to keep cross-run comparability.
    parts: List[Tuple[float, Optional[float]]] = [
        (0.30, _safe_num(row.get("pred_pIC50_ens"))),  # transformed below
        (0.10, _safe_num(row.get("QED"))),
        (0.20, _safe_num(row.get("mpo_score"))),
        (0.20, _safe_num(row.get("docking_score"))),  # transformed below
        (0.10, _safe_num(row.get("sa_score"))),  # transformed below
        (0.10, _safe_num(row.get("admet_risk"))),  # transformed below
    ]

    # Transform to "higher is better" [0,1].
    transformed: List[Tuple[float, Optional[float]]] = []
    for i, (w, v) in enumerate(parts):
        if v is None:
            transformed.append((w, None))
            continue
        if i == 0:  # pIC50
            transformed.append((w, _clip01((v - 5.0) / 6.0)))
        elif i == 3:  # docking: more negative is better
            transformed.append((w, _clip01(((-v) - 4.0) / 6.0)))
        elif i == 4:  # SA: lower is better, typical useful range ~[1,6]
            transformed.append((w, _clip01((6.0 - v) / 5.0)))
        elif i == 5:  # admet_risk: lower is better
            transformed.append((w, _clip01(1.0 - v)))
        else:
            transformed.append((w, _clip01(v)))

    w_sum = sum(w for w, v in transformed if v is not None)
    if w_sum <= 1e-9:
        return None
    score = sum(w * float(v) for w, v in transformed if v is not None) / w_sum
    return float(score)


def _add_fast_score(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    out = df.copy()
    out["fast_score"] = out.apply(_fast_score_row, axis=1)
    return out


def build_payload(
    cfg: Dict[str, Any],
    run_root: Path,
    layout: Dict[str, Path],
    stage4_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    run_root = Path(run_root)
    s4 = Path(stage4_dir) if stage4_dir else layout[STAGE4]
    s2 = layout[STAGE2]
    s3 = layout[STAGE3]
    s1 = layout[STAGE1]
    logs_dir = layout.get("logs", run_root / "logs")

    log_path = find_latest_full_log(logs_dir)
    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path else ""

    disease = (cfg.get("target_discovery", {}) or {}).get("disease", "") or ""
    disease = str(disease).strip() or run_root.name.replace("_", " ")

    ranked_path = s1 / "ranked_targets.csv"
    targets_df = pd.DataFrame()
    if ranked_path.is_file():
        targets_df = pd.read_csv(ranked_path)

    target_tid = ""
    target_sym = "—"
    log_tid, log_sym = _parse_selected_target(log_text)
    if log_tid:
        target_tid = log_tid
    if log_sym:
        target_sym = log_sym
    if not targets_df.empty and "chembl_id" in targets_df.columns:
        target_tid = target_tid or str(targets_df.iloc[0].get("chembl_id", "") or "")
    if not targets_df.empty and "symbol" in targets_df.columns:
        target_sym = target_sym if target_sym != "—" else str(targets_df.iloc[0].get("symbol", "") or target_sym)

    lo = cfg.get("lead_optimization", {}) or {}
    recs = sorted(s4.glob("*_receptor.pdbqt"), key=lambda p: p.stat().st_mtime, reverse=True)
    pdb_id = recs[0].name.replace("_receptor.pdbqt", "") if recs else ""
    if not pdb_id:
        pdb_id = str(lo.get("pdb_id", "") or "").strip()
    if not pdb_id:
        pdb_id = "—"

    funnel_vals, funnel_ok = _funnel_from_log(log_text, s2)
    n_leads = _csv_rows(s3 / "final_lead_candidates.csv")
    n_reinvent = _parse_reinvent(log_text) or 0

    opt_path = s4 / "optimized_leads.csv"
    df_opt = pd.read_csv(opt_path) if opt_path.is_file() else pd.DataFrame()
    df_opt = _add_fast_score(df_opt)

    bench_path = s4 / "benchmark_drugs.csv"
    df_bench = pd.read_csv(bench_path) if bench_path.is_file() else pd.DataFrame()
    df_bench = _add_fast_score(df_bench)

    cv_r2 = _parse_cv_r2(log_text)
    dock_ok, dock_tot = _parse_docking_ok(log_text)
    if dock_tot is None:
        dock_tot = len(df_opt) if not df_opt.empty else None
    if dock_ok is None:
        dock_ok = dock_tot

    admet_flagged = _parse_admet_risk(log_text) or 0

    n_approved_leads = 0
    if not df_opt.empty and "is_approved" in df_opt.columns:
        n_approved_leads = int(df_opt["is_approved"].fillna(False).astype(bool).sum())

    issues = parse_full_log(log_path) if log_path else []

    n_tgts = len(targets_df) if not targets_df.empty else 5
    n_hits = int(funnel_vals[5]) if funnel_ok else _csv_rows(s2 / "final_hit_candidates.csv")

    return {
        "disease": disease,
        "target_chembl_id": target_tid,
        "target_symbol": target_sym,
        "pdb_id": pdb_id,
        "log_name": log_path.name if log_path else "",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cv_r2": cv_r2,
        "dock_ok": dock_ok,
        "dock_tot": dock_tot,
        "admet_flagged": admet_flagged,
        "n_approved_leads": n_approved_leads,
        "n_optimized": len(df_opt),
        "stage_counts": {
            "s1": n_tgts,
            "s2": n_hits,
            "s3": n_leads,
            "s3_extra": n_reinvent,
            "s4": len(df_opt),
        },
        "funnel_vals": funnel_vals,
        "funnel_ok": funnel_ok,
        "targets": _df_records(targets_df),
        "leads": _df_records(df_opt),
        "benchmark": _df_records(df_bench),
        "issues": issues_to_json_serializable(issues),
        "i18n": DASHBOARD_I18N,
    }


def render_html(payload: Dict[str, Any]) -> str:
    data_json = json.dumps(payload, ensure_ascii=False, default=str)
    template = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>T2Lead dashboard</title>
<link href="https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=Instrument+Serif&family=Noto+Sans+SC:wght@400;500;600&family=Space+Grotesk:wght@400;500;600&display=swap" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
:root{--bg:#0a0e17;--surface:#111827;--border:#1e2d4a;--text:#e8edf5;--muted:#8899b4;--a1:#00e5a0;--a2:#00b8d4;--warn:#ff6b6b;--gold:#fbbf24;--purple:#a78bfa}
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:Space Grotesk,Noto Sans SC,sans-serif;background:var(--bg);color:var(--text);padding:28px 32px;line-height:1.45}
.hdr{display:flex;justify-content:space-between;align-items:flex-start;gap:16px;flex-wrap:wrap}
h1{font-family:Instrument Serif,serif;font-size:2.4rem;font-weight:400;background:linear-gradient(135deg,var(--a1),var(--a2));-webkit-background-clip:text;-webkit-text-fill-color:transparent}
.sub{color:var(--muted);font-size:.88rem;font-family:DM Mono,monospace;margin-top:6px;max-width:960px}
.lang{display:flex;background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:3px;gap:2px}
.lang button{border:none;background:transparent;color:var(--muted);padding:8px 16px;border-radius:8px;cursor:pointer;font-family:DM Mono,monospace;font-size:.8rem}
.lang button.on{background:linear-gradient(135deg,var(--a1),var(--a2));color:var(--bg);font-weight:600}
.flow{display:flex;align-items:center;gap:0;flex-wrap:wrap;margin:24px 0}
.node{flex:1;min-width:120px;background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:14px;text-align:center}
.node .n{font-size:.62rem;color:var(--a1);letter-spacing:.1em;text-transform:uppercase;font-family:DM Mono}
.node .v{font-size:1.6rem;font-weight:600;color:var(--a2)}
.arr{color:var(--muted);padding:0 8px}
.grid4{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-bottom:22px}
@media(max-width:900px){.grid4{grid-template-columns:1fr 1fr}}
.card{background:var(--surface);border:1px solid var(--border);border-radius:14px;padding:18px;position:relative}
.card::before{content:'';position:absolute;top:0;left:0;right:0;height:2px;background:linear-gradient(90deg,var(--a1),var(--a2));opacity:.65}
.ct{font-size:.68rem;text-transform:uppercase;letter-spacing:.1em;color:var(--muted);font-family:DM Mono;margin-bottom:6px}
.val{font-size:1.9rem;font-weight:600}
.sec{font-family:Instrument Serif,serif;font-size:1.45rem;margin:28px 0 12px}
.rules{color:var(--muted);font-size:.88rem;margin:10px 0 16px;line-height:1.6}
table{width:100%;border-collapse:collapse;font-size:.78rem;font-family:DM Mono,monospace}
th,td{padding:8px 6px;border-bottom:1px solid rgba(30,45,74,.45);text-align:left}
th{color:var(--muted);font-size:.65rem;text-transform:uppercase;letter-spacing:.06em;position:sticky;top:0;background:var(--surface);z-index:1}
tr:hover td{background:#1a2236}
.badge{display:inline-block;padding:2px 8px;border-radius:6px;font-size:.65rem;font-weight:600}
.b-bench{background:rgba(0,229,160,.18);color:var(--a1);border:1px solid var(--a1)}
.b-lead{background:rgba(0,184,212,.12);color:var(--a2)}
.row-bench td:first-child{box-shadow:inset 3px 0 0 var(--a1)}
a.chembl{color:var(--a2);text-decoration:none;border-bottom:1px dotted var(--a2)}
a.chembl:hover{color:var(--a1)}
.chart{height:280px;margin-top:8px;position:relative}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:16px}
@media(max-width:800px){.grid2{grid-template-columns:1fr}}
.issue{border-left:3px solid var(--warn);padding:12px 14px;margin-bottom:10px;border-radius:0 8px 8px 0;background:#1a2236}
.issue.e{border-color:#f87171}
.issue.w{border-color:var(--gold)}
.issue.i{border-color:var(--a2)}
.issue h4{font-size:.9rem;margin-bottom:4px}
.issue p{color:var(--muted);font-size:.8rem}
details.benign{margin-top:14px;border:1px solid var(--border);border-radius:10px;padding:10px 14px;background:var(--surface)}
details.benign summary{cursor:pointer;color:var(--a2);font-size:.85rem}
.sm{max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:.68rem}
</style>
</head>
<body>
<div class="hdr">
  <div>
    <h1 id="title">T2Lead</h1>
    <p class="sub" id="subtitle"></p>
  </div>
  <div class="lang"><button type="button" class="on" id="btn-en">EN</button><button type="button" id="btn-zh">中文</button></div>
</div>
<div class="flow" id="flow"></div>
<div class="grid4" id="metrics"></div>
<div class="sec" id="sec-tg"></div>
<div class="card" style="margin-bottom:20px"><div style="max-height:360px;overflow:auto"><table id="tbl-targets"><thead></thead><tbody></tbody></table></div></div>
<div class="sec" id="sec-fun"></div>
<div class="card" style="margin-bottom:20px"><div class="chart"><canvas id="cFunnel"></canvas></div></div>
<div class="sec" id="sec-ch"></div>
<div class="rules" id="chart-rules"></div>
<div class="rules" id="score-note"></div>
<div class="grid2">
  <div class="card"><div class="ct" id="c1t"></div><div class="chart"><canvas id="c1"></canvas></div></div>
  <div class="card"><div class="ct" id="c2t"></div><div class="chart"><canvas id="c2"></canvas></div></div>
</div>
<div class="grid2">
  <div class="card"><div class="ct" id="c3t"></div><div class="chart"><canvas id="c3"></canvas></div></div>
  <div class="card"><div class="ct" id="c4t"></div><div class="chart"><canvas id="c4"></canvas></div></div>
</div>
<div class="sec" id="sec-leads"></div>
<div class="rules" id="method-note" style="background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:14px 18px;margin-bottom:14px;line-height:1.7"></div>
<div class="card" style="margin-bottom:24px"><div style="max-height:480px;overflow:auto"><table id="tbl-leads"><thead></thead><tbody></tbody></table></div></div>
<div class="sec" id="sec-iss"></div>
<div id="issues"></div>
<script>
window.__T2LEAD__ = __DATA_JSON__;
</script>
<script>
const DATA = window.__T2LEAD__;
let lang = 'en';
let charts = {};

function T(k, rep) {
  const d = (DATA.i18n && DATA.i18n[lang]) || {};
  let s = d[k] || k;
  if (rep) for (const [a,b] of Object.entries(rep)) s = s.split('{' + a + '}').join(String(b));
  return s;
}
function num(v) {
  if (v === null || v === undefined) return null;
  const x = Number(v);
  return Number.isFinite(x) ? x : null;
}
function fmt(v, digits=3) {
  const x = num(v);
  if (x === null) return '—';
  if (Math.abs(x) > 1e6) return x.toExponential(2);
  return Number(x.toFixed(digits));
}
function bounds(values, pad=0.08, minSpan=0.05, hardMin=null, hardMax=null) {
  const vals = values.map(num).filter(v => v !== null);
  if (!vals.length) return [0, 1];
  let mn = Math.min(...vals), mx = Math.max(...vals);
  if (mx - mn < minSpan) {
    const m = (mx + mn) / 2;
    mn = m - minSpan / 2;
    mx = m + minSpan / 2;
  }
  let span = mx - mn;
  mn -= span * pad;
  mx += span * pad;
  if (hardMin !== null) mn = Math.max(hardMin, mn);
  if (hardMax !== null) mx = Math.min(hardMax, mx);
  return [mn, mx];
}
function clip01(x) { return Math.max(0, Math.min(1, x)); }
function fastScore(row) {
  const preset = num(row.fast_score);
  if (preset !== null) return preset;
  const parts = [
    [0.30, num(row.pred_pIC50_ens) === null ? null : clip01((num(row.pred_pIC50_ens) - 5.0) / 6.0)],
    [0.10, num(row.QED) === null ? null : clip01(num(row.QED))],
    [0.20, num(row.mpo_score) === null ? null : clip01(num(row.mpo_score))],
    [0.20, num(row.docking_score) === null ? null : clip01(((-num(row.docking_score)) - 4.0) / 6.0)],
    [0.10, num(row.sa_score) === null ? null : clip01((6.0 - num(row.sa_score)) / 5.0)],
    [0.10, num(row.admet_risk) === null ? null : clip01(1.0 - num(row.admet_risk))]
  ];
  let ws = 0, ss = 0;
  for (const [w, v] of parts) {
    if (v === null) continue;
    ws += w;
    ss += w * v;
  }
  return ws > 1e-9 ? ss / ws : null;
}
function setLang(l) {
  lang = l;
  document.getElementById('btn-en').classList.toggle('on', l === 'en');
  document.getElementById('btn-zh').classList.toggle('on', l === 'zh');
  renderAll();
}
document.getElementById('btn-en').onclick = () => setLang('en');
document.getElementById('btn-zh').onclick = () => setLang('zh');

function renderAll() {
  const D = DATA;
  document.getElementById('title').textContent = T('title_tpl', { disease: D.disease });
  document.getElementById('subtitle').textContent = T('subtitle_tpl', {
    sym: D.target_symbol, tid: D.target_chembl_id || '—', pdb: D.pdb_id, log: D.log_name || '—'
  });
  const sc = D.stage_counts || {};
  const extra = sc.s3_extra ? (' + ' + Number(sc.s3_extra).toLocaleString() + ' REINVENT4') : '';
  document.getElementById('flow').innerHTML = `
    <div class="node"><div class="n">${T('stg1')}</div><div>${T('stg1name')}</div><div class="v">${sc.s1||0}</div></div><span class="arr">→</span>
    <div class="node"><div class="n">${T('stg2')}</div><div>${T('stg2name')}</div><div class="v">${sc.s2||0}</div><div style="font-size:.68rem;color:var(--muted)">${T('stg2detail_tpl',{n:Number((D.funnel_vals&&D.funnel_vals[0])||0).toLocaleString()})}</div></div><span class="arr">→</span>
    <div class="node"><div class="n">${T('stg3')}</div><div>${T('stg3name')}</div><div class="v">${sc.s3||0}</div><div style="font-size:.68rem;color:var(--muted)">${extra}</div></div><span class="arr">→</span>
    <div class="node"><div class="n">${T('stg4')}</div><div>${T('stg4name')}</div><div class="v">${sc.s4||0}</div><div style="font-size:.68rem;color:var(--muted)">${T('stg4detail')}</div></div>`;

  const cv = D.cv_r2 != null ? Number(D.cv_r2).toFixed(3) : '—';
  document.getElementById('metrics').innerHTML = `
    <div class="card"><div class="ct">${T('metric_cv')}</div><div class="val" style="color:var(--a1)">${cv}</div><div style="color:var(--muted);font-size:.78rem;margin-top:4px">${T('metric_cv_sub')}</div></div>
    <div class="card"><div class="ct">${T('metric_dock')}</div><div class="val" style="color:var(--a2)">${D.dock_ok!=null&&D.dock_tot!=null ? (D.dock_ok+'/'+D.dock_tot) : '—'}</div><div style="color:var(--muted);font-size:.78rem;margin-top:4px">${T('metric_dock_sub_tpl',{ok:D.dock_ok||0,tot:D.dock_tot||0})}</div></div>
    <div class="card"><div class="ct">${T('metric_admet')}</div><div class="val" style="color:var(--purple)">${D.admet_flagged ?? 0}</div><div style="color:var(--muted);font-size:.78rem;margin-top:4px">${T('metric_admet_sub_tpl',{n:D.admet_flagged ?? 0})}</div></div>
    <div class="card"><div class="ct">${T('metric_appr')}</div><div class="val" style="color:var(--gold)">${D.n_approved_leads||0}/${D.n_optimized||0}</div><div style="color:var(--muted);font-size:.78rem;margin-top:4px">${T('metric_appr_sub_tpl',{n:D.n_approved_leads||0,tot:D.n_optimized||0})}</div></div>`;

  document.getElementById('sec-tg').textContent = T('sec_targets');
  document.getElementById('sec-fun').textContent = T('sec_funnel');
  document.getElementById('sec-ch').textContent = T('sec_charts');
  document.getElementById('sec-leads').textContent = T('sec_leads');
  document.getElementById('method-note').innerHTML = T('method_note');
  document.getElementById('sec-iss').textContent = T('sec_issues');
  const rulesEl = document.getElementById('chart-rules');
  rulesEl.innerHTML = T('chart_rules').replace(/green/g,'<b style="color:#00e5a0">■ green</b>')
    .replace(/cyan/g,'<b style="color:#00b8d4">■ cyan</b>')
    .replace(/amber/g,'<b style="color:#fbbf24">■ amber</b>')
    .replace(/red/g,'<b style="color:#ff6b6b">■ red</b>')
    .replace(/绿色/g,'<b style="color:#00e5a0">■ 绿色</b>')
    .replace(/青色/g,'<b style="color:#00b8d4">■ 青色</b>')
    .replace(/黄色/g,'<b style="color:#fbbf24">■ 黄色</b>')
    .replace(/红色/g,'<b style="color:#ff6b6b">■ 红色</b>') + ' ' + T('axis_zoom_note');
  document.getElementById('score-note').textContent = T('score_note');

  const th = `<tr><th>${T('th_rank')}</th><th>${T('th_symbol')}</th><th>ChEMBL</th><th>${T('th_name')}</th><th>${T('th_score')}</th></tr>`;
  document.querySelector('#tbl-targets thead').innerHTML = th;
  const tb = document.querySelector('#tbl-targets tbody');
  tb.innerHTML = '';
  (D.targets||[]).forEach((r,i) => {
    const tr = document.createElement('tr');
    tr.innerHTML = `<td>${i+1}</td><td><strong>${r.symbol||''}</strong></td><td>${r.chembl_id||''}</td><td>${(r.name||'').slice(0,80)}</td><td>${typeof r.score==='number'?r.score.toFixed(3):r.score}</td>`;
    tb.appendChild(tr);
  });

  document.getElementById('c1t').textContent = T('chart_dock_mpo');
  document.getElementById('c2t').textContent = T('chart_opt');
  document.getElementById('c3t').textContent = T('chart_pic50');
  document.getElementById('c4t').textContent = T('chart_md');
  renderIssues();
  renderLeadsTable();
  renderCharts();
}

function renderIssues() {
  const D = DATA;
  const el = document.getElementById('issues');
  const iss = D.issues || [];
  const high = iss.filter(x => x.severity >= 2);
  const low = iss.filter(x => x.severity <= 1);
  let html = '';
  function card(it) {
    const cls = it.severity >= 3 ? 'e' : it.severity === 2 ? 'w' : 'i';
    const title = lang==='zh' ? it.title_zh : it.title_en;
    const desc = lang==='zh' ? it.desc_zh : it.desc_en;
    const cnt = it.count > 1 ? (' ' + T('issues_count_tpl',{n:it.count})) : '';
    return `<div class="issue ${cls}"><h4>${title}${cnt}</h4><p>${desc}</p></div>`;
  }
  high.forEach(it => { html += card(it); });
  if (low.length) {
    html += `<details class="benign"><summary>${T('issues_expand')} (${low.length})</summary>`;
    low.forEach(it => { html += card(it); });
    html += `</details>`;
  }
  el.innerHTML = html || '<p style="color:var(--muted)">—</p>';
}

function chemblCell(r) {
  const id = r.chembl_id || '';
  const url = r.chembl_url || '';
  if (!id) return '—';
  if (url) return `<a class="chembl" target="_blank" rel="noopener" href="${url}">${id}</a>`;
  return id;
}

function renderLeadsTable() {
  const D = DATA;
  const head = `<tr><th>${T('th_type')}</th><th>${T('th_smiles')}</th><th>${T('th_pic50')}</th><th>${T('th_qed')}</th><th>${T('th_mpo')}</th><th>${T('th_dock')}</th><th>${T('th_sa')}</th><th>${T('th_md')}</th><th>${T('th_rmsd')}</th><th>${T('th_fast')}</th><th>${T('th_opt')}</th><th>${T('th_chembl')}</th></tr>`;
  document.querySelector('#tbl-leads thead').innerHTML = head;
  const tb = document.querySelector('#tbl-leads tbody');
  tb.innerHTML = '';
  const bench = D.benchmark || [];
  const leads = D.leads || [];
  if (!bench.length && !leads.length) {
    tb.innerHTML = `<tr><td colspan="12" style="color:var(--muted)">—</td></tr>`;
    return;
  }
  bench.forEach(r => {
    const tr = document.createElement('tr');
    tr.className = 'row-bench';
    const nm = r.pref_name || '';
    const smLabel = nm ? nm.slice(0,34) : ((r.canonical_smiles||'').slice(0,36) + ((r.canonical_smiles||'').length>36?'…':''));
    tr.innerHTML = `<td><span class="badge b-bench">${T('badge_benchmark')}</span></td><td class="sm" title="${(r.canonical_smiles||'').replace(/"/g,'&quot;')}" style="${nm?'color:var(--a1);font-weight:500':''}">${smLabel}</td>
      <td>${fmt(r.pred_pIC50_ens)}</td><td>${fmt(r.QED)}</td><td>${fmt(r.mpo_score)}</td><td>${fmt(r.docking_score)}</td><td>${fmt(r.sa_score)}</td>
      <td>${fmt(r.md_binding_energy)}</td><td>${fmt(r.md_rmsd_mean)}</td><td><strong>${fmt(fastScore(r))}</strong></td><td>${fmt(r.opt_score)}</td><td>${chemblCell(r)}</td>`;
    tb.appendChild(tr);
  });
  if (!bench.length) {
    const tr = document.createElement('tr');
    tr.innerHTML = `<td colspan="12" style="color:var(--muted);font-size:.8rem">${T('no_benchmark')}</td>`;
    tb.appendChild(tr);
  }
  leads.forEach(r => {
    const tr = document.createElement('tr');
    const sm = (r.canonical_smiles||'').slice(0,36) + ((r.canonical_smiles||'').length>36?'…':'');
    tr.innerHTML = `<td><span class="badge b-lead">${T('badge_lead')}</span></td><td class="sm" title="${(r.canonical_smiles||'').replace(/"/g,'&quot;')}">${sm}</td>
      <td>${fmt(r.pred_pIC50_ens)}</td><td>${fmt(r.QED)}</td><td>${fmt(r.mpo_score)}</td><td>${fmt(r.docking_score)}</td><td>${fmt(r.sa_score)}</td>
      <td>${fmt(r.md_binding_energy)}</td><td>${fmt(r.md_rmsd_mean)}</td><td><strong>${fmt(fastScore(r))}</strong></td><td>${fmt(r.opt_score)}</td><td>${chemblCell(r)}</td>`;
    tb.appendChild(tr);
  });
}

function destroyCharts() {
  Object.values(charts).forEach(c => c && c.destroy && c.destroy());
  charts = {};
}
function colorByRmsd(r, isBench) {
  if (isBench) return '#00e5a0';
  const v = num(r.md_rmsd_mean);
  if (v === null) return '#00b8d4';
  if (v > 9) return '#ff6b6b';
  if (v > 6) return '#fbbf24';
  return '#00b8d4';
}
function renderCharts() {
  destroyCharts();
  const D = DATA;
  Chart.defaults.color = '#8899b4';
  Chart.defaults.borderColor = 'rgba(30,45,74,.45)';
  Chart.defaults.font.family = 'Space Grotesk';

  const fv = D.funnel_vals || [0,0,0,0,0,0];
  const fl = [T('funnel_lib'), T('funnel_p1'), T('funnel_qed'), T('funnel_pains'), T('funnel_admet'), T('funnel_hits')];
  charts.f = new Chart(document.getElementById('cFunnel'), {
    type: 'bar',
    data: { labels: fl, datasets: [{ data: fv, backgroundColor: ['rgba(167,139,250,.55)','rgba(0,184,212,.55)','rgba(0,229,160,.55)','rgba(251,191,36,.55)','rgba(0,229,160,.65)','rgba(0,229,160,.9)'], borderRadius: 4 }] },
    options: { indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: false }, tooltip: { callbacks: { label: function(ctx) { return Number(ctx.raw).toLocaleString() + ' molecules'; } } } },
      scales: { x: { type: 'logarithmic', title: { display: true, text: lang==='zh'?'分子数量（对数刻度）':'Molecule count (log scale)' }, grid: { color: 'rgba(30,45,74,.3)' } }, y: { grid: { display: false } } }
    }
  });

  const bench = D.benchmark || [];
  const leads = D.leads || [];
  const rows = ([]).concat(bench, leads);
  const isBenchN = bench.length;
  function rlabel(i) {
    if (i < isBenchN) {
      const nm = bench[i].pref_name || bench[i].chembl_id || '';
      return 'Bench ' + (i+1) + (nm ? (' (' + nm.slice(0,18) + ')') : '');
    }
    return 'Lead ' + (i - isBenchN + 1);
  }
  const labels = rows.map((_, i) => rlabel(i));

  const dockVals = rows.map(r => num(r.docking_score));
  const mpoVals = rows.map(r => num(r.mpo_score));
  const fastVals = rows.map(r => num(fastScore(r)));
  const optVals = rows.map(r => num(r.opt_score));
  const picVals = rows.map(r => num(r.pred_pIC50_ens));
  const rmsdVals = rows.map(r => num(r.md_rmsd_mean));
  const saVals = rows.map(r => num(r.sa_score));

  const [dockMin, dockMax] = bounds(dockVals, 0.08, 0.6);
  const [mpoMin, mpoMax] = bounds(mpoVals, 0.08, 0.05, 0, 1);
  const [fastMin, fastMax] = bounds(fastVals, 0.15, 0.03, 0, 1);
  const [picMin, picMax] = bounds(picVals, 0.08, 0.4);
  const [rmsdMin, rmsdMax] = bounds(rmsdVals, 0.08, 1.0, 0, null);
  const [saMin, saMax] = bounds(saVals, 0.08, 0.5, 0, 10);

  charts.c1 = new Chart(document.getElementById('c1'), {
    type: 'scatter',
    data: { datasets: [{
      data: rows.map((r, i) => ({ x: num(r.docking_score), y: num(r.mpo_score), idx: i })),
      backgroundColor: rows.map((r, i) => colorByRmsd(r, i < isBenchN)),
      pointRadius: 9,
      pointHoverRadius: 13
    }] },
    options: {
      responsive: true, maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: function(ctx) {
              const d = ctx.raw;
              return rlabel(d.idx) + ': Dock=' + fmt(d.x) + '  MPO=' + fmt(d.y);
            }
          }
        }
      },
      scales: {
        x: { min: dockMin, max: dockMax, title: { display: true, text: 'Dock (kcal/mol)' }, grid: { color: 'rgba(30,45,74,.3)' } },
        y: { min: mpoMin, max: mpoMax, title: { display: true, text: 'MPO' }, grid: { color: 'rgba(30,45,74,.3)' } }
      }
    }
  });

  charts.c2 = new Chart(document.getElementById('c2'), {
    type: 'bar',
    data: {
      labels,
      datasets: [{
        data: fastVals,
        backgroundColor: rows.map((r, i) => colorByRmsd(r, i < isBenchN)),
        borderRadius: 6,
        barThickness: 24
      }]
    },
    options: {
      indexAxis: 'y',
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: function(ctx) {
              const i = ctx.dataIndex;
              return 'Fast=' + fmt(ctx.raw, 4) + '  |  Opt(struct)=' + fmt(optVals[i], 4);
            }
          }
        }
      },
      scales: {
        x: {
          min: fastMin,
          max: fastMax,
          title: { display: true, text: 'Fast Composite Score (zoomed)' },
          grid: { color: 'rgba(30,45,74,.3)' }
        },
        y: { grid: { display: false } }
      }
    }
  });

  charts.c3 = new Chart(document.getElementById('c3'), {
    type: 'bar',
    data: { labels, datasets: [{ data: picVals, backgroundColor: 'rgba(167,139,250,.65)', borderColor: '#a78bfa', borderWidth: 1, borderRadius: 4 }] },
    options: {
      responsive: true, maintainAspectRatio: false, plugins: { legend: { display: false } },
      scales: {
        y: { min: picMin, max: picMax, title: { display: true, text: 'pred pIC50 (zoomed)' }, grid: { color: 'rgba(30,45,74,.3)' } },
        x: { grid: { display: false } }
      }
    }
  });

  charts.c4 = new Chart(document.getElementById('c4'), {
    type: 'bar',
    data: {
      labels,
      datasets: [
        {
          label: 'MD RMSD (Å)',
          data: rmsdVals,
          backgroundColor: rows.map((r, i) => colorByRmsd(r, i < isBenchN)),
          borderRadius: 4,
          yAxisID: 'y'
        },
        {
          label: 'SA Score',
          data: saVals,
          type: 'line',
          borderColor: '#a78bfa',
          backgroundColor: 'rgba(167,139,250,.1)',
          pointRadius: 6,
          pointBackgroundColor: '#a78bfa',
          tension: .3,
          yAxisID: 'y1'
        }
      ]
    },
    options: {
      responsive: true, maintainAspectRatio: false,
      plugins: { legend: { position: 'top', labels: { boxWidth: 12, padding: 16 } } },
      scales: {
        y: { min: rmsdMin, max: rmsdMax, title: { display: true, text: 'RMSD (Å)' }, position: 'left', grid: { color: 'rgba(30,45,74,.3)' } },
        y1: { min: saMin, max: saMax, title: { display: true, text: 'SA Score' }, position: 'right', grid: { display: false } },
        x: { grid: { display: false } }
      }
    }
  });
}

renderAll();
</script>
</body>
</html>
"""
    return template.replace("__DATA_JSON__", data_json)


def _dashboard_artifact_dir(
    run_root: Path,
    layout: Dict[str, Path],
    stage4_dir: Optional[Path],
) -> Path:
    """Where to write dashboard HTML, payload JSON, and per-run score snapshot.

    Variant-driven runs pass *stage4_dir* like ``.../stage4_optimization/EGFR_L858R``;
    that differs from ``layout[STAGE4]`` (the stage root), so artifacts go under the
    variant folder and do not overwrite ``<run_root>/dashboard.html`` used by
    non-variant or legacy outputs.
    """
    base_s4 = layout[STAGE4]
    if stage4_dir is not None:
        try:
            a = Path(stage4_dir).resolve()
            b = Path(base_s4).resolve()
            if a != b:
                return a
        except OSError:
            pass
    return Path(run_root).resolve()


def _write_score_cache(
    cfg: Dict[str, Any],
    run_root: Path,
    payload: Dict[str, Any],
    *,
    snapshot_dir: Optional[Path] = None,
) -> None:
    """Persist per-run and cross-run score cache keyed by chembl_id/smiles hash."""
    rows: List[Dict[str, Any]] = []
    disease = payload.get("disease", "")
    target = payload.get("target_chembl_id", "")
    ts = payload.get("generated_at", datetime.now(timezone.utc).isoformat())

    def push(source: str, rec: Dict[str, Any]) -> None:
        smi = str(rec.get("canonical_smiles") or "")
        cid = str(rec.get("chembl_id") or "")
        key = cid if cid else hashlib.sha1(smi.encode("utf-8")).hexdigest()[:16]
        rows.append(
            {
                "generated_at": ts,
                "run_root": str(run_root),
                "disease": disease,
                "target_chembl_id": target,
                "source": source,
                "compound_key": key,
                "chembl_id": cid,
                "pref_name": rec.get("pref_name", ""),
                "canonical_smiles": smi,
                "pred_pIC50_ens": rec.get("pred_pIC50_ens"),
                "QED": rec.get("QED"),
                "mpo_score": rec.get("mpo_score"),
                "docking_score": rec.get("docking_score"),
                "md_binding_energy": rec.get("md_binding_energy"),
                "md_rmsd_mean": rec.get("md_rmsd_mean"),
                "fast_score": rec.get("fast_score"),
                "opt_score": rec.get("opt_score"),
            }
        )

    for r in payload.get("benchmark", []):
        push("benchmark", r)
    for r in payload.get("leads", []):
        push("lead", r)
    if not rows:
        return

    snap_root = Path(snapshot_dir) if snapshot_dir is not None else Path(run_root)
    run_csv = snap_root / "compound_scores_snapshot.csv"
    pd.DataFrame(rows).to_csv(run_csv, index=False)

    base_out = Path((cfg.get("pipeline", {}) or {}).get("out_dir", "./data")).expanduser().resolve()
    cache_dir = base_out / "shared" / "report_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_csv = cache_dir / "compound_scores.csv"

    if cache_csv.is_file():
        old = pd.read_csv(cache_csv)
        merged = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    else:
        merged = pd.DataFrame(rows)
    merged = merged.drop_duplicates(
        subset=["run_root", "source", "compound_key"],
        keep="last",
    )
    merged.to_csv(cache_csv, index=False)


def write_dashboard(
    cfg: Dict[str, Any],
    run_root: Path,
    layout: Dict[str, Path],
    stage4_dir: Optional[Path] = None,
) -> Path:
    payload = build_payload(cfg, run_root, layout, stage4_dir)
    html = render_html(payload)
    art_dir = _dashboard_artifact_dir(run_root, layout, stage4_dir)
    art_dir.mkdir(parents=True, exist_ok=True)
    out = art_dir / "dashboard.html"
    out.write_text(html, encoding="utf-8")
    (art_dir / "report_payload.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    _write_score_cache(cfg, Path(run_root), payload, snapshot_dir=art_dir)
    return out
