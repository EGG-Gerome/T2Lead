"""Stage-scoped output directory layout under each run folder.
各次运行目录下、按阶段划分的输出路径布局（与 pipeline 编排器约定一致）。
"""

from __future__ import annotations

import datetime
import re
from pathlib import Path
from typing import Any, Dict, Optional


def disease_slug(disease: str) -> str:
    """Filesystem-safe directory name for a disease label.
    将疾病标签转为可用于文件系统的目录名（小写、空白与非法字符规范化）。
    """
    slug = disease.strip().lower()
    slug = re.sub(r"[^\w\s-]", "", slug)
    slug = re.sub(r"[\s]+", "_", slug)
    return slug or "default"


def _token_slug(text: str, fallback: str) -> str:
    """Normalize free-text token for directory names."""
    tok = str(text or "").strip().lower()
    tok = re.sub(r"[^\w\s-]", "", tok)
    tok = re.sub(r"[\s]+", "_", tok)
    return tok or fallback


def _strip_known_suffixes(name: str) -> str:
    """Drop common FASTQ/VCF extensions from a filename."""
    low = name.lower()
    for suf in (".fastq.gz", ".fq.gz", ".vcf.gz", ".fastq", ".fq", ".vcf"):
        if low.endswith(suf):
            return name[: -len(suf)]
    return Path(name).stem


def _variant_sample_id(cfg: Dict[str, Any]) -> str:
    """Resolve patient/sample ID for variant-isolated run folders."""
    va = (cfg.get("variant_analysis", {}) or {})
    explicit = str(va.get("sample_id", "") or "").strip()
    if explicit:
        return _token_slug(explicit, "sample")

    vcf_path = str(va.get("vcf_path", "") or "").strip()
    if vcf_path:
        return _token_slug(_strip_known_suffixes(Path(vcf_path).name), "sample")

    tumor_r1 = str(va.get("tumor_fastq_r1", "") or "").strip()
    if tumor_r1:
        return _token_slug(_strip_known_suffixes(Path(tumor_r1).name), "sample")
    return "sample"


def _variant_run_id(cfg: Dict[str, Any]) -> str:
    """Stable run id for this process (explicit override or timestamp)."""
    layout = (cfg.get("pipeline", {}) or {}).get("output_layout", {}) or {}
    explicit = str(layout.get("variant_run_id", "") or "").strip()
    if explicit:
        return _token_slug(explicit, "run")

    runtime = cfg.setdefault("_runtime", {})
    cached = str(runtime.get("variant_run_id", "") or "").strip()
    if cached:
        return cached
    rid = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    runtime["variant_run_id"] = rid
    return rid


def _use_variant_isolated_runs(cfg: Dict[str, Any]) -> bool:
    """Whether variant path should create unique run folders (default: true)."""
    va = (cfg.get("variant_analysis", {}) or {})
    if not bool(va.get("enabled", False)):
        return False
    layout = (cfg.get("pipeline", {}) or {}).get("output_layout", {}) or {}
    return bool(layout.get("variant_isolated_runs", True))


def run_root_for_config(cfg: Dict[str, Any]) -> Path:
    """Directory root for this run.

    Standard path:
      ``out_dir`` or ``out_dir/<disease_slug>``
    Variant path (when enabled + isolated runs):
      ``.../<disease_slug>/variant_runs/<sample_id>/<run_id>``

    本次运行根目录：
    - 标准路径：``out_dir`` 或 ``out_dir/<disease_slug>``
    - 变异路径（启用隔离）：``.../<disease_slug>/variant_runs/<sample_id>/<run_id>``
    """
    from drugpipe.config import get_out_dir

    base = get_out_dir(cfg)
    d = (cfg.get("target_discovery", {}) or {}).get("disease", "") or ""
    d = str(d).strip()
    parent = (base / disease_slug(d)) if d else base

    if _use_variant_isolated_runs(cfg):
        sample_id = _variant_sample_id(cfg)
        run_id = _variant_run_id(cfg)
        r = parent / "variant_runs" / sample_id / run_id
        r.mkdir(parents=True, exist_ok=True)
        return r

    parent.mkdir(parents=True, exist_ok=True)
    return parent


# Subdirectory names (per-disease or default run root).
# 子目录名（按疾病分目录或默认运行根下共用）。
STAGE1 = "stage1_targets"
STAGE2 = "stage2_hits"
STAGE3 = "stage3_leads"
STAGE4 = "stage4_optimization"


def use_stage_subdirs(cfg: Dict[str, Any]) -> bool:
    """Return True when outputs should use stage/* subfolders.
    为 True 时各阶段写入独立子目录；为 False 时阶段输出合并到同一目录。
    """
    layout = cfg.get("pipeline", {}).get("output_layout", {})
    return bool(layout.get("use_stage_subdirs", True))


def stage_paths(run_root: Path, cfg: Dict[str, Any]) -> Dict[str, Path]:
    """Return output paths for each stage under *run_root*.

    When ``use_stage_subdirs`` is False, every stage maps to *run_root*.
    返回 *run_root* 下各阶段输出路径；关闭分目录时所有阶段键均指向 *run_root*。
    """
    run_root = Path(run_root)
    if not use_stage_subdirs(cfg):
        p = run_root
        return {
            "run_root": run_root,
            STAGE1: p,
            STAGE2: p,
            STAGE3: p,
            STAGE4: p,
            "logs": run_root / "logs",
        }
    d1 = run_root / STAGE1
    d2 = run_root / STAGE2
    d3 = run_root / STAGE3
    d4 = run_root / STAGE4
    for p in (d1, d2, d3, d4):
        p.mkdir(parents=True, exist_ok=True)
    log_dir = run_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return {
        "run_root": run_root,
        STAGE1: d1,
        STAGE2: d2,
        STAGE3: d3,
        STAGE4: d4,
        "logs": log_dir,
    }


def crawl_and_stage2_dir(cfg: Dict[str, Any]) -> Path:
    """Per-run Stage 2 directory (hits, models, screening outputs).

    When ``target_to_hit.shared_library_dir`` is set, ChEMBL crawl files and
    ``fp_cache/`` live under :func:`shared_chembl_library_dir` instead; this
    path still holds disease-local artifacts (``dataset_*``, ``model_cache/``,
    ``scored_candidates.csv``, …).

    单次运行的阶段二目录。若配置了 ``shared_library_dir``，ChEMBL 与指纹缓存
    在共享目录；此处仍写疾病/靶点专属产物。
    """
    return stage_paths(run_root_for_config(cfg), cfg)[STAGE2]


def shared_chembl_library_dir(cfg: Dict[str, Any]) -> Optional[Path]:
    """Directory for shared ChEMBL crawl + Morgan ``fp_cache`` (optional).

    If ``target_to_hit.shared_library_dir`` is empty, returns ``None`` (legacy:
    crawl and fingerprints use ``stage2_hits``).

    Relative paths are resolved under :func:`drugpipe.config.get_out_dir`
    (pipeline base ``out_dir``, not the per-disease run root).

    跨疾病共享的 ChEMBL 爬取与指纹缓存目录；未配置则返回 ``None``（沿用
    ``stage2_hits``）。相对路径相对于 ``pipeline.out_dir``。
    """
    from drugpipe.config import get_out_dir

    raw = (cfg.get("target_to_hit", {}) or {}).get("shared_library_dir", "") or ""
    raw = str(raw).strip()
    if not raw:
        return None
    p = Path(raw)
    if not p.is_absolute():
        p = (get_out_dir(cfg) / raw).resolve()
    else:
        p = p.resolve()
    p.mkdir(parents=True, exist_ok=True)
    return p


def chembl_library_root(cfg: Dict[str, Any], stage2_local: Path) -> Path:
    """Resolve crawl + fp_cache root: shared dir if configured, else *stage2_local*."""
    shared = shared_chembl_library_dir(cfg)
    return shared if shared is not None else Path(stage2_local)


def resolve_variant_stage4_out(run_root: Path, cfg: Dict[str, Any], gene: str, mutation: str) -> Path:
    """Per-variant Stage 4 directory (e.g. stage4_optimization/EGFR_L858R/).
    按突变划分的阶段四输出目录，例如 ``stage4_optimization/EGFR_L858R/``。
    """
    base = stage_paths(run_root, cfg)[STAGE4]
    safe_mut = mutation.replace(":", "_").replace("/", "_")
    sub = base / f"{gene}_{safe_mut}"
    sub.mkdir(parents=True, exist_ok=True)
    return sub
