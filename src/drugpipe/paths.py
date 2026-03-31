"""Stage-scoped output directory layout under each run folder.
各次运行目录下、按阶段划分的输出路径布局（与 pipeline 编排器约定一致）。
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict


def disease_slug(disease: str) -> str:
    """Filesystem-safe directory name for a disease label.
    将疾病标签转为可用于文件系统的目录名（小写、空白与非法字符规范化）。
    """
    slug = disease.strip().lower()
    slug = re.sub(r"[^\w\s-]", "", slug)
    slug = re.sub(r"[\s]+", "_", slug)
    return slug or "default"


def run_root_for_config(cfg: Dict[str, Any]) -> Path:
    """Directory root for this run (``out_dir`` or ``out_dir/<disease_slug>``).
    本次运行的根目录：无疾病名时为 ``out_dir``，有疾病名时为 ``out_dir/<disease_slug>``（必要时创建）。
    """
    from drugpipe.config import get_out_dir

    base = get_out_dir(cfg)
    d = (cfg.get("target_discovery", {}) or {}).get("disease", "") or ""
    d = str(d).strip()
    if d:
        r = base / disease_slug(d)
        r.mkdir(parents=True, exist_ok=True)
        return r
    return base


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
    """ChEMBL crawl + Stage 2 artifacts directory.
    ChEMBL 爬取与阶段二产物所在目录（虚拟筛选、模型缓存等）。
    """
    return stage_paths(run_root_for_config(cfg), cfg)[STAGE2]


def resolve_variant_stage4_out(run_root: Path, cfg: Dict[str, Any], gene: str, mutation: str) -> Path:
    """Per-variant Stage 4 directory (e.g. stage4_optimization/EGFR_L858R/).
    按突变划分的阶段四输出目录，例如 ``stage4_optimization/EGFR_L858R/``。
    """
    base = stage_paths(run_root, cfg)[STAGE4]
    safe_mut = mutation.replace(":", "_").replace("/", "_")
    sub = base / f"{gene}_{safe_mut}"
    sub.mkdir(parents=True, exist_ok=True)
    return sub
