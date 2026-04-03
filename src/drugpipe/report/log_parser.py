"""Parse pipeline full logs into grouped, severity-ranked issues for the HTML report."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


class Severity(IntEnum):
    ERROR = 3
    WARNING = 2
    INFO = 1


@dataclass
class LogIssue:
    severity: Severity
    key: str
    title_en: str
    desc_en: str
    title_zh: str
    desc_zh: str
    count: int = 1
    sample_lines: List[str] = field(default_factory=list)


_IMPORTANT_PATTERNS: List[Tuple[str, re.Pattern, Severity, str, str, str, str]] = [
    (
        "nan_energy",
        re.compile(r"ΔG=nan|complex=.*e\+1[5-9]|complex=.*e\+[2-9]\d", re.I),
        Severity.WARNING,
        "MM-GBSA energy overflow / NaN",
        "At least one molecule produced NaN or extreme complex energy. The affected molecule's MD fields are empty; composite scoring fell back to other metrics.",
        "MM-GBSA 能量溢出 / NaN",
        "个别分子出现 NaN 或极大复合物能量，MD 字段为空，综合评分已回退到其他指标。",
    ),
    (
        "positive_dg",
        re.compile(r"MM-GBSA:.*ΔG=\s*[\d]{3,}[\d.]*\s*±"),
        Severity.INFO,
        "Some MM-GBSA ΔG values near zero or slightly positive",
        "Non-covalent MM-GBSA binding energies near zero are expected for covalent inhibitors evaluated without covalent bond terms. See methodology note.",
        "部分 MM-GBSA ΔG 接近零或略为正值",
        "对共价抑制剂进行非共价评估时，MM-GBSA 结合能接近零属预期结果。详见方法论说明。",
    ),
    (
        "high_rmsd",
        re.compile(r"traj_RMSD=\s*1[0-9]\.|traj_RMSD=\s*[89]\."),
        Severity.WARNING,
        "High MD trajectory RMSD (8–15+ Å)",
        "Some ligands show large displacement during MD — binding mode may be unstable for those molecules.",
        "部分 MD 轨迹 RMSD 偏高（8–15+ Å）",
        "部分配体在 MD 中位移较大，结合模式可能不稳定。",
    ),
    (
        "admet_drop",
        re.compile(r"ADMET hard filter: dropping\s+(\d+)\s+compounds"),
        Severity.WARNING,
        "ADMET hard filter removed compounds",
        "",
        "ADMET 硬过滤剔除了化合物",
        "",
    ),
]

_BENIGN_PATTERNS: List[Tuple[str, re.Pattern, str, str, str, str]] = [
    (
        "unl_lig_residue",
        re.compile(r"Did not recognize residue (UNL|LIG)"),
        "OpenMM UNL/LIG residue notice",
        "OpenMM first tries a generic residue name then generates a GAFF-2.11 template — expected during ligand parameterisation. No action needed.",
        "OpenMM UNL/LIG 残基提示",
        "OpenMM 先尝试通用残基名再生成 GAFF 模板，配体参数化时正常出现，无需处理。",
    ),
    (
        "pint_redefine",
        re.compile(r"pint\.util.*Redefining"),
        "Pint unit redefinition warning",
        "Harmless library warning from pint / OpenMM unit handling.",
        "Pint 单位重定义警告",
        "来自 pint 与 OpenMM 单位系统的无害警告。",
    ),
    (
        "torchvision_warn",
        re.compile(r"Failed to load image Python extension|torchvision"),
        "Torchvision image extension warning",
        "Torchvision image IO extension not loaded — irrelevant to the drug discovery pipeline.",
        "Torchvision 图像扩展警告",
        "Torchvision 图像 IO 扩展未加载，与药物发现流程无关。",
    ),
]


def parse_full_log(log_path: Path) -> List[LogIssue]:
    """Read *log_path* and return issues sorted by severity (error first)."""
    if not log_path.is_file():
        return []
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []

    buckets: Dict[str, LogIssue] = {}

    for line in text.splitlines():
        is_error = "[ERROR]" in line
        is_warning = "[WARNING]" in line

        if is_error:
            fk = "error_generic"
            if fk not in buckets:
                buckets[fk] = LogIssue(
                    severity=Severity.ERROR, key=fk,
                    title_en="Pipeline error", desc_en=line.strip()[:500],
                    title_zh="流水线错误", desc_zh=line.strip()[:500],
                )
            else:
                buckets[fk].count += 1
            continue

        for key, pat, sev, ten, den, tzh, dzh in _IMPORTANT_PATTERNS:
            if pat.search(line):
                if key not in buckets:
                    desc_en = den
                    desc_zh = dzh
                    if key == "admet_drop":
                        m = re.search(r"dropping\s+(\d+)\s+compounds\s*\(([^)]+)\)", line)
                        if m:
                            desc_en = f"{m.group(1)} compounds dropped ({m.group(2)}) before MD. See log for details."
                            desc_zh = f"{m.group(1)} 个化合物因 {m.group(2)} 在 MD 前被剔除。详见日志。"
                    buckets[key] = LogIssue(
                        severity=sev, key=key,
                        title_en=ten, desc_en=desc_en,
                        title_zh=tzh, desc_zh=desc_zh,
                        sample_lines=[line.strip()[:300]],
                    )
                else:
                    buckets[key].count += 1
                    if len(buckets[key].sample_lines) < 3:
                        buckets[key].sample_lines.append(line.strip()[:300])
                break

        if is_warning:
            matched_benign = False
            for key, pat, ten, den, tzh, dzh in _BENIGN_PATTERNS:
                if pat.search(line):
                    matched_benign = True
                    if key not in buckets:
                        buckets[key] = LogIssue(
                            severity=Severity.INFO, key=key,
                            title_en=ten, desc_en=den,
                            title_zh=tzh, desc_zh=dzh,
                            sample_lines=[line.strip()[:300]],
                        )
                    else:
                        buckets[key].count += 1
                    break
            if not matched_benign and "warning_generic" not in [b.key for b in buckets.values()]:
                pass

    issues = list(buckets.values())
    issues.sort(key=lambda x: (-int(x.severity), x.key))
    return issues


def issues_to_json_serializable(issues: List[LogIssue]) -> List[Dict[str, Any]]:
    """Serialize log issues for embedding in HTML JSON."""
    out: List[Dict[str, Any]] = []
    for it in issues:
        out.append(
            {
                "severity": int(it.severity),
                "severity_name": it.severity.name,
                "title_en": it.title_en,
                "desc_en": it.desc_en,
                "title_zh": it.title_zh,
                "desc_zh": it.desc_zh,
                "count": it.count,
                "collapsible": it.severity == Severity.INFO,
            }
        )
    return out


def find_latest_full_log(logs_dir: Path) -> Optional[Path]:
    """Pick the best ``*_full.log`` under *logs_dir*.

    Prefers logs that contain "PIPELINE COMPLETE" (successful runs).
    Falls back to the most recently modified log.
    """
    if not logs_dir.is_dir():
        return None
    all_logs = sorted(logs_dir.glob("*_full.log"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not all_logs:
        return None
    for candidate in all_logs:
        try:
            tail = candidate.read_bytes()[-2000:]
            if b"PIPELINE COMPLETE" in tail:
                return candidate
        except OSError:
            continue
    return all_logs[0]
