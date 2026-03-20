"""I/O helpers: CSV append, state persistence, directory management."""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 输入输出辅助：CSV 追加、状态持久化、目录管理。

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict

import pandas as pd


# EN: ensure_dir core behavior and intent.
# 中文：ensure_dir 的核心行为与设计意图。
def ensure_dir(path: str | Path) -> Path:
    """Create directory if not exists; return Path. / 若不存在则创建目录并返回 Path。"""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------------------
# Checkpoint / state persistence
# 断点 / 状态持久化
# ---------------------------------------------------------------------------

# EN: load_state core behavior and intent.
# 中文：load_state 的核心行为与设计意图。
def load_state(path: str | Path) -> Dict[str, Any]:
    """Load JSON state from file; return {} if missing. / 从文件加载 JSON 状态，不存在则返回 {}。"""
    path = Path(path)
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


# EN: save_state core behavior and intent.
# 中文：save_state 的核心行为与设计意图。
def save_state(state: Dict[str, Any], path: str | Path) -> None:
    """Save state to JSON (atomic write). / 将状态保存为 JSON（原子写入）。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# CSV helpers
# CSV 辅助函数
# ---------------------------------------------------------------------------

# EN: append_csv core behavior and intent.
# 中文：append_csv 的核心行为与设计意图。
def append_csv(path: str | Path, df: pd.DataFrame) -> None:
    """Append *df* to a CSV file, writing a header only if the file is new."""
    # 将 df 追加到 CSV，仅当文件为新文件时写入表头。
    path = Path(path)
    header = not path.exists()
    df.to_csv(path, mode="a", header=header, index=False)


# EN: read_csv_safe core behavior and intent.
# 中文：read_csv_safe 的核心行为与设计意图。
def read_csv_safe(path: str | Path) -> pd.DataFrame:
    """Read CSV or raise FileNotFoundError. / 读取 CSV，不存在则抛出 FileNotFoundError。"""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    return pd.read_csv(path)
