"""Configuration loader: YAML file + environment variable overrides."""
# EN: Configuration loader: YAML file + environment variable overrides.
# 中文：说明模块职责、上下游关系与维护注意事项。

# 配置加载：YAML 文件 + 环境变量覆盖。

from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict

import yaml
from dotenv import load_dotenv

_DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "configs" / "default_config.yaml"


def _deep_merge(base: dict, override: dict) -> dict:
    """EN: Recursively merge `override` into `base` in place.
    中文：将 `override` 递归合并到 `base`，会原地修改 `base`。
    """
    for k, v in override.items():
        if k in base and isinstance(base[k], dict) and isinstance(v, dict):
            _deep_merge(base[k], v)
        else:
            base[k] = v
    return base


def _apply_env_overrides(cfg: dict, prefix: str = "DP") -> dict:
    """EN: Map environment variables with a prefix into nested config keys.
    中文：将指定前缀的环境变量映射为嵌套配置字段。

    EN: `DP_TARGET_TO_HIT__CHEMBL__PAGE_LIMIT=500`
    maps to `cfg["target_to_hit"]["chembl"]["page_limit"] = 500`.
    中文：例如 `DP_TARGET_TO_HIT__CHEMBL__PAGE_LIMIT=500` 会映射到
    `cfg["target_to_hit"]["chembl"]["page_limit"] = 500`。
    """
    # 扫描以 prefix 开头的环境变量并注入配置。约定：DP_TARGET_TO_HIT__CHEMBL__PAGE_LIMIT=500 映射到 cfg["target_to_hit"]["chembl"]["page_limit"]。
    for key, val in os.environ.items():
        if not key.startswith(prefix + "_"):
            continue
        parts = key[len(prefix) + 1 :].lower().split("__")
        node = cfg
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        leaf = parts[-1]
        # Best-effort type cast / 尽力类型转换
        node[leaf] = _auto_cast(val)
    return cfg


def _auto_cast(val: str) -> Any:
    """EN: Cast env string to bool/int/float/str when possible.
    中文：尽可能将环境变量字符串转换为 bool/int/float/str。
    """
    if val.lower() in ("true", "yes"):
        return True
    if val.lower() in ("false", "no"):
        return False
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return val


def load_config(
    config_path: str | Path | None = None,
    overrides: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """EN: Load pipeline configuration from default/user/env/overrides.
    中文：按默认配置、用户配置、环境变量和运行时覆盖项加载最终配置。

    EN: Priority (highest -> lowest):
      1) `overrides` dict
      2) environment variables (`DP_*`)
      3) user YAML (`config_path`)
      4) built-in `configs/default_config.yaml`
    中文：优先级（高 -> 低）：
      1）`overrides` 字典
      2）环境变量（`DP_*`）
      3）用户 YAML（`config_path`）
      4）内置 `configs/default_config.yaml`
    """
    load_dotenv()

    with open(_DEFAULT_CONFIG, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if config_path is not None:
        with open(config_path, "r", encoding="utf-8") as f:
            user_cfg = yaml.safe_load(f) or {}
        cfg = _deep_merge(cfg, user_cfg)

    _apply_env_overrides(cfg)

    if overrides:
        cfg = _deep_merge(cfg, overrides)

    return cfg


def get_out_dir(cfg: Dict[str, Any]) -> Path:
    """EN: Return output directory and create it when missing.
    中文：返回输出目录，并在不存在时自动创建。
    """
    out = Path(cfg["pipeline"]["out_dir"])
    out.mkdir(parents=True, exist_ok=True)
    return out
