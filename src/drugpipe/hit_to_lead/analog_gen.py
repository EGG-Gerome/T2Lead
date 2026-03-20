"""Analog generation using CReM (chemical space exploration via fragment mutation)."""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 使用 CReM（通过片段突变的化学空间探索）生成类似物。

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import pandas as pd

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)

# CReM is optional / CReM 为可选依赖
try:
    from crem.crem import mutate_mol, grow_mol
    _CREM_OK = True
except ImportError:
    _CREM_OK = False


# EN: AnalogGenerator core behavior and intent.
# 中文：AnalogGenerator 的核心行为与设计意图。
class AnalogGenerator:
    """
    Generate molecular analogs for each hit compound using CReM
    (Context-aware Replacement of Environment on Molecules).

    If CReM or its fragment database is not available, this step is
    gracefully skipped and only the original hits pass through.
    """
    # 用 CReM 为每个 hit 生成分子类似物；无 CReM 或片段库时跳过，仅保留原 hit。

    # EN: __init__ core behavior and intent.
    # 中文：__init__ 的核心行为与设计意图。
    def __init__(self, cfg: Dict[str, Any]):
        h2l = cfg.get("hit_to_lead", {})
        ag = h2l.get("analog_gen", {})
        self.enabled = bool(ag.get("enabled", True))
        self.db_path = ag.get("crem_db_path", "") or ""
        self.max_per_hit = int(ag.get("max_analogs_per_hit", 50))
        self.radius = int(ag.get("crem_radius", 3))
        self.min_size = int(ag.get("crem_min_size", 0))
        self.max_size = int(ag.get("crem_max_size", 10))

    # ------------------------------------------------------------------
    # EN: generate core behavior and intent.
    # 中文：generate 的核心行为与设计意图。
    def generate(self, df_hits: pd.DataFrame) -> pd.DataFrame:
        """
        For each hit, attempt to generate structural analogs.

        Returns a DataFrame with columns:
          ``canonical_smiles``, ``source_smiles``, ``origin``
        where ``origin`` is ``"hit"`` or ``"crem_analog"``.
        """
        # 返回列：canonical_smiles、source_smiles、origin（"hit" 或 "crem_analog"）。
        rows: List[Dict[str, Any]] = []

        for _, row in df_hits.iterrows():
            smi = row["canonical_smiles"]
            rows.append({
                "canonical_smiles": smi,
                "source_smiles": smi,
                "origin": "hit",
                **{k: row[k] for k in row.index if k not in ("canonical_smiles",)},
            })

        if not self.enabled:
            logger.info("Analog generation disabled in config.")
            return pd.DataFrame(rows)

        if not _CREM_OK:
            logger.warning(
                "CReM not installed (pip install crem). Skipping analog generation."
            )
            return pd.DataFrame(rows)

        if not self.db_path:
            logger.warning(
                "crem_db_path not set — CReM needs a fragment DB. "
                "Download one from https://github.com/DrrDom/crem#databases . "
                "Skipping analog generation."
            )
            return pd.DataFrame(rows)

        logger.info(
            "Generating analogs for %d hits (max %d/hit, radius=%d) ...",
            len(df_hits), self.max_per_hit, self.radius,
        )

        for _, row in df_hits.iterrows():
            smi = row["canonical_smiles"]
            analogs = self._mutate(smi)
            for asmi in analogs:
                rows.append({
                    "canonical_smiles": asmi,
                    "source_smiles": smi,
                    "origin": "crem_analog",
                })

        df = pd.DataFrame(rows).drop_duplicates(subset=["canonical_smiles"])
        n_new = (df["origin"] == "crem_analog").sum()
        logger.info("Analog generation complete: %d new analogs + %d original hits", n_new, len(df_hits))
        return df

    # ------------------------------------------------------------------
    # EN: _mutate core behavior and intent.
    # 中文：_mutate 的核心行为与设计意图。
    def _mutate(self, smiles: str) -> List[str]:
        mol = safe_mol(smiles)
        if mol is None:
            return []
        try:
            results = list(mutate_mol(
                mol,
                db_name=self.db_path,
                radius=self.radius,
                min_size=self.min_size,
                max_size=self.max_size,
                return_mol=False,
            ))
            seen = {smiles}
            out: List[str] = []
            for smi in results:
                if smi not in seen and safe_mol(smi) is not None:
                    seen.add(smi)
                    out.append(smi)
                    if len(out) >= self.max_per_hit:
                        break
            return out
        except Exception as exc:
            logger.debug("CReM mutate failed for %s: %s", smiles, exc)
            return []
