"""Morgan fingerprint featurizer with on-disk caching."""
# Morgan 指纹特征化器（支持磁盘缓存）。

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from drugpipe.utils.chem import batch_fps

logger = logging.getLogger(__name__)


def _fp_cache_key(smiles_list: List[str], radius: int, nbits: int) -> str:
    """Deterministic hash of SMILES list + FP parameters."""
    h = hashlib.sha256()
    h.update(f"r={radius},n={nbits},len={len(smiles_list)}\n".encode())
    for smi in smiles_list:
        h.update(smi.encode())
        h.update(b"\n")
    return h.hexdigest()[:16]


class MorganFeaturizer:
    """Compute Morgan (ECFP-like) fingerprint bit vectors from SMILES.

    When *cache_dir* is set, fingerprint matrices are saved to ``.npy``
    files keyed by a hash of (SMILES list + radius + nbits).  Subsequent
    calls with the same inputs load the cached file in < 1 second instead
    of re-computing (~11 min for 2.8M molecules).
    """

    def __init__(self, cfg: Dict[str, Any], cache_dir: Optional[Path] = None):
        feat_cfg = cfg.get("target_to_hit", {}).get("featurizer", {})
        self.radius = int(feat_cfg.get("fp_radius", 2))
        self.nbits = int(feat_cfg.get("fp_nbits", 2048))
        self.cache_dir = cache_dir

    def transform(self, smiles_list: List[str]) -> np.ndarray:
        """Return shape ``(N, nbits)`` float32 matrix, with caching."""
        if self.cache_dir is not None:
            cached = self._load_cache(smiles_list)
            if cached is not None:
                return cached

        X = batch_fps(smiles_list, radius=self.radius, nbits=self.nbits)

        if self.cache_dir is not None:
            self._save_cache(smiles_list, X)

        return X

    # ------------------------------------------------------------------
    def _cache_path(self, smiles_list: List[str]) -> Path:
        key = _fp_cache_key(smiles_list, self.radius, self.nbits)
        return self.cache_dir / f"morgan_r{self.radius}_b{self.nbits}_{key}.npy"

    def _load_cache(self, smiles_list: List[str]) -> Optional[np.ndarray]:
        path = self._cache_path(smiles_list)
        if not path.exists():
            return None
        X = np.load(path, mmap_mode="r")
        if X.shape == (len(smiles_list), self.nbits):
            logger.info("Loaded cached fingerprints: %s (%d × %d)", path.name, *X.shape)
            return np.array(X, dtype=np.float32)
        logger.warning("Cache shape mismatch, recomputing.")
        return None

    def _save_cache(self, smiles_list: List[str], X: np.ndarray) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self._cache_path(smiles_list)
        np.save(path, X)
        size_mb = path.stat().st_size / (1024 * 1024)
        logger.info("Saved fingerprint cache: %s (%.1f MB)", path.name, size_mb)
