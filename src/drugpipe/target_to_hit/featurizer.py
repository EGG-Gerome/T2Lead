"""Morgan fingerprint featurizer."""
# Morgan 指纹特征化器。

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from drugpipe.utils.chem import batch_fps


class MorganFeaturizer:
    """Compute Morgan (ECFP-like) fingerprint bit vectors from SMILES."""
    # 从 SMILES 计算 Morgan（ECFP 类）指纹位向量。

    def __init__(self, cfg: Dict[str, Any]):
        feat_cfg = cfg.get("target_to_hit", {}).get("featurizer", {})
        self.radius = int(feat_cfg.get("fp_radius", 2))
        self.nbits = int(feat_cfg.get("fp_nbits", 2048))

    def transform(self, smiles_list: List[str]) -> np.ndarray:
        """Return shape ``(N, nbits)`` float32 matrix."""
        # 返回形状 (N, nbits) 的 float32 矩阵。
        return batch_fps(smiles_list, radius=self.radius, nbits=self.nbits)
