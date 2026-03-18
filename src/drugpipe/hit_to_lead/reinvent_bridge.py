"""Optional REINVENT4 integration for RL-based molecular optimization.

REINVENT4 (https://github.com/MolecularAI/REINVENT4) is a powerful generative
model that uses reinforcement learning to optimize molecules toward desired
multi-objective profiles.  This module generates a REINVENT4 configuration,
launches it as a subprocess, and collects the output SMILES.

Requirements:
  - REINVENT4 installed (``pip install reinvent``)
  - A prior model file (e.g. ``pubchem_ecfp4_…_reinvent4_dict_voc.prior``)
  - Paths set in ``configs/default_config.yaml`` → ``hit_to_lead.reinvent4``
"""
# 可选 REINVENT4 集成：基于强化学习的分子优化；生成配置、以子进程启动并收集输出 SMILES。

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)


_TOML_TEMPLATE = """\
[parameters]
run_type = "reinforcement_learning"
model_file = "{prior_path}"
output_csv = "{output_csv}"
batch_size = {batch_size}
use_checkpoint = false

[learning_strategy]
type = "dap"
sigma = 128
rate = 0.0001

[stage.0]
max_score = 1.0
min_steps = 10
max_steps = {n_steps}
chkpt_file = "{chkpt_file}"
termination = "plateau"

[stage.0.scoring]
type = "geometric_mean"

[stage.0.scoring.component.0]
custom_alerts.endpoint = [
    "[*;$([F,Cl,Br,I])]",
    "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
]
custom_alerts.weight = 0.5

[stage.0.scoring.component.1]
QED.endpoint = ["QED"]
QED.weight = 1.0

[stage.0.scoring.component.2]
TanimotoSimilarity.endpoint = ["{reference_smiles}"]
TanimotoSimilarity.weight = 0.6
TanimotoSimilarity.params.radius = 2
TanimotoSimilarity.params.use_counts = false

{predictive_model_block}

[diversity_filter]
type = "IdenticalMurckoScaffold"
bucket_size = 25
minscore = 0.4
minsimilarity = 0.4
"""


_PREDICTIVE_MODEL_BLOCK = """\
[stage.0.scoring.component.3]
ExternalProcess.endpoint = ["{scoring_script}"]
ExternalProcess.weight = 0.8
"""


class Reinvent4Bridge:
    """Generate molecules via REINVENT4 RL optimization."""
    # 通过 REINVENT4 强化学习优化生成分子。

    def __init__(self, cfg: Dict[str, Any]):
        r4 = cfg.get("hit_to_lead", {}).get("reinvent4", {})
        self.enabled = bool(r4.get("enabled", False))
        self.reinvent_path = r4.get("reinvent_path", "") or ""
        self.prior_path = r4.get("prior_path", "") or ""
        self.batch_size = int(r4.get("batch_size", 128))
        self.n_steps = int(r4.get("n_steps", 100))

    # ------------------------------------------------------------------
    def run(
        self,
        df_hits: pd.DataFrame,
        out_dir: Path,
        model_predict_fn=None,
        featurizer_fn=None,
    ) -> pd.DataFrame:
        """
        Run REINVENT4 optimization seeded with top hit SMILES.

        If *model_predict_fn* and *featurizer_fn* are provided, a scoring
        script is generated so REINVENT4 can use the pIC50 predictor as an
        RL reward component.

        Returns a DataFrame of generated SMILES (``canonical_smiles``,
        ``origin="reinvent4"``).  Returns empty DataFrame on skip/failure.
        """
        if not self.enabled:
            logger.info("REINVENT4 bridge disabled.")
            return pd.DataFrame()

        if not self.reinvent_path or not Path(self.reinvent_path).exists():
            logger.warning("reinvent_path not found: '%s'. Skipping.", self.reinvent_path)
            return pd.DataFrame()

        if not self.prior_path or not Path(self.prior_path).exists():
            logger.warning("prior_path not found: '%s'. Skipping.", self.prior_path)
            return pd.DataFrame()

        ref_smi = df_hits["canonical_smiles"].iloc[0] if len(df_hits) > 0 else "c1ccccc1"

        pred_block = ""
        if model_predict_fn is not None and featurizer_fn is not None:
            scoring_script = self._write_scoring_script(out_dir, model_predict_fn, featurizer_fn)
            if scoring_script:
                pred_block = _PREDICTIVE_MODEL_BLOCK.format(scoring_script=str(scoring_script))

        output_csv = out_dir / "reinvent4_output.csv"
        chkpt_file = out_dir / "reinvent4_agent.ckpt"
        toml_content = _TOML_TEMPLATE.format(
            prior_path=self.prior_path,
            output_csv=str(output_csv),
            batch_size=self.batch_size,
            n_steps=self.n_steps,
            chkpt_file=str(chkpt_file),
            reference_smiles=ref_smi,
            predictive_model_block=pred_block,
        )

        toml_path = out_dir / "reinvent4_config.toml"
        toml_path.write_text(toml_content, encoding="utf-8")

        logger.info("Launching REINVENT4 (steps=%d, batch=%d) ...", self.n_steps, self.batch_size)
        try:
            result = subprocess.run(
                [self.reinvent_path, str(toml_path)],
                capture_output=True,
                text=True,
                timeout=3600,
            )
            if result.returncode != 0:
                logger.error("REINVENT4 failed (rc=%d): %s", result.returncode, result.stderr[:500])
                return pd.DataFrame()
        except FileNotFoundError:
            logger.error("REINVENT4 executable not found at '%s'", self.reinvent_path)
            return pd.DataFrame()
        except subprocess.TimeoutExpired:
            logger.error("REINVENT4 timed out after 1 hour.")
            return pd.DataFrame()

        return self._parse_output(output_csv)

    # ------------------------------------------------------------------
    @staticmethod
    def _write_scoring_script(
        out_dir: Path,
        model_predict_fn,
        featurizer_fn,
    ) -> Optional[Path]:
        """Write a thin Python script that REINVENT4 ExternalProcess can call.

        The script reads SMILES from stdin, computes predicted pIC50 via the
        serialised model artifacts, and writes scores to stdout.
        """
        try:
            import joblib
            model_pkl = out_dir / "reinvent4_scorer_model.pkl"
            feat_pkl = out_dir / "reinvent4_scorer_feat.pkl"
            joblib.dump(model_predict_fn, model_pkl)
            joblib.dump(featurizer_fn, feat_pkl)
        except Exception as exc:
            logger.warning("Could not serialise scoring model for REINVENT4: %s", exc)
            return None

        script = out_dir / "reinvent4_scorer.py"
        script.write_text(
            f"""#!/usr/bin/env python
\"\"\"REINVENT4 external scoring script — predicted pIC50.\"\"\"
import sys, joblib, numpy as np
model = joblib.load("{model_pkl}")
feat  = joblib.load("{feat_pkl}")
smiles = [line.strip() for line in sys.stdin if line.strip()]
if not smiles:
    sys.exit(0)
X = feat(smiles)
preds = model(X)
preds_norm = np.clip((preds - 5.0) / 5.0, 0.0, 1.0)
for v in preds_norm:
    print(f"{{v:.4f}}")
""",
            encoding="utf-8",
        )
        script.chmod(0o755)
        logger.info("REINVENT4 scoring script written: %s", script)
        return script

    # ------------------------------------------------------------------
    @staticmethod
    def _parse_output(csv_path: Path) -> pd.DataFrame:
        if not csv_path.exists():
            logger.warning("REINVENT4 output file not found: %s", csv_path)
            return pd.DataFrame()

        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            logger.error("Failed to read REINVENT4 output: %s", exc)
            return pd.DataFrame()

        smi_col = None
        for col in ("SMILES", "smiles", "Molecule", "canonical_smiles"):
            if col in df.columns:
                smi_col = col
                break
        if smi_col is None and len(df.columns) > 0:
            smi_col = df.columns[0]

        if smi_col is None:
            return pd.DataFrame()

        valid = []
        for smi in df[smi_col]:
            if safe_mol(str(smi)) is not None:
                valid.append(str(smi))

        out = pd.DataFrame({
            "canonical_smiles": valid,
            "origin": "reinvent4",
        }).drop_duplicates(subset=["canonical_smiles"])

        logger.info("REINVENT4 produced %d valid unique molecules.", len(out))
        return out
