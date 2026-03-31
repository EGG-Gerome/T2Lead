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
# Optional REINVENT4 integration for RL-based molecular optimization.
# 说明模块职责、上下游关系与维护注意事项。

# 可选 REINVENT4 集成：基于强化学习的分子优化；生成配置、以子进程启动并收集输出 SMILES。

from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)


def _resolve_reinvent_device(cfg_device: str) -> str:
    """Resolve ``pipeline.device`` to a PyTorch device string for REINVENT4.

    Priority: cuda > mps > cpu (same semantics as pipeline-wide device).
    REINVENT4 passes this string directly to ``torch.device()``.
    """
    cfg_device = cfg_device.lower().strip()

    try:
        import torch
    except ImportError:
        return "cpu"

    if cfg_device in ("cuda", "auto"):
        if torch.cuda.is_available():
            return "cuda:0"
        if cfg_device == "cuda":
            logger.warning("CUDA requested but unavailable; REINVENT4 falling back to CPU.")

    if cfg_device in ("mps", "auto"):
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
        if cfg_device == "mps":
            logger.warning("MPS requested but unavailable; REINVENT4 falling back to CPU.")

    return "cpu"


_TOML_TEMPLATE = """\
run_type = "staged_learning"
device = "{device}"

[parameters]
prior_file = "{prior_path}"
agent_file = "{prior_path}"
summary_csv_prefix = "{output_prefix}"
batch_size = {batch_size}
use_checkpoint = false
unique_sequences = true
randomize_smiles = true

[learning_strategy]
type = "dap"
sigma = 128
rate = 0.0001

[diversity_filter]
type = "IdenticalMurckoScaffold"
bucket_size = 25
minscore = 0.4
minsimilarity = 0.4

[[stage]]
max_score = 1.0
min_steps = 10
max_steps = {n_steps}
chkpt_file = "{chkpt_file}"
termination = "simple"

[stage.scoring]
type = "geometric_mean"

[[stage.scoring.component]]
[stage.scoring.component.QED]
[[stage.scoring.component.QED.endpoint]]
name = "QED"
weight = 1.0

[[stage.scoring.component]]
[stage.scoring.component.TanimotoSimilarity]
[[stage.scoring.component.TanimotoSimilarity.endpoint]]
name = "Tanimoto to seed hit"
weight = 0.6
params.smiles = ["{reference_smiles}"]
params.radius = 2
params.use_counts = false
params.use_features = false

{predictive_model_block}
"""


_PREDICTIVE_MODEL_BLOCK = """\
[[stage.scoring.component]]
[stage.scoring.component.ExternalProcess]
[[stage.scoring.component.ExternalProcess.endpoint]]
name = "pIC50 predictor"
weight = 0.8
params.executable = "{python_exe}"
params.args = "{scoring_script}"
params.property = "pIC50"
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
        self.log_level = (r4.get("log_level") or "warning").strip().lower()

        pipeline_device = cfg.get("pipeline", {}).get("device", "auto")
        self.device = _resolve_reinvent_device(pipeline_device)

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
                pred_block = _PREDICTIVE_MODEL_BLOCK.format(
                    python_exe=sys.executable,
                    scoring_script=str(scoring_script),
                )

        output_prefix = str(out_dir / "reinvent4")
        chkpt_file = out_dir / "reinvent4_agent.ckpt"
        toml_content = _TOML_TEMPLATE.format(
            device=self.device,
            prior_path=self.prior_path,
            output_prefix=output_prefix,
            batch_size=self.batch_size,
            n_steps=self.n_steps,
            chkpt_file=str(chkpt_file),
            reference_smiles=ref_smi,
            predictive_model_block=pred_block,
        )

        toml_path = out_dir / "reinvent4_config.toml"
        toml_path.write_text(toml_content, encoding="utf-8")

        logger.info(
            "Launching REINVENT4 (device=%s, steps=%d, batch=%d, log_level=%s) ...",
            self.device,
            self.n_steps,
            self.batch_size,
            self.log_level,
        )
        try:
            env = os.environ.copy()
            conda_lib = str(Path(sys.executable).parent.parent / "lib")
            env["LD_LIBRARY_PATH"] = conda_lib + ":" + env.get("LD_LIBRARY_PATH", "")
            result = subprocess.run(
                [self.reinvent_path, "--log-level", self.log_level, str(toml_path)],
                capture_output=True,
                text=True,
                timeout=3600,
                env=env,
            )
            if result.returncode != 0:
                logger.error(
                    "REINVENT4 failed (rc=%d):\nSTDOUT: %s\nSTDERR: %s",
                    result.returncode,
                    result.stdout[-1000:] if result.stdout else "",
                    result.stderr[-2000:] if result.stderr else "",
                )
                return pd.DataFrame()
        except FileNotFoundError:
            logger.error("REINVENT4 executable not found at '%s'", self.reinvent_path)
            return pd.DataFrame()
        except subprocess.TimeoutExpired:
            logger.error("REINVENT4 timed out after 1 hour.")
            return pd.DataFrame()

        output_csv = out_dir / "reinvent4_1.csv"
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
import sys, json, joblib, numpy as np
model = joblib.load("{model_pkl}")
feat  = joblib.load("{feat_pkl}")
smiles = [line.strip() for line in sys.stdin if line.strip()]
if not smiles:
    print(json.dumps({{"payload": {{"pIC50": []}}}})  )
    sys.exit(0)
X = feat(smiles)
preds = model(X)
scores = np.clip((preds - 5.0) / 5.0, 0.0, 1.0).tolist()
print(json.dumps({{"payload": {{"pIC50": scores}}}}))
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
