"""Three-stage pipeline orchestrator: Target Discovery → T2H → H2L."""
# 三阶段流水线编排器：靶点发现 → 靶点到先导化合物 → 先导化合物优化。

from __future__ import annotations

import logging
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from drugpipe.config import get_out_dir, load_config

logger = logging.getLogger(__name__)


def _set_seed(seed: int) -> None:
    """Set global and numpy/torch RNG seed. / 设置全局及 numpy/torch 随机种子。"""
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


# ======================================================================
# Stage runners
# 各阶段执行函数
# ======================================================================

def run_target_discovery(cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Stage 1: identify therapeutic targets for a disease."""
    # 阶段一：识别疾病的治疗靶点。
    from drugpipe.target_discovery.target_ranker import TargetRanker

    ranker = TargetRanker(cfg)
    targets = ranker.run()
    logger.info("Stage 1 complete — %d targets identified.", len(targets))
    return targets


def run_target_to_hit(
    cfg: Dict[str, Any],
    target_chembl_id: Optional[str] = None,
) -> pd.DataFrame:
    """Stage 2: ChEMBL crawl → ML training → virtual screening → ADMET filter."""
    # 阶段二：ChEMBL 爬取 → ML 训练 → 虚拟筛选 → ADMET 过滤。
    from drugpipe.target_to_hit.chembl_api import ChEMBLCrawler
    from drugpipe.target_to_hit.dataset import DatasetBuilder
    from drugpipe.target_to_hit.featurizer import MorganFeaturizer
    from drugpipe.target_to_hit.filters import ADMETFilter
    from drugpipe.target_to_hit.models import ModelTrainer
    from drugpipe.target_to_hit.screener import VirtualScreener

    out_dir = get_out_dir(cfg)
    seed = int(cfg.get("pipeline", {}).get("seed", 42))
    _set_seed(seed)

    # Resolve target ID / 解析靶点 ID
    tid = target_chembl_id or cfg.get("target_to_hit", {}).get("target_chembl_id", "") or None

    # 1. Crawl / 爬取
    crawler = ChEMBLCrawler(cfg, out_dir)
    mol_csv = crawler.crawl_molecules()
    act_csv = crawler.crawl_activities()

    # 2. Dataset / 数据集构建
    builder = DatasetBuilder(cfg, out_dir)
    df_ds = builder.build(mol_csv, act_csv, tid)

    # 3. Featurize / 特征化
    feat = MorganFeaturizer(cfg)
    X = feat.transform(df_ds["canonical_smiles"].tolist())
    y = df_ds["pIC50_median"].astype(float).values

    # 4. Train / 训练
    trainer = ModelTrainer(cfg)
    metrics = trainer.train(X, y)
    logger.info("Training metrics: %s", metrics)

    # 5. Virtual screening on all molecules / 对全库分子虚拟筛选
    df_mol = pd.read_csv(mol_csv).dropna(subset=["canonical_smiles"]).drop_duplicates("molecule_chembl_id")
    screener = VirtualScreener(cfg, out_dir)
    df_scored = screener.run(df_mol, trainer)

    # 6. ADMET filter / ADMET 过滤
    filt = ADMETFilter(cfg, out_dir)
    df_hits = filt.run(df_scored)

    logger.info("Stage 2 complete — %d hit candidates.", len(df_hits))

    # Stash trainer/featurizer on the returned frame for Stage 3 re-use
    # 将 trainer/featurizer 挂在返回的 DataFrame 上供阶段三复用
    df_hits.attrs["_trainer"] = trainer
    df_hits.attrs["_featurizer"] = feat

    return df_hits


def run_hit_to_lead(
    cfg: Dict[str, Any],
    df_hits: pd.DataFrame,
    model_predict_fn=None,
    featurizer_fn=None,
) -> pd.DataFrame:
    """Stage 3: scaffold analysis → clustering → analog gen → MPO → leads."""
    # 阶段三：骨架分析 → 聚类 → 类似物生成 → MPO → 先导输出。
    from drugpipe.hit_to_lead.lead_ranker import LeadRanker
    from drugpipe.hit_to_lead.reinvent_bridge import Reinvent4Bridge

    out_dir = get_out_dir(cfg)
    ranker = LeadRanker(cfg, out_dir)
    df_leads = ranker.run(df_hits, model_predict_fn=model_predict_fn, featurizer_fn=featurizer_fn)

    # Optional REINVENT4 enrichment / 可选 REINVENT4  enrichment
    r4 = Reinvent4Bridge(cfg)
    if r4.enabled:
        df_r4 = r4.run(df_hits, out_dir)
        if not df_r4.empty:
            logger.info("Merging %d REINVENT4 molecules into lead pool.", len(df_r4))
            df_leads = pd.concat([df_leads, df_r4], ignore_index=True).drop_duplicates(
                subset=["canonical_smiles"]
            )

    logger.info("Stage 3 complete — %d lead candidates.", len(df_leads))
    return df_leads


# ======================================================================
# Full pipeline
# 完整流水线
# ======================================================================

def run_pipeline(cfg: Dict[str, Any]) -> None:
    """Execute the stages listed in ``pipeline.stages``."""
    # 执行配置中 pipeline.stages 所列阶段。
    stages = cfg.get("pipeline", {}).get("stages", [])
    seed = int(cfg.get("pipeline", {}).get("seed", 42))
    _set_seed(seed)

    out_dir = get_out_dir(cfg)
    logger.info("Pipeline output directory: %s", out_dir)

    target_chembl_id: Optional[str] = None
    df_hits: Optional[pd.DataFrame] = None
    trainer = None
    featurizer = None

    # Stage 1
    if "target_discovery" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 1: Target Discovery")
        logger.info("=" * 60)
        targets = run_target_discovery(cfg)
        if targets:
            target_chembl_id = targets[0]["chembl_id"]
            logger.info("Selected primary target: %s", target_chembl_id)
        else:
            logger.warning("No targets found. Stage 2 will use config target_chembl_id or auto-select.")

    # Stage 2
    if "target_to_hit" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 2: Target to Hit")
        logger.info("=" * 60)
        df_hits = run_target_to_hit(cfg, target_chembl_id=target_chembl_id)
        trainer = df_hits.attrs.get("_trainer")
        featurizer = df_hits.attrs.get("_featurizer")

    # Stage 3
    if "hit_to_lead" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 3: Hit to Lead")
        logger.info("=" * 60)
        if df_hits is None:
            hit_csv = out_dir / "final_hit_candidates.csv"
            if hit_csv.exists():
                df_hits = pd.read_csv(hit_csv)
                logger.info("Loaded hits from %s (%d rows)", hit_csv, len(df_hits))
            else:
                logger.error("No hit candidates available. Run Stage 2 first.")
                return

        predict_fn = trainer.predict if trainer else None
        feat_fn = featurizer.transform if featurizer else None
        run_hit_to_lead(cfg, df_hits, model_predict_fn=predict_fn, featurizer_fn=feat_fn)

    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETE")
    logger.info("=" * 60)


# ======================================================================
# CLI entry point
# 命令行入口
# ======================================================================

def main() -> None:
    """Entry point for ``drugpipe`` console script."""
    # drugpipe 控制台脚本入口。
    import argparse

    parser = argparse.ArgumentParser(
        description="DrugPipe: end-to-end drug discovery pipeline",
    )
    parser.add_argument(
        "-c", "--config",
        default=None,
        help="Path to a custom YAML config file (merges with default).",
    )
    parser.add_argument(
        "--stages",
        nargs="*",
        default=None,
        help="Override pipeline.stages (e.g. --stages target_to_hit hit_to_lead).",
    )
    parser.add_argument(
        "--disease",
        default=None,
        help="Disease name for Stage 1 target discovery.",
    )
    parser.add_argument(
        "--target",
        default=None,
        help="ChEMBL target ID for Stage 2 (skip Stage 1).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable DEBUG logging.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    overrides: Dict[str, Any] = {}
    if args.stages is not None:
        overrides.setdefault("pipeline", {})["stages"] = args.stages
    if args.disease:
        overrides.setdefault("target_discovery", {})["disease"] = args.disease
    if args.target:
        overrides.setdefault("target_to_hit", {})["target_chembl_id"] = args.target

    cfg = load_config(config_path=args.config, overrides=overrides)
    run_pipeline(cfg)


if __name__ == "__main__":
    main()
