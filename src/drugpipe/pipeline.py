"""Four-stage pipeline orchestrator: Target Discovery → T2H → H2L → Lead Optimization."""
# 四阶段流水线编排器：靶点发现 → 靶点到苗头 → 苗头到先导 → 先导优化。

from __future__ import annotations

import logging
import random
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from rdkit import RDLogger

from drugpipe.config import get_out_dir, load_config

RDLogger.DisableLog("rdApp.warning")

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
# Target viability check
# 靶点可用性检查
# ======================================================================

def _pick_viable_target(
    cfg: Dict[str, Any],
    targets: List[Dict[str, Any]],
    out_dir: Path,
) -> str:
    """Pick the first target from *targets* that has enough IC50 records.

    Falls back to the top-ranked target if activity data is not yet available
    (i.e. first run before crawl).
    """
    act_csv = out_dir / "activities_ic50.csv"
    min_samples = int(
        cfg.get("target_to_hit", {}).get("dataset", {}).get("min_train_samples", 200)
    )

    if not act_csv.exists():
        picked = targets[0]["chembl_id"]
        logger.info("Activity data not yet crawled; using top-ranked target %s", picked)
        return picked

    df_act = pd.read_csv(act_csv, usecols=["target_chembl_id"])
    counts = df_act["target_chembl_id"].value_counts()

    for t in targets:
        cid = t["chembl_id"]
        n = counts.get(cid, 0)
        symbol = t.get("symbol", "")
        logger.info("  %s (%s): %d IC50 records", cid, symbol, n)
        if n >= min_samples:
            logger.info(
                "Selected target %s (%s) — %d samples >= minimum %d",
                cid, symbol, n, min_samples,
            )
            return cid

    picked = targets[0]["chembl_id"]
    logger.warning(
        "No target met the minimum %d samples. Falling back to top-ranked %s.",
        min_samples, picked,
    )
    return picked


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
    crawl_out_dir: Optional[Path] = None,
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
    shared_dir = crawl_out_dir or out_dir
    seed = int(cfg.get("pipeline", {}).get("seed", 42))
    _set_seed(seed)

    # Resolve target ID / 解析靶点 ID
    tid = target_chembl_id or cfg.get("target_to_hit", {}).get("target_chembl_id", "") or None

    # 1. Crawl (shared directory) / 爬取（共享目录）
    crawler = ChEMBLCrawler(cfg, shared_dir)
    mol_csv = crawler.crawl_molecules()
    act_csv = crawler.crawl_activities()

    # 2. Dataset (per-disease directory) / 数据集构建（按疾病子目录）
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

    # Optional REINVENT4 enrichment / 可选 REINVENT4 强化学习 enrichment
    r4 = Reinvent4Bridge(cfg)
    if r4.enabled:
        df_r4 = r4.run(df_hits, out_dir, model_predict_fn=model_predict_fn, featurizer_fn=featurizer_fn)
        if not df_r4.empty:
            logger.info("Merging %d REINVENT4 molecules into lead pool.", len(df_r4))
            df_leads = pd.concat([df_leads, df_r4], ignore_index=True).drop_duplicates(
                subset=["canonical_smiles"]
            )

    logger.info("Stage 3 complete — %d lead candidates.", len(df_leads))
    return df_leads


def run_lead_optimization(
    cfg: Dict[str, Any],
    df_leads: pd.DataFrame,
) -> pd.DataFrame:
    """Stage 4: protein prep → docking → enhanced ADMET → MD → ranking."""
    # 阶段四：蛋白准备 → 对接 → 增强 ADMET → MD 模拟 → 综合排序。
    from drugpipe.lead_optimization.lead_optimizer import LeadOptimizer

    out_dir = get_out_dir(cfg)
    optimizer = LeadOptimizer(cfg, out_dir)
    df_optimized = optimizer.run(df_leads)

    logger.info("Stage 4 complete — %d optimized lead candidates.", len(df_optimized))
    return df_optimized


# ======================================================================
# Full pipeline
# 完整流水线
# ======================================================================

def _disease_slug(disease: str) -> str:
    """Sanitize disease name into a filesystem-safe directory name."""
    import re
    slug = disease.strip().lower()
    slug = re.sub(r"[^\w\s-]", "", slug)
    slug = re.sub(r"[\s]+", "_", slug)
    return slug or "default"


def run_pipeline(cfg: Dict[str, Any]) -> None:
    """Execute the stages listed in ``pipeline.stages``."""
    # 执行配置中 pipeline.stages 所列阶段。
    stages = cfg.get("pipeline", {}).get("stages", [])
    seed = int(cfg.get("pipeline", {}).get("seed", 42))
    _set_seed(seed)

    base_out_dir = get_out_dir(cfg)

    disease = cfg.get("target_discovery", {}).get("disease", "")
    if disease:
        run_out_dir = base_out_dir / _disease_slug(disease)
        run_out_dir.mkdir(parents=True, exist_ok=True)
    else:
        run_out_dir = base_out_dir

    logger.info("Pipeline base output: %s", base_out_dir)
    logger.info("Run-specific output : %s", run_out_dir)

    target_chembl_id: Optional[str] = None
    df_hits: Optional[pd.DataFrame] = None
    df_leads: Optional[pd.DataFrame] = None
    trainer = None
    featurizer = None

    # Stage 1
    all_targets: List[Dict[str, Any]] = []
    if "target_discovery" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 1: Target Discovery")
        logger.info("=" * 60)
        all_targets = run_target_discovery(cfg)
        if all_targets:
            target_chembl_id = all_targets[0]["chembl_id"]
            logger.info("Top-ranked target: %s", target_chembl_id)
        else:
            logger.warning("No targets found. Stage 2 will use config target_chembl_id or auto-select.")

    # Stage 2
    if "target_to_hit" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 2: Target to Hit")
        logger.info("=" * 60)

        if all_targets and "target_discovery" in stages:
            # Crawl first so _pick_viable_target can check actual IC50 counts
            act_csv = base_out_dir / "activities_ic50.csv"
            if not act_csv.exists():
                logger.info("Crawling ChEMBL data before target selection ...")
                from drugpipe.target_to_hit.chembl_api import ChEMBLCrawler
                crawler = ChEMBLCrawler(cfg, base_out_dir)
                crawler.crawl_molecules()
                crawler.crawl_activities()
            target_chembl_id = _pick_viable_target(cfg, all_targets, base_out_dir)

        # Crawl data is shared (base_out_dir); per-disease outputs go to run_out_dir
        cfg_s2 = _override_out_dir(cfg, run_out_dir)
        df_hits = run_target_to_hit(cfg_s2, target_chembl_id=target_chembl_id,
                                     crawl_out_dir=base_out_dir)
        trainer = df_hits.attrs.get("_trainer")
        featurizer = df_hits.attrs.get("_featurizer")

    # Stage 3
    if "hit_to_lead" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 3: Hit to Lead")
        logger.info("=" * 60)
        if df_hits is None:
            hit_csv = run_out_dir / "final_hit_candidates.csv"
            if hit_csv.exists():
                df_hits = pd.read_csv(hit_csv)
                logger.info("Loaded hits from %s (%d rows)", hit_csv, len(df_hits))
            else:
                logger.error("No hit candidates available. Run Stage 2 first.")
                return

        cfg_s3 = _override_out_dir(cfg, run_out_dir)
        predict_fn = trainer.predict if trainer else None
        feat_fn = featurizer.transform if featurizer else None
        df_leads = run_hit_to_lead(cfg_s3, df_hits, model_predict_fn=predict_fn, featurizer_fn=feat_fn)

    # Stage 4
    if "lead_optimization" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 4: Lead Optimization")
        logger.info("=" * 60)
        if df_leads is None:
            leads_csv = run_out_dir / "final_lead_candidates.csv"
            if leads_csv.exists():
                df_leads = pd.read_csv(leads_csv)
                logger.info("Loaded leads from %s (%d rows)", leads_csv, len(df_leads))
            else:
                logger.error("No lead candidates available. Run Stage 3 first.")
                return

        cfg_s4 = _override_out_dir(cfg, run_out_dir)
        run_lead_optimization(cfg_s4, df_leads)

    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETE")
    logger.info("=" * 60)


def _override_out_dir(cfg: Dict[str, Any], out_dir: Path) -> Dict[str, Any]:
    """Return a shallow config copy with pipeline.out_dir overridden."""
    import copy
    cfg2 = copy.deepcopy(cfg)
    cfg2.setdefault("pipeline", {})["out_dir"] = str(out_dir)
    return cfg2


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
