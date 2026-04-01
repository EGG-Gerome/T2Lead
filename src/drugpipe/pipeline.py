"""Four-stage pipeline orchestrator: Target Discovery → T2H → H2L → Lead Optimization."""
# Four-stage pipeline orchestrator: Target Discovery → T2H → H2L → Lead Optimization.
# 说明模块职责、上下游关系与维护注意事项。

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
from drugpipe.paths import (
    STAGE1,
    STAGE2,
    STAGE3,
    STAGE4,
    chembl_library_root,
    disease_slug,
    resolve_variant_stage4_out,
    run_root_for_config,
    stage_paths,
)

RDLogger.DisableLog("rdApp.warning")

logger = logging.getLogger(__name__)


def _set_seed(seed: int) -> None:
    """Set deterministic seeds for Python, NumPy and Torch.
    为 Python、NumPy 和 Torch 设置可复现随机种子。
    """
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
    """Pick the first target with enough IC50 training records.
    选择第一个满足 IC50 样本量阈值的候选靶点。

    If activity data has not been crawled yet, fall back to the top-ranked
    target so the first run can proceed.
    若活性数据尚未抓取（首次运行），则回退到排序第一的靶点以继续流程。
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
    """Stage 1 - identify disease-relevant therapeutic targets.
    阶段 1：识别与疾病相关的治疗靶点。
    """
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
    """Stage 2 - data acquisition, model training, screening and filtering.
    阶段 2：数据获取、模型训练、虚拟筛选与 ADMET 过滤。

    Supported modes:
      1) ChEMBL mode (default): crawl molecules and IC50 activities from ChEMBL.
      2) External data mode: use user IC50 CSV via `target_to_hit.external_activities_csv`.
      3) Docking-only mode: skip ML and pass drug-like candidates to Stage 4 docking.
    支持三种模式：
      1）ChEMBL 模式（默认）：抓取分子与 IC50 活性数据；
      2）外部数据模式：通过 `target_to_hit.external_activities_csv` 使用用户数据；
      3）仅对接模式：跳过 ML，按成药性筛选后交给阶段 4 对接评分。
    """
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

    t2h_cfg = cfg.get("target_to_hit", {})
    tid = target_chembl_id or t2h_cfg.get("target_chembl_id", "") or None
    docking_only = bool(t2h_cfg.get("docking_only", False))
    ext_act_csv = t2h_cfg.get("external_activities_csv", "") or ""
    ext_lib_csv = t2h_cfg.get("screening_library_csv", "") or ""

    # --- Acquire molecule library ---
    if ext_lib_csv and Path(ext_lib_csv).exists():
        mol_csv = Path(ext_lib_csv)
        logger.info("Using user-supplied screening library: %s", mol_csv)
    else:
        crawler = ChEMBLCrawler(cfg, shared_dir)
        mol_csv = crawler.crawl_molecules()

    # --- Docking-only mode: skip ML, filter by drug-likeness only ---
    if docking_only:
        logger.info("Docking-only mode: skipping ML model training.")
        logger.info("Candidates will be selected by drug-likeness (QED + ADMET rules).")
        df_mol = pd.read_csv(mol_csv).dropna(subset=["canonical_smiles"]).drop_duplicates("molecule_chembl_id")
        filt = ADMETFilter(cfg, out_dir)
        screener = VirtualScreener(cfg, out_dir)
        screener._add_properties(df_mol)
        df_mol["pred_pIC50_ens"] = np.nan
        df_mol["pred_IC50_nM_ens"] = np.nan
        df_hits = filt.run(df_mol)
        logger.info("Stage 2 (docking-only) complete — %d candidates for docking.", len(df_hits))
        return df_hits

    # --- Acquire activity data ---
    if ext_act_csv and Path(ext_act_csv).exists():
        act_csv = Path(ext_act_csv)
        logger.info("Using user-supplied activity data: %s", act_csv)
    else:
        crawler = ChEMBLCrawler(cfg, shared_dir)
        act_csv = crawler.crawl_activities()

    # --- ML-based screening pipeline ---
    fp_cache_dir = shared_dir / "fp_cache"
    builder = DatasetBuilder(cfg, out_dir)
    df_ds = builder.build(mol_csv, act_csv, tid)

    feat = MorganFeaturizer(cfg, cache_dir=fp_cache_dir)
    X = feat.transform(df_ds["canonical_smiles"].tolist())
    y = df_ds["pIC50_median"].astype(float).values

    model_cache_dir = out_dir / "model_cache"
    trainer = ModelTrainer(cfg)
    metrics = trainer.train(X, y, cache_dir=model_cache_dir)
    logger.info("Training metrics: %s", metrics)

    df_mol = pd.read_csv(mol_csv).dropna(subset=["canonical_smiles"]).drop_duplicates("molecule_chembl_id")
    screener = VirtualScreener(cfg, out_dir, fp_cache_dir=fp_cache_dir)
    df_scored = screener.run(df_mol, trainer)

    filt = ADMETFilter(cfg, out_dir)
    df_hits = filt.run(df_scored)

    logger.info("Stage 2 complete — %d hit candidates.", len(df_hits))

    df_hits.attrs["_trainer"] = trainer
    df_hits.attrs["_featurizer"] = feat

    return df_hits


def run_hit_to_lead(
    cfg: Dict[str, Any],
    df_hits: pd.DataFrame,
    model_predict_fn=None,
    featurizer_fn=None,
) -> pd.DataFrame:
    """Stage 3 - scaffold analysis, clustering, analog generation and MPO ranking.
    阶段 3：骨架分析、聚类、类似物生成与 MPO 排序，输出先导分子。
    """
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
            df_merged = pd.concat([df_leads, df_r4], ignore_index=True).drop_duplicates(
                subset=["canonical_smiles"]
            )

            from drugpipe.hit_to_lead.mpo import MPOScorer
            mpo = MPOScorer(cfg)
            df_merged = mpo.score(
                df_merged,
                model_predict_fn=model_predict_fn,
                featurizer_fn=featurizer_fn,
            )

            top_n = int(
                cfg.get("hit_to_lead", {}).get("output", {}).get("top_n_leads", 50)
            )
            df_merged = df_merged.dropna(subset=["mpo_score"])
            df_merged = df_merged.sort_values("mpo_score", ascending=False)
            df_leads = df_merged.head(top_n).reset_index(drop=True)

            leads_csv = out_dir / "final_lead_candidates.csv"
            keep_cols = [
                "canonical_smiles", "source_smiles", "origin",
                "pred_pIC50_ens", "QED", "mpo_score",
                "admet_pass", "HasAlert", "scaffold_smi", "cluster_id",
            ]
            keep_cols = [c for c in keep_cols if c in df_leads.columns]
            df_leads[keep_cols].to_csv(leads_csv, index=False)
            logger.info(
                "Re-ranked leads after REINVENT4 merge: %s (%d rows)",
                leads_csv, len(df_leads),
            )

    logger.info("Stage 3 complete — %d lead candidates.", len(df_leads))
    return df_leads


def run_lead_optimization(
    cfg: Dict[str, Any],
    df_leads: pd.DataFrame,
) -> pd.DataFrame:
    """Stage 4 - protein prep, docking, ADMET, MD and final ranking.
    阶段 4：蛋白准备、对接、ADMET、MD 评估与最终排序。
    """
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

def _run_variant_analysis(cfg: Dict[str, Any], out_dir: Path) -> List[Dict[str, Any]]:
    """Execute the variant-analysis pathway: VCF parsing → mutant sequences → structures.

    Returns a list of structure records (gene, mutation, pdb_path, method).
    """
    from drugpipe.variant_analysis.vcf_parser import VCFParser
    from drugpipe.variant_analysis.mutant_sequence import MutantSequenceBuilder
    from drugpipe.variant_analysis.structure_bridge import StructureBridge

    va_cfg = cfg.get("variant_analysis", {})
    vcf_path = va_cfg.get("vcf_path", "")

    if not vcf_path:
        from drugpipe.variant_analysis.sarek_runner import SarekRunner
        runner = SarekRunner(cfg)
        t_r1 = va_cfg.get("tumor_fastq_r1", "")
        t_r2 = va_cfg.get("tumor_fastq_r2", "")
        n_r1 = va_cfg.get("normal_fastq_r1", "")
        n_r2 = va_cfg.get("normal_fastq_r2", "")
        if all([t_r1, t_r2, n_r1, n_r2]):
            logger.info("Running nf-core/sarek for somatic variant calling...")
            vcf_result = runner.run(t_r1, t_r2, n_r1, n_r2, out_dir / "variant_calling")
            if vcf_result:
                vcf_path = str(vcf_result)
            else:
                logger.error("Sarek variant calling failed.")
                return []
        else:
            logger.error("variant_analysis enabled but no VCF or FASTQs provided.")
            return []

    va_dir = out_dir / "variant_analysis"
    va_dir.mkdir(parents=True, exist_ok=True)

    parser = VCFParser(cfg)
    variants = parser.parse(vcf_path)
    if not variants:
        logger.warning("No protein-altering variants found in %s", vcf_path)
        return []

    for v in variants:
        logger.info("  Variant: %s", v.short_label)

    builder = MutantSequenceBuilder(cfg)
    mutants = builder.build(variants)
    if not mutants:
        logger.warning("Could not build any mutant sequences.")
        return []

    fasta_dir = va_dir / "mutant_fastas"
    builder.write_fastas(mutants, fasta_dir)

    bridge = StructureBridge(cfg)
    structures = bridge.resolve_structures(mutants, va_dir / "structures")
    return structures


def run_pipeline(cfg: Dict[str, Any]) -> None:
    """Execute stages listed in `pipeline.stages` with shared state handoff.
    按 `pipeline.stages` 顺序执行各阶段，并在阶段间传递中间结果。
    """
    stages = cfg.get("pipeline", {}).get("stages", [])
    seed = int(cfg.get("pipeline", {}).get("seed", 42))
    _set_seed(seed)

    base_out_dir = get_out_dir(cfg)
    run_out_dir = run_root_for_config(cfg)
    layout = stage_paths(run_out_dir, cfg)

    logger.info("Pipeline base output: %s", base_out_dir)
    logger.info("Run-specific output : %s", run_out_dir)

    target_chembl_id: Optional[str] = None
    df_hits: Optional[pd.DataFrame] = None
    df_leads: Optional[pd.DataFrame] = None
    trainer = None
    featurizer = None

    # Variant-driven path: if VCF/FASTQs provided, bypass Stages 1–3
    va_cfg = cfg.get("variant_analysis", {})
    if va_cfg.get("enabled", False):
        logger.info("=" * 60)
        logger.info("VARIANT ANALYSIS: Somatic mutation-driven pathway")
        logger.info("=" * 60)
        mutant_structures = _run_variant_analysis(cfg, run_out_dir)
        if mutant_structures:
            logger.info("Variant analysis produced %d mutant structures.", len(mutant_structures))
            if "lead_optimization" in stages:
                leads_csv = layout[STAGE3] / "final_lead_candidates.csv"
                if not leads_csv.exists():
                    logger.error(
                        "Variant-driven Stage 4 requires lead candidates at %s "
                        "(run Stages 2–3 first or copy a lead CSV there).",
                        leads_csv,
                    )
                    return
                df_vcf_leads = pd.read_csv(leads_csv)
                logger.info("Loaded %d leads from %s", len(df_vcf_leads), leads_csv)
                for ms in mutant_structures:
                    logger.info(
                        "Running Stage 4 against %s %s (%s)",
                        ms["gene"], ms["mutation"], ms["method"],
                    )
                    s4_out = resolve_variant_stage4_out(
                        run_out_dir, cfg, str(ms["gene"]), str(ms["mutation"]),
                    )
                    cfg_mut = _override_out_dir(cfg, s4_out)
                    lo = cfg_mut.setdefault("lead_optimization", {})
                    lo["pdb_id"] = ""
                    lo["protein_sequence"] = ""
                    lo["_mutant_pdb_path"] = ms["pdb_path"]
                    run_lead_optimization(cfg_mut, df_vcf_leads)
            logger.info("=" * 60)
            logger.info("VARIANT-DRIVEN PIPELINE COMPLETE")
            logger.info("=" * 60)
            return

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
            try:
                pd.DataFrame(all_targets).to_csv(
                    layout[STAGE1] / "ranked_targets.csv",
                    index=False,
                )
            except Exception as exc:
                logger.debug("Could not write ranked_targets.csv: %s", exc)
        else:
            logger.warning("No targets found. Stage 2 will use config target_chembl_id or auto-select.")

    # Stage 2
    if "target_to_hit" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 2: Target to Hit")
        logger.info("=" * 60)

        s2_local = layout[STAGE2]
        library_root = chembl_library_root(cfg, s2_local)
        if library_root != s2_local:
            logger.info("Shared ChEMBL library + fp_cache: %s", library_root)
            logger.info("Stage 2 disease-local outputs: %s", s2_local)

        if all_targets and "target_discovery" in stages:
            # Crawl first so _pick_viable_target can check actual IC50 counts
            act_csv = library_root / "activities_ic50.csv"
            if not act_csv.exists():
                logger.info("Crawling ChEMBL data before target selection ...")
                from drugpipe.target_to_hit.chembl_api import ChEMBLCrawler
                crawler = ChEMBLCrawler(cfg, library_root)
                crawler.crawl_molecules()
                crawler.crawl_activities()
            target_chembl_id = _pick_viable_target(cfg, all_targets, library_root)

        cfg_s2 = _override_out_dir(cfg, s2_local)
        df_hits = run_target_to_hit(
            cfg_s2,
            target_chembl_id=target_chembl_id,
            crawl_out_dir=library_root,
        )
        trainer = df_hits.attrs.get("_trainer")
        featurizer = df_hits.attrs.get("_featurizer")

    # Stage 3
    if "hit_to_lead" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 3: Hit to Lead")
        logger.info("=" * 60)
        if df_hits is None:
            hit_csv = layout[STAGE2] / "final_hit_candidates.csv"
            if hit_csv.exists():
                df_hits = pd.read_csv(hit_csv)
                logger.info("Loaded hits from %s (%d rows)", hit_csv, len(df_hits))
            else:
                logger.error("No hit candidates available. Run Stage 2 first.")
                return

        cfg_s3 = _override_out_dir(cfg, layout[STAGE3])
        predict_fn = trainer.predict if trainer else None
        feat_fn = featurizer.transform if featurizer else None
        df_leads = run_hit_to_lead(cfg_s3, df_hits, model_predict_fn=predict_fn, featurizer_fn=feat_fn)

    # Stage 4
    if "lead_optimization" in stages:
        logger.info("=" * 60)
        logger.info("STAGE 4: Lead Optimization")
        logger.info("=" * 60)
        if df_leads is None:
            leads_csv = layout[STAGE3] / "final_lead_candidates.csv"
            if leads_csv.exists():
                df_leads = pd.read_csv(leads_csv)
                logger.info("Loaded leads from %s (%d rows)", leads_csv, len(df_leads))
            else:
                logger.error("No lead candidates available. Run Stage 3 first.")
                return

        cfg_s4 = _override_out_dir(cfg, layout[STAGE4])
        run_lead_optimization(cfg_s4, df_leads)

    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETE")
    logger.info("=" * 60)


def _override_out_dir(cfg: Dict[str, Any], out_dir: Path) -> Dict[str, Any]:
    """Return a deep-copied config with `pipeline.out_dir` overridden.
    返回深拷贝后的配置，并覆写 `pipeline.out_dir`。
    """
    import copy
    cfg2 = copy.deepcopy(cfg)
    cfg2.setdefault("pipeline", {})["out_dir"] = str(out_dir)
    return cfg2


# ======================================================================
# CLI entry point
# 命令行入口
# ======================================================================

def main() -> None:
    """CLI entry point that parses args, builds config, and runs pipeline.
    命令行入口：解析参数、构建配置并启动流水线。
    """
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
        "--activities-csv",
        default=None,
        help="Path to user-supplied IC50 CSV for novel targets not in ChEMBL.",
    )
    parser.add_argument(
        "--screening-library",
        default=None,
        help="Path to user-supplied SMILES library CSV for screening.",
    )
    parser.add_argument(
        "--docking-only",
        action="store_true",
        help="Skip ML model, use docking scores only (for targets with no IC50 data).",
    )
    parser.add_argument(
        "--protein-sequence",
        default=None,
        help="Amino acid sequence for ESMFold structure prediction (no PDB needed).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable DEBUG logging.",
    )

    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"

    logging.basicConfig(level=log_level, format=log_fmt)

    overrides: Dict[str, Any] = {}
    if args.stages is not None:
        overrides.setdefault("pipeline", {})["stages"] = args.stages
    if args.disease:
        overrides.setdefault("target_discovery", {})["disease"] = args.disease
    if args.target:
        overrides.setdefault("target_to_hit", {})["target_chembl_id"] = args.target
    if args.activities_csv:
        overrides.setdefault("target_to_hit", {})["external_activities_csv"] = args.activities_csv
    if args.screening_library:
        overrides.setdefault("target_to_hit", {})["screening_library_csv"] = args.screening_library
    if args.docking_only:
        overrides.setdefault("target_to_hit", {})["docking_only"] = True
    if args.protein_sequence:
        overrides.setdefault("lead_optimization", {})["protein_sequence"] = args.protein_sequence
        if not args.target:
            overrides.setdefault("lead_optimization", {})["pdb_id"] = ""

    cfg = load_config(config_path=args.config, overrides=overrides)

    # --- Set up file logging ---
    _setup_file_logging(cfg, log_level, log_fmt)

    run_pipeline(cfg)


# ======================================================================
# File logging
# ======================================================================

def _setup_file_logging(
    cfg: Dict[str, Any], level: int, fmt: str,
) -> None:
    """Attach file handlers for detailed and summary logs.
    为详细日志与摘要日志绑定文件处理器。

    Output files:
      - `<out_dir>/logs/<timestamp>_..._full.log`
      - `<out_dir>/logs/<timestamp>_..._summary.log`
    输出文件：
      - `<out_dir>/logs/<时间戳>_..._full.log`
      - `<out_dir>/logs/<时间戳>_..._summary.log`
    """
    import datetime

    run_out = run_root_for_config(cfg)
    log_dir = stage_paths(run_out, cfg)["logs"]
    log_dir.mkdir(parents=True, exist_ok=True)

    disease = cfg.get("target_discovery", {}).get("disease", "")
    target = cfg.get("target_to_hit", {}).get("target_chembl_id", "")
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    parts = [ts]
    if disease:
        parts.append(disease_slug(disease))
    if target:
        parts.append(target)
    stem = "_".join(parts)

    root_logger = logging.getLogger()

    full_path = log_dir / f"{stem}_full.log"
    fh = logging.FileHandler(full_path, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(fmt))
    root_logger.addHandler(fh)

    summary_path = log_dir / f"{stem}_summary.log"
    sh = logging.FileHandler(summary_path, encoding="utf-8")
    sh.setLevel(logging.WARNING)
    sh.setFormatter(logging.Formatter(fmt))

    class _InfoKeywordFilter(logging.Filter):
        _KEYWORDS = (
            "Stage", "complete", "Selected target", "PIPELINE",
            "Training metrics", "Fold", "RMSE", "Docking progress",
            "Docking complete", "MD simulation", "MM-GBSA",
            "Hit candidates", "Lead candidates", "Optimized leads",
            "saved", "Loaded cached", "ESMFold",
        )
        def filter(self, record: logging.LogRecord) -> bool:
            if record.levelno >= logging.WARNING:
                return True
            return any(kw in record.getMessage() for kw in self._KEYWORDS)

    sh.addFilter(_InfoKeywordFilter())
    root_logger.addHandler(sh)

    logger.info("Full log: %s", full_path)
    logger.info("Summary log: %s", summary_path)


if __name__ == "__main__":
    main()
