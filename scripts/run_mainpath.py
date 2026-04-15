#!/usr/bin/env python3
"""Convenience runner for variant-driven main path (no-overwrite defaults).

Examples:
  # VCF entry
  python scripts/run_mainpath.py \
    --disease "breast cancer" \
    --sample-id patient_001 \
    --vcf-path /data/patient_001.ann.vcf.gz \
    -v

  # FASTQ entry
  python scripts/run_mainpath.py \
    --disease "lung cancer" \
    --sample-id patient_002 \
    --tumor-r1 /data/tumor_R1.fastq.gz \
    --tumor-r2 /data/tumor_R2.fastq.gz \
    --normal-r1 /data/normal_R1.fastq.gz \
    --normal-r2 /data/normal_R2.fastq.gz \
    -v
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Any, Dict

from drugpipe.config import load_config
from drugpipe.paths import run_root_for_config
from drugpipe.pipeline import _setup_file_logging, run_pipeline


def _build_overrides(args: argparse.Namespace) -> Dict[str, Any]:
    overrides: Dict[str, Any] = {
        "pipeline": {
            "out_dir": args.out_dir,
            "stages": ["lead_optimization"],
            "output_layout": {
                "variant_isolated_runs": True,
            },
        },
        "target_discovery": {
            "disease": args.disease,
        },
        "variant_analysis": {
            "enabled": True,
            "sample_id": args.sample_id or "",
            "auto_stage23": not args.no_auto_stage23,
            "sarek": {
                "profile": args.sarek_profile,
            },
        },
    }

    if args.device:
        overrides.setdefault("pipeline", {})["device"] = args.device
    if args.vcf_path:
        overrides["variant_analysis"]["vcf_path"] = args.vcf_path
    else:
        overrides["variant_analysis"]["vcf_path"] = ""
        overrides["variant_analysis"]["tumor_fastq_r1"] = args.tumor_r1
        overrides["variant_analysis"]["tumor_fastq_r2"] = args.tumor_r2
        overrides["variant_analysis"]["normal_fastq_r1"] = args.normal_r1
        overrides["variant_analysis"]["normal_fastq_r2"] = args.normal_r2

    if args.no_share_chembl_cache:
        overrides.setdefault("target_to_hit", {})["shared_library_dir"] = ""

    md_cfg = overrides.setdefault("lead_optimization", {}).setdefault("md_simulation", {})
    if args.parallel_workers is not None:
        md_cfg["parallel_workers"] = args.parallel_workers
    if args.schedule_tasks_across_gpus is not None:
        md_cfg["schedule_tasks_across_gpus"] = args.schedule_tasks_across_gpus
    if args.gpu_device_indices:
        md_cfg["gpu_device_indices"] = args.gpu_device_indices

    if args.xena_sample_strategy or args.xena_sample_id:
        xa = overrides.setdefault("variant_analysis", {}).setdefault("xena", {})
        if args.xena_sample_strategy:
            xa["sample_strategy"] = args.xena_sample_strategy
        if args.xena_sample_id:
            xa["sample_id"] = args.xena_sample_id
    return overrides


def _validate_args(args: argparse.Namespace) -> None:
    if args.vcf_path:
        return
    required_fastq = [args.tumor_r1, args.tumor_r2, args.normal_r1, args.normal_r2]
    if not all(required_fastq):
        raise SystemExit(
            "Either provide --vcf-path OR provide all FASTQs: "
            "--tumor-r1/--tumor-r2/--normal-r1/--normal-r2."
        )


def _parse_gpu_device_indices(raw: str) -> list[int]:
    if not raw.strip():
        return []
    try:
        values = [int(x.strip()) for x in raw.split(",") if x.strip()]
    except ValueError as exc:
        raise SystemExit("--gpu-device-indices must be a comma-separated integer list, e.g. 0,1,2") from exc
    if any(v < 0 for v in values):
        raise SystemExit("--gpu-device-indices only accepts non-negative integers.")
    return values


def main() -> None:
    env_parallel_workers = os.getenv("PARALLEL_WORKERS", "").strip()
    default_parallel_workers = None
    if env_parallel_workers:
        try:
            default_parallel_workers = int(env_parallel_workers)
        except ValueError:
            raise SystemExit("Invalid PARALLEL_WORKERS env; expected integer (e.g. 1).")

    parser = argparse.ArgumentParser(
        description="Run variant-driven main path with isolated outputs by default.",
    )
    parser.add_argument("--disease", required=True, help='e.g. "breast cancer"')
    parser.add_argument(
        "--sample-id",
        default="",
        help="Patient/sample ID used in output isolation path.",
    )
    parser.add_argument(
        "--out-dir",
        default="/root/autodl-fs/T2Lead_mainpath",
        help="Output root for main path runs.",
    )
    parser.add_argument(
        "--device",
        default="",
        choices=["", "auto", "cuda", "cpu", "mps"],
        help="Override pipeline.device (auto|cuda|cpu|mps).",
    )
    parser.add_argument(
        "--vcf-path",
        default="",
        help="VEP VCF, or Xena somatic mutation TSV (auto-converted when input_auto_detect is on).",
    )
    parser.add_argument("--tumor-r1", default="", help="Tumor FASTQ R1")
    parser.add_argument("--tumor-r2", default="", help="Tumor FASTQ R2")
    parser.add_argument("--normal-r1", default="", help="Normal FASTQ R1")
    parser.add_argument("--normal-r2", default="", help="Normal FASTQ R2")
    parser.add_argument(
        "--sarek-profile",
        default="singularity",
        choices=["docker", "singularity", "conda"],
        help="Nextflow profile for Sarek FASTQ calling.",
    )
    parser.add_argument(
        "--no-auto-stage23",
        action="store_true",
        help="Disable per-gene automatic Stage 2-3; reuse stage3 leads instead.",
    )
    parser.add_argument(
        "--share-chembl-cache",
        action="store_true",
        help="(兼容旧参数) 默认已经共享 ChEMBL 缓存，可忽略此参数。",
    )
    parser.add_argument(
        "--no-share-chembl-cache",
        action="store_true",
        help="Disable shared ChEMBL crawl/fp_cache and use run-local Stage2 cache.",
    )
    parser.add_argument(
        "--xena-sample-strategy",
        default="",
        choices=["", "max_variants", "explicit"],
        help="When --vcf-path is a Xena TSV: max_variants (default in YAML) or explicit.",
    )
    parser.add_argument(
        "--xena-sample-id",
        default="",
        help="When strategy is explicit: TCGA sample id (must exist in TSV sample column).",
    )
    parser.add_argument(
        "--parallel-workers",
        type=int,
        default=default_parallel_workers,
        help=(
            "Override lead_optimization.md_simulation.parallel_workers. "
            "If unset, keep YAML value. Also supports PARALLEL_WORKERS env compatibility."
        ),
    )
    parser.add_argument(
        "--schedule-tasks-across-gpus",
        dest="schedule_tasks_across_gpus",
        action="store_true",
        help="Enable MD task round-robin across multiple GPUs.",
    )
    parser.add_argument(
        "--no-schedule-tasks-across-gpus",
        dest="schedule_tasks_across_gpus",
        action="store_false",
        help="Disable MD task round-robin across GPUs.",
    )
    parser.set_defaults(schedule_tasks_across_gpus=None)
    parser.add_argument(
        "--gpu-device-indices",
        default="",
        help="Comma-separated CUDA logical devices for MD (e.g. 0,1,2).",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    args.gpu_device_indices = _parse_gpu_device_indices(args.gpu_device_indices)

    _validate_args(args)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    logging.basicConfig(level=log_level, format=log_fmt)

    overrides = _build_overrides(args)
    cfg = load_config(config_path="/root/T2Lead/configs/default_config.yaml", overrides=overrides)
    _setup_file_logging(cfg, log_level, log_fmt)

    run_root = run_root_for_config(cfg)
    print(f"Resolved run root: {Path(run_root)}")
    run_pipeline(cfg)


if __name__ == "__main__":
    main()

