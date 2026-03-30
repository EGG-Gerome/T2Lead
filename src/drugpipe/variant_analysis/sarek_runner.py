"""Wrapper for nf-core/sarek somatic variant calling pipeline.

Generates a Nextflow sample-sheet and launches ``nf-core/sarek`` with
Mutect2 + VEP annotation.  This module is optional — users can also
supply a pre-annotated VCF directly.

Requirements:
  - Nextflow >= 24.04
  - Singularity or Docker for containerised execution
  - Reference genome data (GRCh38 igenomes)
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

_SAREK_DEFAULT_VERSION = "3.8.1"


class SarekRunner:
    """Launch nf-core/sarek for tumor/normal somatic variant calling."""

    def __init__(self, cfg: Optional[Dict[str, Any]] = None):
        va = (cfg or {}).get("variant_analysis", {})
        sarek = va.get("sarek", {})
        self.version = sarek.get("version", _SAREK_DEFAULT_VERSION)
        self.genome = sarek.get("genome", "GATK.GRCh38")
        self.profile = sarek.get("profile", "docker")
        self.extra_args = sarek.get("extra_args", "")
        self._max_cpus = int(sarek.get("max_cpus", 16))
        self._max_memory = sarek.get("max_memory", "32.GB")

    @staticmethod
    def is_available() -> bool:
        """Check whether Nextflow is on PATH."""
        return shutil.which("nextflow") is not None

    def run(
        self,
        tumor_fastq_r1: str | Path,
        tumor_fastq_r2: str | Path,
        normal_fastq_r1: str | Path,
        normal_fastq_r2: str | Path,
        out_dir: str | Path,
        *,
        patient_id: str = "patient1",
        tumor_sample: str = "tumor",
        normal_sample: str = "normal",
    ) -> Optional[Path]:
        """Run Sarek and return path to the annotated VCF.

        Returns None if Nextflow is not installed or the pipeline fails.
        """
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        if not self.is_available():
            logger.error(
                "Nextflow not found on PATH. Install: "
                "curl -s https://get.nextflow.io | bash"
            )
            return None

        samplesheet = self._write_samplesheet(
            out_dir,
            patient_id=patient_id,
            tumor_sample=tumor_sample,
            normal_sample=normal_sample,
            tumor_r1=str(tumor_fastq_r1),
            tumor_r2=str(tumor_fastq_r2),
            normal_r1=str(normal_fastq_r1),
            normal_r2=str(normal_fastq_r2),
        )

        cmd = self._build_command(samplesheet, out_dir)
        logger.info("Launching Sarek: %s", " ".join(cmd))

        try:
            result = subprocess.run(
                cmd,
                cwd=str(out_dir),
                capture_output=True,
                text=True,
                timeout=72 * 3600,
            )
            if result.returncode != 0:
                logger.error("Sarek failed (exit %d):\n%s", result.returncode, result.stderr[-2000:])
                return None
        except subprocess.TimeoutExpired:
            logger.error("Sarek timed out after 72 hours.")
            return None
        except Exception as exc:
            logger.error("Sarek execution error: %s", exc)
            return None

        vcf = self._find_output_vcf(out_dir, tumor_sample)
        if vcf:
            logger.info("Sarek output VCF: %s", vcf)
        else:
            logger.warning("Could not locate annotated VCF in Sarek output.")
        return vcf

    def _write_samplesheet(
        self,
        out_dir: Path,
        **kwargs,
    ) -> Path:
        """Write the nf-core/sarek input CSV samplesheet."""
        ss_path = out_dir / "samplesheet.csv"
        with open(ss_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "patient", "sample", "lane", "fastq_1", "fastq_2", "status",
            ])
            writer.writerow([
                kwargs["patient_id"],
                kwargs["normal_sample"],
                "lane1",
                kwargs["normal_r1"],
                kwargs["normal_r2"],
                "0",
            ])
            writer.writerow([
                kwargs["patient_id"],
                kwargs["tumor_sample"],
                "lane1",
                kwargs["tumor_r1"],
                kwargs["tumor_r2"],
                "1",
            ])
        return ss_path

    def _build_command(self, samplesheet: Path, out_dir: Path) -> list[str]:
        """Construct the nextflow run command."""
        cmd = [
            "nextflow", "run",
            f"nf-core/sarek",
            "-r", self.version,
            "-profile", self.profile,
            "--input", str(samplesheet),
            "--outdir", str(out_dir / "sarek_results"),
            "--genome", self.genome,
            "--tools", "mutect2",
            "--annotate_tools", "vep",
            "--max_cpus", str(self._max_cpus),
            "--max_memory", self._max_memory,
        ]
        if self.extra_args:
            cmd.extend(self.extra_args.split())
        return cmd

    @staticmethod
    def _find_output_vcf(out_dir: Path, tumor_sample: str) -> Optional[Path]:
        """Locate the VEP-annotated VCF in Sarek's output tree."""
        patterns = [
            f"sarek_results/annotation/{tumor_sample}/*mutect2*vep*.vcf.gz",
            f"sarek_results/annotation/*/{tumor_sample}*mutect2*vep*.vcf.gz",
            "sarek_results/annotation/**/*mutect2*vep*.vcf.gz",
        ]
        for pat in patterns:
            matches = sorted(out_dir.glob(pat))
            if matches:
                return matches[0]
        return None
