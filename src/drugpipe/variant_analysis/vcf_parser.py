"""Parse VEP/SnpEff-annotated VCF files to extract actionable somatic mutations.

Extracts missense (and other protein-altering) variants from a filtered,
annotated VCF produced by GATK Mutect2 + Ensembl VEP.  Each variant is
mapped to its gene, transcript, and protein-level consequence (HGVSp)
so downstream modules can build mutant protein sequences.

Dependencies:
  - pysam (conda install -c bioconda pysam)  — optional, used for indexed VCFs
  - cyvcf2 (pip install cyvcf2)              — optional, faster VCF parsing
  Falls back to a lightweight built-in parser when neither is installed.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

logger = logging.getLogger(__name__)

# VEP CSQ sub-field order (Ensembl default from --vcf output)
_DEFAULT_CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
    "SYMBOL_SOURCE", "HGNC_ID",
]

PROTEIN_ALTERING = frozenset({
    "missense_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
})

HIGH_MODERATE_IMPACT = frozenset({"HIGH", "MODERATE"})


@dataclass
class SomaticVariant:
    """A single protein-altering somatic variant."""

    chrom: str
    pos: int
    ref: str
    alt: str
    gene_symbol: str
    transcript_id: str
    hgvsp: str
    consequence: str
    impact: str
    protein_position: Optional[int] = None
    ref_aa: Optional[str] = None
    alt_aa: Optional[str] = None
    extra: Dict[str, str] = field(default_factory=dict)

    @property
    def short_label(self) -> str:
        """Human-readable label like 'EGFR L858R'."""
        if self.ref_aa and self.alt_aa and self.protein_position:
            return f"{self.gene_symbol} {self.ref_aa}{self.protein_position}{self.alt_aa}"
        return f"{self.gene_symbol} {self.hgvsp}"


class VCFParser:
    """Extract protein-altering somatic variants from a VEP-annotated VCF."""

    def __init__(
        self,
        cfg: Optional[Dict[str, Any]] = None,
        *,
        driver_genes: Optional[Sequence[str]] = None,
        min_impact: str = "MODERATE",
    ):
        va = (cfg or {}).get("variant_analysis", {})
        self.driver_genes: Optional[frozenset[str]] = None
        gene_list = va.get("driver_genes", None) or driver_genes
        if gene_list:
            self.driver_genes = frozenset(g.upper() for g in gene_list)

        self.min_impact = min_impact.upper()
        self._accepted_impacts = {"HIGH"}
        if self.min_impact == "MODERATE":
            self._accepted_impacts.add("MODERATE")
        elif self.min_impact == "LOW":
            self._accepted_impacts |= {"MODERATE", "LOW"}

    def parse(self, vcf_path: str | Path) -> List[SomaticVariant]:
        """Parse a VCF file and return protein-altering variants."""
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")

        variants = self._parse_vcf_text(vcf_path)
        logger.info(
            "Parsed %d protein-altering variants from %s", len(variants), vcf_path.name,
        )
        return variants

    def _parse_vcf_text(self, vcf_path: Path) -> List[SomaticVariant]:
        """Lightweight VCF parser that extracts VEP CSQ annotations."""
        csq_fields: Optional[List[str]] = None
        variants: List[SomaticVariant] = []

        open_fn = self._open_fn(vcf_path)
        with open_fn(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = self._parse_csq_header(line)
                    continue
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t", 8)
                if len(parts) < 8:
                    continue

                chrom, pos_str, _, ref, alt, _, filt, info = parts[:8]
                if filt not in ("PASS", ".", ""):
                    continue

                csq_str = self._extract_info_field(info, "CSQ")
                if not csq_str:
                    continue

                fields = csq_fields or _DEFAULT_CSQ_FIELDS
                for transcript_csq in csq_str.split(","):
                    var = self._parse_one_csq(
                        chrom, int(pos_str), ref, alt, transcript_csq, fields,
                    )
                    if var is not None:
                        variants.append(var)

        return variants

    def _parse_one_csq(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        csq_str: str,
        fields: List[str],
    ) -> Optional[SomaticVariant]:
        """Parse a single VEP CSQ annotation entry."""
        vals = csq_str.split("|")
        csq = {fields[i]: vals[i] if i < len(vals) else "" for i in range(len(fields))}

        impact = csq.get("IMPACT", "")
        if impact not in self._accepted_impacts:
            return None

        consequence = csq.get("Consequence", "")
        if not any(c in PROTEIN_ALTERING for c in consequence.split("&")):
            return None

        gene = csq.get("SYMBOL", "")
        if self.driver_genes and gene.upper() not in self.driver_genes:
            return None

        hgvsp = csq.get("HGVSp", "")
        amino_acids = csq.get("Amino_acids", "")
        prot_pos_str = csq.get("Protein_position", "")

        ref_aa, alt_aa = None, None
        if "/" in amino_acids:
            parts = amino_acids.split("/")
            ref_aa, alt_aa = parts[0], parts[1] if len(parts) > 1 else None

        prot_pos: Optional[int] = None
        if prot_pos_str and prot_pos_str != "-":
            try:
                prot_pos = int(prot_pos_str.split("-")[0])
            except ValueError:
                pass

        return SomaticVariant(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            gene_symbol=gene,
            transcript_id=csq.get("Feature", ""),
            hgvsp=hgvsp,
            consequence=consequence,
            impact=impact,
            protein_position=prot_pos,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            extra={k: v for k, v in csq.items() if v and k not in {
                "SYMBOL", "Feature", "HGVSp", "Consequence", "IMPACT",
                "Amino_acids", "Protein_position",
            }},
        )

    @staticmethod
    def _parse_csq_header(header_line: str) -> List[str]:
        """Extract field names from the VEP CSQ meta-info line."""
        match = re.search(r'Format:\s*([^"]+)', header_line)
        if match:
            return match.group(1).strip().rstrip(">").rstrip('"').split("|")
        return _DEFAULT_CSQ_FIELDS

    @staticmethod
    def _extract_info_field(info: str, key: str) -> Optional[str]:
        """Extract a specific key=value from the VCF INFO column."""
        prefix = f"{key}="
        for part in info.split(";"):
            if part.startswith(prefix):
                return part[len(prefix):]
        return None

    @staticmethod
    def _open_fn(path: Path):
        """Return gzip.open for .gz files, otherwise builtins.open."""
        if path.suffix == ".gz":
            import gzip
            return gzip.open
        return open
