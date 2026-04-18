"""Convert UCSC Xena / TCGA tabular mutation files to minimal VEP-style VCF.

The pipeline expects VCF lines with a ``CSQ=`` INFO field compatible with
:class:`drugpipe.variant_analysis.vcf_parser.VCFParser`.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from drugpipe.variant_analysis import mutant_sequence as mutseq
from drugpipe.variant_analysis.vcf_parser import PROTEIN_ALTERING

logger = logging.getLogger(__name__)

# Typical Xena TCGA mutation table (e.g. *somaticmutation_wxs.tsv)
_XENA_REQUIRED_COLS = frozenset({
    "sample", "gene", "chrom", "start", "end", "ref", "alt",
    "Amino_Acid_Change", "effect",
})


def is_xena_somatic_mutation_tsv(path: Path) -> bool:
    """Return True if *path* looks like a Xena somatic mutation TSV."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline()
    except OSError:
        return False
    if "\t" not in header:
        return False
    cols = set(header.strip().split("\t"))
    return _XENA_REQUIRED_COLS.issubset(cols)


def _impact_for_consequence(consequence: str) -> str:
    c = consequence.lower()
    if c in ("stop_gained", "stop_lost", "start_lost", "frameshift_variant"):
        return "HIGH"
    return "MODERATE"


def _pick_primary_consequence(effect_field: str) -> Optional[str]:
    """Pick first VEP-like consequence token that is protein-altering."""
    if not effect_field or not str(effect_field).strip():
        return None
    for tok in str(effect_field).replace("|", ";").split(";"):
        t = tok.strip()
        if not t:
            continue
        if t in PROTEIN_ALTERING:
            return t
    return None


# p.H1047R  /  p.His1047Arg  /  p.G12D  /  p.Q173*
_AA_CHG_3 = re.compile(
    r"^p\.(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|\*|\?)$",
)
_AA_CHG_1 = re.compile(
    r"^p\.(?P<ref>[A-Z\*])(?P<pos>\d+)(?P<alt>[A-Z\*])(?:/|$)?",
)


def _aa3_to_1(aa3: str) -> Optional[str]:
    if len(aa3) == 1:
        return aa3.upper()
    key = aa3[:1].upper() + aa3[1:].lower()
    return mutseq._AA3TO1.get(key) if aa3 else None  # noqa: SLF001


def _parse_amino_acid_change(raw: str) -> Optional[Tuple[str, int, str]]:
    """Return (ref_aa, protein_position, alt_aa) single-letter, or None."""
    s = (raw or "").strip()
    if not s or s in (".", "-"):
        return None
    m = _AA_CHG_1.match(s)
    if m:
        ref_a = m.group("ref").upper()
        alt_a = m.group("alt").upper()
        pos = int(m.group("pos"))
        if ref_a == "*":
            ref_a = "*"
        return ref_a, pos, alt_a
    m3 = _AA_CHG_3.match(s)
    if m3:
        ref_a = _aa3_to_1(m3.group("ref"))
        alt_s = m3.group("alt")
        alt_a = _aa3_to_1(alt_s) if len(alt_s) == 3 else alt_s.upper()
        pos = int(m3.group("pos"))
        if ref_a and alt_a:
            return ref_a, pos, alt_a
    return None


def _normalize_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


def _build_csq(
    alt_allele: str,
    consequence: str,
    gene: str,
    hgvsp_short: str,
    ref_aa: str,
    alt_aa: str,
    protein_position: int,
) -> str:
    """One CSQ transcript block (23 fields; matches ``VCFParser`` default order)."""
    sym = gene.strip() or "UNKNOWN"
    impact = _impact_for_consequence(consequence)
    aa_field = f"{ref_aa}/{alt_aa}"
    fields = [
        alt_allele,
        consequence,
        impact,
        sym,
        "",
        "Transcript",
        "",
        "protein_coding",
        "",
        "",
        "",
        hgvsp_short,
        "",
        "",
        str(protein_position),
        aa_field,
        "",
        "",
        "",
        "",
        "",
        "HGNC",
        "",
    ]
    if len(fields) != 23:
        raise RuntimeError(f"CSQ field count mismatch: {len(fields)}")
    return "|".join(fields)


def pick_xena_sample(
    df: pd.DataFrame,
    strategy: str,
    explicit_id: str,
) -> str:
    """Resolve which ``sample`` column value to use."""
    eid = (explicit_id or "").strip()
    if eid:
        if eid not in set(df["sample"].astype(str)):
            raise ValueError(
                f"xena sample_id={eid!r} not found in TSV sample column.",
            )
        return eid
    st = (strategy or "max_variants").strip().lower()
    if st == "max_variants":
        vc = df["sample"].astype(str).value_counts()
        chosen = str(vc.index[0])
        logger.info(
            "Xena: selected sample %s with %d mutation rows (max_variants).",
            chosen,
            int(vc.iloc[0]),
        )
        return chosen
    raise ValueError(f"Unknown variant_analysis.xena.sample_strategy: {strategy!r}")


def convert_xena_tsv_to_vep_vcf(
    tsv_path: Path,
    out_dir: Path,
    cfg: Dict[str, Any],
) -> Path:
    """Read Xena TSV, subset to one sample, emit minimal VCF with CSQ.

    Parameters
    ----------
    tsv_path
        Path to ``*.tsv`` (must pass :func:`is_xena_somatic_mutation_tsv`).
    out_dir
        Directory for ``xena_converted.<sample>.vep.vcf``.
    cfg
        Full pipeline config; reads ``variant_analysis.xena``.
    """
    va = cfg.get("variant_analysis", {}) or {}
    xa = va.get("xena", {}) or {}
    strategy = str(xa.get("sample_strategy", "max_variants"))
    explicit = str(xa.get("sample_id", "") or "").strip()

    tsv_path = Path(tsv_path).resolve()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    missing = _XENA_REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Xena TSV missing columns: {sorted(missing)}")

    sample_id = pick_xena_sample(df, strategy, explicit)
    sub = df[df["sample"].astype(str) == sample_id].copy()
    logger.info(
        "Xena: %d rows for sample %s (from %s).",
        len(sub),
        sample_id,
        tsv_path.name,
    )

    lines: List[str] = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations (minimal VEP-like). '
        'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|'
        'HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|'
        'Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">',
        '##INFO=<ID=XENA_SOURCE,Number=1,Type=String,Description="Original Xena TSV path">',
        '##INFO=<ID=XENA_SAMPLE,Number=1,Type=String,Description="Xena sample id">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]

    src_tag = str(tsv_path)
    n_out = 0
    skipped = 0

    for _, row in sub.iterrows():
        cons = _pick_primary_consequence(str(row.get("effect", "")))
        if not cons:
            skipped += 1
            continue

        aa_raw = row.get("Amino_Acid_Change", "")
        parsed = _parse_amino_acid_change(str(aa_raw) if pd.notna(aa_raw) else "")
        if not parsed:
            skipped += 1
            continue
        ref_aa, prot_pos, alt_aa = parsed

        ref = str(row.get("ref", "")).strip()
        alt = str(row.get("alt", "")).strip()
        if not ref and not alt:
            skipped += 1
            continue

        chrom = _normalize_chrom(str(row.get("chrom", "")))
        try:
            pos = int(row["start"])
        except (TypeError, ValueError):
            skipped += 1
            continue

        gene = str(row.get("gene", "")).strip()
        alt_allele = alt if alt not in (".", "-", "") else ""
        if not alt_allele:
            # deletion-style: VEP CSQ Allele often still encodes alt base or '-'
            alt_allele = alt if alt else "-"

        hgvsp_short = f"p.{ref_aa}{prot_pos}{alt_aa}"
        csq = _build_csq(
            alt_allele,
            cons,
            gene,
            hgvsp_short,
            ref_aa,
            alt_aa,
            prot_pos,
        )

        info = (
            f"CSQ={csq};"
            f"XENA_SOURCE={_vcf_info_escape(src_tag)};"
            f"XENA_SAMPLE={_vcf_info_escape(sample_id)}"
        )

        # Minimal VCF: REF/ALT must be non-empty for SNVs; for indels keep as in table
        if not ref:
            ref = "."

        line = "\t".join([
            chrom,
            str(pos),
            ".",
            ref,
            alt if alt not in ("", ".") else ".",
            ".",
            "PASS",
            info,
        ])
        lines.append(line)
        n_out += 1

    out_path = out_dir / f"xena_converted.{_safe_slug(sample_id)}.vep.vcf"
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info(
        "Wrote Xena-derived VCF: %s (%d variants, %d rows skipped).",
        out_path,
        n_out,
        skipped,
    )
    return out_path


def _vcf_info_escape(s: str) -> str:
    return str(s).replace("%", "%25").replace(";", "%3B").replace(" ", "_")


def _safe_slug(sample_id: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", sample_id)[:200] or "sample"
