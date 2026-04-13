"""Tests for Xena TSV → minimal VEP VCF conversion."""

from __future__ import annotations

from pathlib import Path

import pytest

from drugpipe.variant_analysis.vcf_parser import VCFParser
from drugpipe.variant_analysis.xena_adapter import (
    convert_xena_tsv_to_vep_vcf,
    is_xena_somatic_mutation_tsv,
    pick_xena_sample,
)


def test_is_xena_tsv_header(tmp_path: Path) -> None:
    p = tmp_path / "m.tsv"
    p.write_text(
        "sample\tgene\tchrom\tstart\tend\tref\talt\t"
        "Tumor_Sample_Barcode\tAmino_Acid_Change\teffect\tcallers\tdna_vaf\n",
        encoding="utf-8",
    )
    assert is_xena_somatic_mutation_tsv(p) is True


def test_convert_subset_and_parse(tmp_path: Path) -> None:
    tsv = tmp_path / "tiny.tsv"
    tsv.write_text(
        "\t".join([
            "sample", "gene", "chrom", "start", "end", "ref", "alt",
            "Tumor_Sample_Barcode", "Amino_Acid_Change", "effect", "callers", "dna_vaf",
        ])
        + "\n"
        + "\t".join([
            "S1", "PIK3CA", "chr3", "179234297", "179234297", "A", "G",
            "bar", "p.H1047R", "missense_variant", "m", "0.5",
        ])
        + "\n",
        encoding="utf-8",
    )
    cfg = {
        "variant_analysis": {
            "xena": {"sample_strategy": "explicit", "sample_id": "S1"},
        },
    }
    out = convert_xena_tsv_to_vep_vcf(tsv, tmp_path / "out", cfg)
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert "CSQ=" in text
    parser = VCFParser(cfg)
    variants = parser.parse(out)
    assert len(variants) == 1
    assert "PIK3CA" in variants[0].short_label


def test_pick_max_variants() -> None:
    import pandas as pd

    df = pd.DataFrame({
        "sample": ["a", "a", "b", "b", "b", "b"],
    })
    assert pick_xena_sample(df, "max_variants", "") == "b"


def test_resolve_md_gpu_device_indices_explicit() -> None:
    from drugpipe.lead_optimization.md_simulation import _parse_int_list, _resolve_md_gpu_device_indices

    assert _parse_int_list([0, 1, 2]) == [0, 1, 2]
    assert _parse_int_list("0, 2 ,4") == [0, 2, 4]
    md = {"gpu_device_indices": [1, 3], "auto_discover_gpu_devices": True}
    assert _resolve_md_gpu_device_indices(md) == [1, 3]
