"""Tests for the variant_analysis module: VCF parsing, mutant sequences, and structure bridge."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from drugpipe.variant_analysis.vcf_parser import (
    VCFParser,
    SomaticVariant,
    PROTEIN_ALTERING,
)
from drugpipe.variant_analysis.mutant_sequence import MutantSequenceBuilder, _AA3TO1


# ======================================================================
# Fixtures
# ======================================================================

VEP_VCF_CONTENT = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    7\t55259515\t.\tT\tG\t.\tPASS\tCSQ=G|missense_variant|MODERATE|EGFR|ENSG00000146648|Transcript|ENST00000275493|protein_coding|21/28||ENST00000275493.6:c.2573T>G|ENSP00000275493.2:p.Leu858Arg|2695|2573|858|L/R|cTg/cGg|COSM6224||1||HGNC|HGNC:3236
    3\t178936091\t.\tA\tG\t.\tPASS\tCSQ=G|missense_variant|MODERATE|PIK3CA|ENSG00000121879|Transcript|ENST00000263967|protein_coding|21/21||ENST00000263967.3:c.3140A>G|ENSP00000263967.3:p.His1047Arg|3207|3140|1047|H/R|cAt/cGt|COSM775|||1||HGNC|HGNC:8975
    1\t100000\t.\tA\tT\t.\tlow_qual\tCSQ=T|synonymous_variant|LOW|FAKEGENE|ENG1|Transcript|ENST999|protein_coding|1/1|||||||||||HGNC|HGNC:0
""")


@pytest.fixture
def vcf_file(tmp_path: Path) -> Path:
    vcf = tmp_path / "test.vcf"
    vcf.write_text(VEP_VCF_CONTENT)
    return vcf


# ======================================================================
# VCF Parser Tests
# ======================================================================

class TestVCFParser:
    def test_parse_extracts_missense(self, vcf_file: Path):
        parser = VCFParser()
        variants = parser.parse(vcf_file)
        assert len(variants) == 2
        genes = {v.gene_symbol for v in variants}
        assert genes == {"EGFR", "PIK3CA"}

    def test_parse_egfr_details(self, vcf_file: Path):
        parser = VCFParser()
        variants = parser.parse(vcf_file)
        egfr = next(v for v in variants if v.gene_symbol == "EGFR")
        assert egfr.protein_position == 858
        assert egfr.ref_aa == "L"
        assert egfr.alt_aa == "R"
        assert egfr.consequence == "missense_variant"
        assert egfr.impact == "MODERATE"
        assert egfr.short_label == "EGFR L858R"

    def test_parse_pik3ca_details(self, vcf_file: Path):
        parser = VCFParser()
        variants = parser.parse(vcf_file)
        pik3ca = next(v for v in variants if v.gene_symbol == "PIK3CA")
        assert pik3ca.protein_position == 1047
        assert pik3ca.ref_aa == "H"
        assert pik3ca.alt_aa == "R"
        assert pik3ca.short_label == "PIK3CA H1047R"

    def test_filters_non_pass(self, vcf_file: Path):
        """low_qual filter should be excluded."""
        parser = VCFParser()
        variants = parser.parse(vcf_file)
        assert all(v.gene_symbol != "FAKEGENE" for v in variants)

    def test_driver_gene_filter(self, vcf_file: Path):
        parser = VCFParser(driver_genes=["EGFR"])
        variants = parser.parse(vcf_file)
        assert len(variants) == 1
        assert variants[0].gene_symbol == "EGFR"

    def test_high_impact_only(self, vcf_file: Path):
        parser = VCFParser(min_impact="HIGH")
        variants = parser.parse(vcf_file)
        assert len(variants) == 0

    def test_missing_file(self):
        parser = VCFParser()
        with pytest.raises(FileNotFoundError):
            parser.parse("/nonexistent/file.vcf")


# ======================================================================
# Amino Acid Conversion Tests
# ======================================================================

class TestAminoAcidConversion:
    def test_three_letter_to_one(self):
        builder = MutantSequenceBuilder()
        assert builder._to_single_letter("Leu") == "L"
        assert builder._to_single_letter("Arg") == "R"
        assert builder._to_single_letter("His") == "H"

    def test_single_letter_passthrough(self):
        builder = MutantSequenceBuilder()
        assert builder._to_single_letter("L") == "L"
        assert builder._to_single_letter("R") == "R"

    def test_none_input(self):
        builder = MutantSequenceBuilder()
        assert builder._to_single_letter(None) is None

    def test_all_standard_aa(self):
        for three, one in _AA3TO1.items():
            if three == "Ter":
                continue
            assert MutantSequenceBuilder._to_single_letter(three) == one


# ======================================================================
# Mutant Sequence Builder (unit tests, no network)
# ======================================================================

class TestMutantSequenceBuilder:
    def test_build_with_mock_sequence(self, vcf_file: Path, monkeypatch):
        """Test sequence building with a mocked wildtype fetch."""
        parser = VCFParser()
        variants = parser.parse(vcf_file)
        builder = MutantSequenceBuilder()

        fake_seq = "M" + "A" * 856 + "L" + "A" * 200
        monkeypatch.setattr(builder, "_fetch_sequence", lambda *a: fake_seq)
        monkeypatch.setattr(builder, "_resolve_uniprot", lambda *a: "P00533")

        egfr_var = next(v for v in variants if v.gene_symbol == "EGFR")
        results = builder.build([egfr_var])
        assert len(results) == 1
        mp = results[0]
        assert mp.mutant_sequence[857] == "R"
        assert mp.wildtype_sequence[857] == "L"
        assert mp.mutation_label == "L858R"

    def test_write_fasta(self, tmp_path: Path):
        from drugpipe.variant_analysis.mutant_sequence import MutantProtein
        mp = MutantProtein(
            variant=SomaticVariant(
                chrom="7", pos=55259515, ref="T", alt="G",
                gene_symbol="EGFR", transcript_id="ENST1",
                hgvsp="p.Leu858Arg", consequence="missense_variant",
                impact="MODERATE", protein_position=858,
                ref_aa="L", alt_aa="R",
            ),
            gene_symbol="EGFR",
            uniprot_id="P00533",
            wildtype_sequence="MAAAL" + "A" * 95,
            mutant_sequence="MAAAR" + "A" * 95,
            mutation_label="L858R",
        )
        builder = MutantSequenceBuilder()
        paths = builder.write_fastas([mp], tmp_path / "fastas")
        assert len(paths) == 1
        assert paths[0].exists()
        content = paths[0].read_text()
        assert ">EGFR_L858R" in content
        assert "MAAAR" in content
