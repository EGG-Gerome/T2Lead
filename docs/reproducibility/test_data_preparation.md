# Test data preparation (WES / WGS)

T2Lead’s `variant_analysis` module expects either:

1. **Paired FASTQ** for tumor and normal samples (passed through to nf-core/sarek when `vcf_path` is empty), or  
2. A **VEP-annotated VCF** with consequence fields compatible with `vcf_parser.py`.

## Public references

- **1000 Genomes / gnomAD / TCGA** — download BAM/FASTQ per journal / dbGaP policies; convert to FASTQ if you only have BAM.
- **Synthetic / minimal VCF** — for CI-style checks, craft a small VCF with one missense variant and PASS filters; see `tests/test_variant_analysis.py` for VEP-like INFO patterns.

## Samplesheet (sarek)

nf-core/sarek expects CSV columns such as `patient,sex,status,sample,lane,fastq_1,fastq_2` (see sarek 3.x docs). Place the CSV path in `make run-sarek INPUT=...`.

## Consent and scope

Use de-identified research data; exome panels should match the **genome build** configured for sarek (`GRCh38` vs `GRCh37`) and for VEP annotations.
