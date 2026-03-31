# Reproduction checklist

Use this to verify T2Lead on a clean machine or for publication supplements.

## 1. Environment

- [ ] Install Conda; `make install` from repo root (or follow `Dockerfile`).
- [ ] `make test` → all tests pass (4 OpenMM-heavy tests may skip without OpenMM).

## 2. Minimal disease path

```bash
export DP_PIPELINE__OUT_DIR=/tmp/t2lead_smoke
python scripts/run_pipeline.py --disease "breast cancer" \
  --stages target_discovery target_to_hit \
  -v
```

Expect `stage1_targets/ranked_targets.csv` and `stage2_hits/` populated (ChEMBL crawl may take a long time; reduce limits in YAML for smoke tests).

## 3. Variant path (integration)

- [ ] Install Nextflow + Java 21; Docker running for sarek.
- [ ] Place tumor/normal FASTQs and a valid sarek **samplesheet** under `data/user_inputs/` (see nf-core/sarek docs).
- [ ] Run `make run-sarek INPUT=...` **or** provide a VEP VCF via `variant_analysis.vcf_path`.
- [ ] Ensure `stage3_leads/final_lead_candidates.csv` exists (run Stages 2–3 first, or copy a small test CSV with `canonical_smiles`).
- [ ] Set `variant_analysis.enabled: true` and run `python scripts/run_pipeline.py -c variant.yaml -v`.
- [ ] Confirm `stage4_optimization/<GENE_MUT>/optimized_leads.csv` appears.

## 4. Synthetic VCF / unit level

- [ ] `python -m pytest tests/test_variant_analysis.py -v` (no network).

## 5. Optional GPU

- [ ] After `make install` with CUDA PyTorch, run Stage 4 with `pipeline.device: cuda` and confirm OpenMM selects CUDA in logs.

## Troubleshooting

- **Empty optimized leads:** check ADMET hard filter removed all rows; relax `sa_score_max` or disable specific `hard_filter` flags temporarily.
- **Sarek / VCF path:** the Makefile writes to `variant_calling/sarek_results/` by default; point `vcf_path` at the real annotated VCF from that run.
- **Dual output roots:** always set `DP_PIPELINE__OUT_DIR` so intermediate and final artifacts share one tree.

For WES preparation (public data, BQSR, panel targeting), follow your institution’s best practices; T2Lead consumes **FASTQ or VCF**, not raw instrument files.
