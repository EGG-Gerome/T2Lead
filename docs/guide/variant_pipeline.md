# Variant-driven workflow (core path)

[中文版 (variant_pipeline_zh.md)](variant_pipeline_zh.md)

## Overview

1. Obtain **tumor/normal** alignments and somatic variants (e.g. **nf-core/sarek** with Mutect2 + VEP), or provide a **ready-made VEP-annotated VCF**.
2. Enable **`variant_analysis`** in config and set `vcf_path` **or** all four FASTQ paths.
3. Run the pipeline: it parses VCF → builds mutant sequences → resolves structures (PDB / ESMFold) → runs **Stage 4** per variant.

## Lead molecule source (choose one)

| Mode | Config | Description |
|------|--------|-------------|
| **Auto Stage 2-3** (recommended) | `auto_stage23: true` | Each variant gene is automatically mapped to a ChEMBL target → screening → lead generation → docking against mutant structure. **No pre-existing CSV needed.** |
| **Reuse existing leads** | `auto_stage23: false` | Uses `stage3_leads/final_lead_candidates.csv` from a prior disease-path run. Fast but molecules may not be specific to the mutant target. |

## Configuration

In `configs/default_config.yaml` (or a custom YAML like `configs/variant_breast.yaml`):

```yaml
pipeline:
  stages:
    - lead_optimization

variant_analysis:
  enabled: true
  vcf_path: "data/user_inputs/vcf/sample.ann.vcf"   # or leave empty and set FASTQs
  sample_id: "patient_001"    # optional but recommended for multi-patient isolation
  driver_genes: []        # empty = all; or ["EGFR", "PIK3CA", ...]
  min_impact: "MODERATE"
  auto_stage23: true      # auto-run Stage 2-3 per variant gene
```

## CLI usage (no file editing needed)

```bash
# Provide VCF + auto Stage 2-3
python scripts/run_pipeline.py \
  --disease "breast cancer" \
  --vcf-path /path/to/your.vcf.gz \
  --auto-stage23 \
  --stages lead_optimization \
  -v

# Or just use a config file
python scripts/run_pipeline.py -c configs/variant_breast.yaml -v
```

## One-command main path runner

For multi-patient isolated runs without long CLI flags, use:

```bash
python scripts/run_mainpath.py \
  --disease "breast cancer" \
  --sample-id patient_001 \
  --vcf-path /path/to/patient_001.ann.vcf.gz \
  -v
```

Defaults in `run_mainpath.py`:
- output root: `/root/autodl-fs/T2Lead_mainpath`
- `auto_stage23: true`
- Sarek profile: `singularity`
- ChEMBL shared cache disabled (isolation first)

Available CLI flags:

| Flag | Description |
|------|-------------|
| `--disease "name"` | Set disease name (determines output subdirectory) |
| `--vcf-path /path` | Auto-enable variant_analysis and set VCF path |
| `--auto-stage23` | Enable per-gene automatic Stage 2-3 |
| `--stages s1 s2 ...` | Override pipeline stages |

## Automatic handoff

When `vcf_path` is set, the Python pipeline parses it directly. When FASTQs are set and `vcf_path` is empty, `SarekRunner` invokes Nextflow; the discovered VCF path is then parsed — no manual copy to a fixed `results/sarek/...` path unless you run sarek standalone.

Important:
- **Samplesheet is only for FASTQ path (Sarek)**.
- **VCF path does not need a samplesheet**; it points directly to one annotated VCF file.
- After VEP annotation, it is still a VCF (`.vcf` / `.vcf.gz`) with `CSQ` in INFO.

## Standalone sarek (Makefile)

```bash
make install-sarek
make run-sarek INPUT=data/user_inputs/samplesheet.csv
```

Point `variant_analysis.vcf_path` at the **annotated** VCF emitted by your Nextflow `outdir` (see sarek docs for exact filenames).

## Outputs per mutation

By default (`variant_isolated_runs: true`) outputs are isolated under:

- `<pipeline.out_dir>/<disease_slug>/variant_runs/<sample_id>/<run_id>/...`

Then per mutation, under `stage4_optimization/<GENE_MUTATION>/`:

- `optimized_leads.csv`
- Receptor PDB/PDBQT, `docking_poses/`, optional `md_trajectories/`

When `auto_stage23: true`, each mutation directory also contains:

- `stage2_hits/` — ChEMBL screening results for that gene target
- `stage3_leads/` — lead candidates specific to that gene target

The mutant structure is copied into `lead_optimization`'s working directory via `lead_optimization._mutant_pdb_path` (set by the orchestrator).

## References

- [somatic_variant_pipeline_feasibility.md](../research/somatic_variant_pipeline_feasibility.md)
- [reproduction_steps.md](../reproducibility/reproduction_steps.md)
