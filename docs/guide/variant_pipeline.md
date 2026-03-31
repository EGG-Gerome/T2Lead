# Variant-driven workflow (core path)

## Overview

1. Obtain **tumor/normal** alignments and somatic variants (e.g. **nf-core/sarek** with Mutect2 + VEP), or provide a **ready-made VEP-annotated VCF**.
2. Enable **`variant_analysis`** in config and set `vcf_path` **or** all four FASTQ paths.
3. Prepare (or reuse) **`stage3_leads/final_lead_candidates.csv`** from Stages 2–3 so Stage 4 has molecules to dock.
4. Run the pipeline: it parses VCF → builds mutant sequences → resolves structures (PDB / ESMFold) → runs **Stage 4** per variant into `stage4_optimization/<GENE_MUTATION>/`.

## Configuration

In `configs/default_config.yaml` (or a custom YAML):

```yaml
pipeline:
  stages:
    - lead_optimization    # often: only Stage 4 after upstream genomics

variant_analysis:
  enabled: true
  vcf_path: "data/user_inputs/vcf/sample.ann.vcf"   # or leave empty and set FASTQs
  # tumor_fastq_r1, tumor_fastq_r2, normal_fastq_r1, normal_fastq_r2
  driver_genes: []        # empty = all; or ["EGFR", "PIK3CA", ...]
  min_impact: "MODERATE"
```

**Automatic handoff:** when `vcf_path` is set, the Python pipeline parses it directly. When FASTQs are set and `vcf_path` is empty, `SarekRunner` invokes Nextflow; the discovered VCF path is then parsed — no manual copy to a fixed `results/sarek/...` path unless you run sarek standalone.

## Standalone sarek (Makefile)

```bash
make install-sarek
make run-sarek INPUT=data/user_inputs/samplesheet.csv
```

Point `variant_analysis.vcf_path` at the **annotated** VCF emitted by your Nextflow `outdir` (see sarek docs for exact filenames).

## Prerequisites for Stage 4 only

You **must** have lead molecules on disk:

`<out_dir>/<disease_slug>/stage3_leads/final_lead_candidates.csv`

(Or run full Stages 2–3 once for the same disease/config.)

## Outputs per mutation

Under `stage4_optimization/<GENE_MUTATION>/`:

- `optimized_leads.csv`
- Receptor PDB/PDBQT, `docking_poses/`, optional `md_trajectories/`

The mutant structure is copied into `lead_optimization`’s working directory via `lead_optimization._mutant_pdb_path` (set by the orchestrator).

## References

- [somatic_variant_pipeline_feasibility.md](../research/somatic_variant_pipeline_feasibility.md)
- [reproduction_steps.md](../reproducibility/reproduction_steps.md)
