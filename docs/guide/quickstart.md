# Quickstart

## Requirements

- Linux or macOS (Windows via WSL2).
- **Conda** (Miniconda / Mambaforge) for RDKit + OpenMM.
- **Java 21+** if you run nf-core/sarek via Nextflow.
- **Docker** (or Singularity) for sarek’s default profile.
- NVIDIA GPU recommended for MLP training and OpenMM MD.

## Install

From the repository root:

```bash
make install
```

This creates the `t2lead` conda environment, installs the scientific stack, editable `drugpipe`, PyTorch (default: CUDA 12.4 wheels), CReM DB download, REINVENT4 + prior, and writes `.env` with `DP_*` paths.

Override PyTorch index, e.g. CPU-only:

```bash
make install TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu
```

Non-Make install mirrors the same steps; see the historical README in git history or follow `Makefile` targets.

## Environment variables

- `DP_PIPELINE__OUT_DIR` — output base directory (recommended on large disks, e.g. AutoDL: `/autodl-fs/data/T2Lead`).
- `DP_*` — nested overrides; see [configuration.md](configuration.md).

## User data locations

Convention (paths are suggestions; point YAML/CLI to any path):

- `data/user_inputs/activities/` — IC50 CSV for custom training data.
- `data/user_inputs/screening_library/` — SMILES screening library.
- `data/user_inputs/fastq/` — tumor/normal FASTQ for variant calling.
- `data/user_inputs/vcf/` — annotated VCF for `variant_analysis.vcf_path`.

See [../../data/user_inputs/README.md](../../data/user_inputs/README.md).

## Smoke tests

```bash
make test
# or
python -m pytest tests/ -v
```

## Next steps

- Variant-centric workflow: [variant_pipeline.md](variant_pipeline.md)
- Disease-centric workflow: [disease_pipeline.md](disease_pipeline.md)
