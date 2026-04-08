# T2Lead

[中文说明 (README_zh.md)](README_zh.md)

**Primary path — Somatic variants → mutant protein → docking / MD:** pair tumor/normal exome data (or a VEP-annotated VCF) with a small-molecule lead library and screen compounds against **patient-specific** protein structures.

**Secondary path — Disease → targets → hits → leads:** classic modular pipeline from a disease string through ChEMBL, ML virtual screening, CReM / REINVENT4, and Stage 4 optimization.

```
Tumor/normal FASTQ or VCF ──► variant_analysis (VCF → mutations → structures)
                                        │
                    Lead CSV (Stage 3) ─► Stage 4: Vina + ADMET gates + MM-GBSA (OpenMM)

Disease name ──► Stages 1–3 (OpenTargets, ChEMBL, MPO / RL) ──► Stage 4
```

## Quick start

```bash
cd T2Lead
make install          # conda env + RDKit/OpenMM + REINVENT4 + .env (large download)
make test             # pytest (no GPU required)

# Secondary workflow (disease-driven)
make run DISEASE="breast cancer"

# Variant-driven Stage 4 (after you have final_lead_candidates.csv — see docs)
# Set variant_analysis.enabled: true and vcf_path or FASTQs in configs/default_config.yaml
python scripts/run_pipeline.py -c my_variant_config.yaml -v
```

Detailed setup: [docs/guide/quickstart.md](docs/guide/quickstart.md).

## CI and PR workflow

- GitHub Actions runs the `CI` workflow on every pull request and every branch push.
- The local equivalent is `make test` or `python -m pytest tests/ -v`.
- Prefer short-lived branches such as `feat/*`, `fix/*`, `docs/*`, and `chore/*`.
- Open pull requests into `main` and wait for required `CI` checks to pass before merge.
- Contributor conventions live in [CONTRIBUTING.md](CONTRIBUTING.md).

## Outputs (per disease / run)

With default `pipeline.output_layout.use_stage_subdirs`, artifacts live under `<pipeline.out_dir>/<disease_slug>/`:

| Directory | Contents |
|-----------|----------|
| `shared/chembl_stage2/` | Default shared ChEMBL molecule + IC50 crawl, `crawl_state.json`, `fp_cache/` (`target_to_hit.shared_library_dir`) |
| `stage1_targets/` | `ranked_targets.csv` |
| `stage2_hits/` | Per-disease training table, `model_cache/`, screening outputs, `final_hit_candidates.csv`, … |
| `stage3_leads/` | `final_lead_candidates.csv`, REINVENT4 files |
| `stage4_optimization/` | PDB/PDBQT, `docking_poses/`, `optimized_leads.csv`, MD trajectories |
| `logs/` | Full and summary logs |

User-owned inputs: [data/user_inputs/README.md](data/user_inputs/README.md).

## Scoring (Stage 4)

- **ADMET hard filter**: hERG liability, Veber failure, and high synthetic accessibility (SA) **drop** compounds before MD.
- **CYP risk**: soft weight only (`lead_optimization.scoring.w_cyp_soft`).
- **Composite rank**: AutoDock Vina (25%), implicit-solvent MM-GBSA-style term (40%), pose / trajectory RMSD stability (35%), plus optional CYP soft term — weights in [configs/default_config.yaml](configs/default_config.yaml).
- **Ensemble MM-GBSA**: short independent MD replicas sample complex energy; reported **mean ± std** of an approximate ΔG (see [docs/guide/configuration.md](docs/guide/configuration.md)).

*Scores rank computational hypotheses; they are not experimental potency.*

## nf-core / sarek (optional upstream)

WES/WGS variant calling is documented in [docs/guide/variant_pipeline.md](docs/guide/variant_pipeline.md). Makefile targets: `make install-sarek`, `make run-sarek INPUT=samplesheet.csv`. **Does not** conflict with `make test` (single project Makefile).

## Container runtime

Stages 1–4 run **without any container** — conda is all you need. Containers are only required for the optional nf-core/sarek upstream (variant calling).

### Docker — local dev / servers with Docker access

For: laptops, Docker Desktop, cloud VMs with sudo.

```bash
# T2Lead image
docker compose build t2lead
docker compose run --rm t2lead python -m pytest /app/tests -q

# sarek
nextflow run nf-core/sarek -profile docker ...
```

GPU: uncomment the `deploy` section in [docker-compose.yml](docker-compose.yml). Default image uses CPU PyTorch.

### Apptainer (Singularity) — AutoDL / HPC / restricted containers

For: AutoDL instances, university HPC clusters, any environment where `dockerd` cannot start. These environments typically lack the kernel network privileges Docker needs (iptables / virtual bridges). Apptainer requires no daemon and no network capabilities — the standard choice on HPC.

```bash
# Install Apptainer (Ubuntu 22.04)
wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer_1.4.5_amd64.deb
sudo apt install -y ./apptainer_1.4.5_amd64.deb

# sarek
nextflow run nf-core/sarek -profile singularity ...
```

> **Note:** nf-core maintains both `-profile docker` and `-profile singularity`; they are functionally identical, differing only in the container backend. If you are working on a cloud GPU server such as AutoDL, use the Singularity path.

## Documentation index

| Doc | Purpose | 中文 |
|-----|---------|------|
| [docs/guide/quickstart.md](docs/guide/quickstart.md) | Install, env vars, first run | [quickstart_zh.md](docs/guide/quickstart_zh.md) |
| [docs/guide/variant_pipeline.md](docs/guide/variant_pipeline.md) | VCF/FASTQ → structures → Stage 4 | [variant_pipeline_zh.md](docs/guide/variant_pipeline_zh.md) |
| [docs/guide/disease_pipeline.md](docs/guide/disease_pipeline.md) | Stages 1–4 from disease name | [disease_pipeline_zh.md](docs/guide/disease_pipeline_zh.md) |
| [docs/guide/configuration.md](docs/guide/configuration.md) | YAML / `DP_*` reference | [configuration_zh.md](docs/guide/configuration_zh.md) |
| [docs/guide/output_reference.md](docs/guide/output_reference.md) | File-by-file outputs | [output_reference_zh.md](docs/guide/output_reference_zh.md) |
| [docs/reproducibility/reproduction_steps.md](docs/reproducibility/reproduction_steps.md) | Reviewer / replicator checklist | [reproduction_steps_zh.md](docs/reproducibility/reproduction_steps_zh.md) |
| [docs/reproducibility/test_data_preparation.md](docs/reproducibility/test_data_preparation.md) | WES/WGS test data | [test_data_preparation_zh.md](docs/reproducibility/test_data_preparation_zh.md) |
| [docs/research/somatic_variant_pipeline_feasibility.md](docs/research/somatic_variant_pipeline_feasibility.md) | Feasibility note | [somatic_variant_pipeline_feasibility_zh.md](docs/research/somatic_variant_pipeline_feasibility_zh.md) |

## Project layout (short)

```
T2Lead/
├── configs/default_config.yaml
├── src/drugpipe/          # pipeline, stages, variant_analysis, paths.py
├── scripts/run_pipeline.py
├── tests/
└── docs/guide/ … docs/reproducibility/
```

## License

This repository is an **internal commercial project** of hupper and is provided under **All rights reserved** by default.

No copying, redistribution, modification, or public release is permitted without prior written authorization from hupper.
