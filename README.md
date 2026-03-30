# T2Lead

[中文说明 (README_zh.md)](README_zh.md)

**Target → Hit → Lead → Optimized Lead**: an end-to-end, modular drug discovery pipeline that takes a disease name as input and outputs ranked, experimentally-validated lead candidate molecules.

```
Disease name ──► Target Discovery ──► Target-to-Hit ──► Hit-to-Lead ──► Lead Optimization ──► Optimized Leads CSV
 (string)        (OpenTargets /       (ChEMBL crawl /   (Scaffold /      (Docking / ADMET /
                  OriGene)             ML + VS + ADMET)   CReM / REINVENT4) MD / MM-GBSA)
```

## Quick I/O Summary

| Stage | Input | Output |
|---|---|---|
| **1 — Target Discovery** | Disease name `str` (e.g. `"lung cancer"`) or EFO ID (e.g. `EFO_0001378`) | Ranked target list: `chembl_id, symbol, name, score, source` |
| **2 — Target to Hit** | ChEMBL target ID `str` (e.g. `"CHEMBL204"`) | `data/<disease>/final_hit_candidates.csv` |
| **3 — Hit to Lead** | Hit CSV from Stage 2 (must contain `canonical_smiles`, `pred_pIC50_ens`) | `data/<disease>/final_lead_candidates.csv` |
| **4 — Lead Optimization** | Lead CSV from Stage 3 + PDB ID for target protein | `data/<disease>/optimized_leads.csv` |

## Important: What Happens Under the Hood

This is **NOT** a pre-trained model you just run. Here's what each stage actually does:

- **Stage 1** queries the OpenTargets public API online (seconds).
- **Stage 2** crawls ChEMBL to download molecule + IC50 data first (can take hours depending on `MAX_MOLECULES` / `MAX_ACTIVITIES` settings), then trains Random Forest + MLP models **on the fly** using the crawled data for your specific target. No pre-built model — it learns from real bioactivity data each run.
- **Stage 3** analyzes scaffolds and clusters (instant), generates analogs via CReM if a fragment DB is provided, scores everything with MPO, and optionally runs REINVENT4 RL-based molecular generation.
- **Stage 4** prepares the target protein structure (PDB), docks all lead candidates with AutoDock Vina, runs enhanced ADMET profiling (SA score, hERG, CYP), and optionally performs MD simulation + MM-GBSA binding energy calculation with OpenMM (GPU-accelerated).

## Features

- **Stage 1 — Target Discovery**: Query OpenTargets GraphQL API by disease name or EFO ID; optionally merge results with an OriGene AI agent for LLM-assisted target recommendation.
- **Stage 2 — Target to Hit**: Crawl ChEMBL molecules + IC50 activities (with checkpoint resume), build pIC50 datasets, train Random Forest + Torch MLP regressors (GPU-accelerated), run virtual screening on the full molecule library, and filter by ADMET rules / QED / PAINS alerts.
- **Stage 3 — Hit to Lead**: Murcko scaffold SAR analysis, Butina diversity clustering, CReM-based analog enumeration, multi-parameter optimization (MPO) scoring, and optional REINVENT4 RL-based molecular generation.
- **Stage 4 — Lead Optimization**: Protein preparation (RCSB PDB fetch + PDBFixer), molecular docking (AutoDock Vina), enhanced ADMET profiling (SA score, hERG toxicophore, CYP inhibition risk, Veber rules), optional MD simulation + MM-GBSA (OpenMM with CUDA GPU support), and composite scoring.
- **GPU Support**: Automatically detects CUDA (NVIDIA), MPS (Apple Silicon), or falls back to CPU. Configurable via `pipeline.device`.
- **Configurable**: Single YAML config controls all parameters; environment variables override any field.
- **Modular**: Each stage can run independently or chain together via the pipeline orchestrator.

## Project Structure

```
T2Lead/
├── configs/default_config.yaml    # All pipeline parameters
├── src/drugpipe/
│   ├── config.py                  # Config loader
│   ├── pipeline.py                # Orchestrator + CLI entry
│   ├── target_discovery/          # Stage 1
│   │   ├── opentargets.py
│   │   ├── origene_client.py
│   │   └── target_ranker.py
│   ├── target_to_hit/             # Stage 2
│   │   ├── chembl_api.py
│   │   ├── dataset.py
│   │   ├── featurizer.py
│   │   ├── models.py
│   │   ├── screener.py
│   │   └── filters.py
│   ├── hit_to_lead/               # Stage 3
│   │   ├── scaffold.py
│   │   ├── clustering.py
│   │   ├── analog_gen.py
│   │   ├── mpo.py
│   │   ├── reinvent_bridge.py
│   │   └── lead_ranker.py
│   ├── lead_optimization/         # Stage 4
│   │   ├── protein_prep.py
│   │   ├── docking.py
│   │   ├── admet_deep.py
│   │   ├── md_simulation.py
│   │   └── lead_optimizer.py
│   └── utils/                     # Shared utilities
│       ├── chem.py
│       ├── io.py
│       └── http.py
├── scripts/
│   ├── run_pipeline.py            # Run full pipeline
│   └── run_stage.py               # Run a single stage
└── data/                          # Output directory (gitignored)
    ├── logs/                      # Full + summary run logs
    ├── fp_cache/                  # Cached Morgan fingerprints
    └── <disease>/                 # Per-disease outputs
```

## Installation

### Quick Start (Makefile)

From the repo root, **one command** creates the `t2lead` conda env (if missing), installs RDKit plus the OpenMM / MD stack from conda-forge, `pip install -e ".[docking,h2l]"`, PyTorch, downloads the default CReM fragment database, clones **REINVENT4** into `./REINVENT4`, fetches the prior checkpoint, runs the optional prior-metadata fix, and writes **`.env`** with `DP_*` paths for REINVENT4 and CReM. Expect a **large download** and **tens of minutes** on first run.

```bash
cd T2Lead

# Default: CUDA 12.4 PyTorch wheels (cu124) — good for RTX 4090/4080 and other Ada GPUs; use a GPU when you have one.
# Override only when needed:
#   make install TORCH_INDEX_URL=https://download.pytorch.org/whl/cu128   # Blackwell (RTX 50)
#   make install TORCH_INDEX_URL=https://download.pytorch.org/whl/cu118   # Ampere (RTX 30)
#   make install TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu     # CPU-only machines
make install

# Run full pipeline
make run DISEASE="breast cancer"
```

Partial targets (e.g. refresh REINVENT4 only): `make install-reinvent4` (alias of `install-reinvent4-full`), `make download-reinvent4-prior`, `make install-env` (re-write `DP_*` in `.env` after manual path changes).

**Makefile extras:** The Makefile prepends `$(HOME)/miniconda3/bin` to `PATH`. If Conda is installed elsewhere (Anaconda, Miniforge, `/opt/conda`, etc.), pass `CONDA_ROOT`, e.g. `make install CONDA_ROOT=/opt/conda`. The `all` optional extra in `pyproject.toml` only groups pip packages (`torch`, CReM, docking); it does **not** replace the **conda-forge** RDKit + OpenMM stack—use `make install` or the manual conda step for MD and full pipeline support.

### Step-by-step Installation

The steps below mirror **`make install`** if you cannot use Make.

### Prerequisites

- Python >= 3.9
- RDKit (install via conda or pip)
- NVIDIA GPU recommended (CUDA) for MLP training + MD simulation

> **GPU Compatibility**: RTX 50-series (Blackwell, e.g. RTX 5090) requires PyTorch with CUDA 12.8+ (`cu128` wheels). RTX 40-series (Ada) works with CUDA 12.4+ (`cu124`). RTX 30-series (Ampere) works with CUDA 11.8+ (`cu118`). See step 4 below.

### Steps

```bash
# 1. Clone / enter the project
cd T2Lead

# 2. Create a conda environment (recommended for RDKit + OpenMM)
conda create -n t2lead python=3.11 -y
conda init && source ~/.bashrc	# Run only after first Conda installation
conda activate t2lead
conda install -c conda-forge rdkit openmm pdbfixer mdtraj openmmforcefields openff-toolkit -y

# 3. Install package + docking + CReM extras (versions in pyproject.toml)
pip install -e ".[docking,h2l]"

# 4. Deep learning — choose ONE line matching your GPU (or CPU):
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu124  # RTX 4090/4080 (Ada); Makefile default
# pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128  # RTX 5090/5080 (Blackwell)
# pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118  # RTX 3090/3080 (Ampere)
# pip install torch torchvision  # CPU only (PyPI)

# 5. CReM fragment database (same default as Makefile)
bash scripts/download_crem_db.sh

# 6. REINVENT4 + prior (default directory matches Makefile: ./REINVENT4)
git clone https://github.com/MolecularAI/REINVENT4.git REINVENT4
pip install -e ./REINVENT4
mkdir -p REINVENT4/priors
curl -L "https://zenodo.org/api/records/15641297/files/reinvent.prior/content" \
  -o REINVENT4/priors/reinvent.prior
python scripts/fix_reinvent_prior_metadata.py REINVENT4/priors/reinvent.prior || true

# 7. .env — same as the end of `make install` (DP_* paths for REINVENT4 + CReM)
cp .env.example .env
python scripts/bootstrap_env.py --env-file .env --conda-prefix "$CONDA_PREFIX" \
  --prior-path "$(realpath REINVENT4/priors/reinvent.prior)" \
  --crem-db "$(realpath data/crem_db/chembl33_sa25_f5.db)"
```

### Environment Variables

**`make install`** copies `.env.example` to `.env` when missing and appends `DP_HIT_TO_LEAD__REINVENT4__REINVENT_PATH`, `DP_HIT_TO_LEAD__REINVENT4__PRIOR_PATH`, and `DP_HIT_TO_LEAD__ANALOG_GEN__CREM_DB_PATH` for the current machine (idempotent; re-run `make install-env` after moving conda or the repo).

Additional overrides (shell or `.env`):

```bash
# Recommended on AutoDL / containers with small root disks — put outputs on the large mount:
export DP_PIPELINE__OUT_DIR=/autodl-fs/data/T2Lead

# Optional: reuse a CReM DB elsewhere instead of the default under data/crem_db/
# export DP_HIT_TO_LEAD__ANALOG_GEN__CREM_DB_PATH=/path/to/chembl33_sa25_f5.db
```

## Usage

### Scenario A — Known Target in ChEMBL (default)

```bash
# Run all four stages for a disease (auto-discovers targets)
python scripts/run_pipeline.py --disease "breast cancer" -v

# Run with a specific ChEMBL target (skip Stage 1)
python scripts/run_pipeline.py --target CHEMBL4005 --stages target_to_hit hit_to_lead lead_optimization

# Custom config + verbose logging
python scripts/run_pipeline.py -c my_config.yaml -v
```

Set the protein structure for docking in `configs/default_config.yaml`:

```yaml
lead_optimization:
  pdb_id: "4JPS"    # experimental crystal structure from RCSB PDB
```

### Single Stage

```bash
# Stage 1: Target Discovery
python scripts/run_stage.py target_discovery --disease "breast cancer"

# Stage 2: Target to Hit
python scripts/run_stage.py target_to_hit --target CHEMBL4005

# Stage 3: Hit to Lead
python scripts/run_stage.py hit_to_lead

# Stage 4: Lead Optimization (reads leads from data/<disease>/final_lead_candidates.csv)
python scripts/run_stage.py lead_optimization
```

### Novel Targets (Not in ChEMBL)

For targets with little or no ChEMBL data, Stage 2 supports three alternative modes:

**Scenario B — User-supplied IC50 data** (target has assay data from literature / BindingDB / your own lab, but not in ChEMBL):

```bash
# Provide your own IC50 CSV (must have: molecule_chembl_id, target_chembl_id,
# standard_value, standard_type, standard_units)
python scripts/run_pipeline.py \
  --target MY_TARGET_ID \
  --activities-csv /path/to/my_ic50_data.csv \
  -v

# Optionally also provide a custom screening library
python scripts/run_pipeline.py \
  --target MY_TARGET_ID \
  --activities-csv /path/to/my_ic50_data.csv \
  --screening-library /path/to/my_compounds.csv \
  -v
```

**Scenario C — Docking-only mode** (truly novel target with zero IC50 activity data):

`--docking-only` skips the ML pipeline (RandomForest + MLP training and virtual screening) because there is no IC50 data to train on. Stage 1 is not needed if you already know your target (use `--target` directly). Stage 2 filters candidates by drug-likeness only. Stage 3 still generates analogs via CReM. Stage 4 docking becomes the primary scoring method.

```bash
# With a known PDB structure:
python scripts/run_pipeline.py \
  --target MY_NOVEL_TARGET \
  --docking-only \
  -v

# With only an amino acid sequence (no PDB needed — ESMFold predicts the 3D structure automatically):
python scripts/run_pipeline.py \
  --target MY_NOVEL_TARGET \
  --docking-only \
  --protein-sequence "MTEYKLVVVGAVGVGKSALT..." \
  -v
```

The pipeline automatically resolves the protein structure in this priority:

| Priority | Source | What you provide |
|----------|--------|------------------|
| 1 | RCSB PDB | `pdb_id: "4JPS"` in config |
| 2 | Local PDB file | Place `.pdb` file in output directory |
| 3 | **ESMFold API** (automatic) | `--protein-sequence "MTEYKLVV..."` or `protein_sequence` in config |

ESMFold (Meta AI) predicts AlphaFold-quality 3D structures from amino acid sequences via a free REST API — no signup, no GPU, fully integrated into the pipeline.

> **Sequence length limit**: ESMFold API works best for sequences under 400 residues. For longer proteins, download from [AlphaFold DB](https://alphafold.ebi.ac.uk/) or use [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb).

### As a Python Library

```python
from drugpipe.config import load_config
from drugpipe.pipeline import (
    run_target_discovery, run_target_to_hit,
    run_hit_to_lead, run_lead_optimization,
)

cfg = load_config(overrides={
    "target_discovery": {"disease": "breast cancer"},
    "pipeline": {"stages": ["target_discovery", "target_to_hit", "hit_to_lead", "lead_optimization"]},
    "lead_optimization": {"pdb_id": "4JPS"},
})

# Stage 1
targets = run_target_discovery(cfg)
target_id = targets[0]["chembl_id"]

# Stage 2
df_hits = run_target_to_hit(cfg, target_chembl_id=target_id)

# Stage 3
trainer = df_hits.attrs.get("_trainer")
featurizer = df_hits.attrs.get("_featurizer")
df_leads = run_hit_to_lead(
    cfg, df_hits,
    model_predict_fn=trainer.predict if trainer else None,
    featurizer_fn=featurizer.transform if featurizer else None,
)

# Stage 4
df_optimized = run_lead_optimization(cfg, df_leads)
```

## Configuration

All parameters live in `configs/default_config.yaml`. Key sections:

| Section | What it controls |
|---|---|
| `pipeline.stages` | Which stages to run (`target_discovery`, `target_to_hit`, `hit_to_lead`, `lead_optimization`) |
| `pipeline.seed` | Global random seed |
| `pipeline.device` | Compute device: `auto`, `cuda`, `mps`, `cpu` |
| `target_discovery.disease` | Input disease name or EFO ID |
| `target_discovery.origene.enabled` | Whether to call a local OriGene service |
| `target_to_hit.chembl.*` | ChEMBL crawl parameters (page_limit, max records, etc.) |
| `target_to_hit.model.*` | ML model hyperparameters (RF, MLP) |
| `target_to_hit.filter.*` | ADMET rule thresholds, pIC50 cutoff, QED minimum |
| `hit_to_lead.analog_gen.*` | CReM analog generation settings |
| `hit_to_lead.mpo.*` | MPO weight distribution |
| `hit_to_lead.reinvent4.*` | REINVENT4 RL optimization (paths, steps, batch size) |
| `lead_optimization.pdb_id` | RCSB PDB ID for target protein structure |
| `lead_optimization.docking.*` | AutoDock Vina docking parameters |
| `lead_optimization.admet_deep.*` | Enhanced ADMET profiling settings |
| `lead_optimization.md_simulation.*` | OpenMM MD simulation settings |
| `lead_optimization.scoring.*` | Composite scoring weights |

Override any value via environment variables with the `DP_` prefix:

```bash
export DP_TARGET_TO_HIT__CHEMBL__MAX_ACTIVITIES=50000
export DP_LEAD_OPTIMIZATION__PDB_ID=4JPS
```

## Stage Details

### Stage 1: Target Discovery

**Input**: disease name (e.g. `"lung cancer"`) or EFO ID (e.g. `EFO_0001378`) — a plain string.

**What happens**:
- Queries OpenTargets Platform GraphQL API for disease-target associations
- Maps Ensembl gene IDs to ChEMBL target IDs
- Optionally merges OriGene AI recommendations (requires separate [OriGene](https://github.com/GENTEL-lab/OriGene) deployment)

**Output**: ranked list of ChEMBL target IDs with association scores (printed + passed to Stage 2).

### Stage 2: Target to Hit

**Input**: a ChEMBL target ID string (from Stage 1 or user-specified, e.g. `"CHEMBL4005"`).

**What happens**:
1. **ChEMBL Crawl** — download molecules (SMILES) and IC50 activity data via REST API, with checkpoint resume
2. **Target Selection** — automatically picks the highest-scored target with enough IC50 data (>= 200 records)
3. **Dataset Build** — join molecules + activities, convert IC50 to pIC50, aggregate per-molecule duplicates
4. **Featurize** — compute 2048-bit Morgan fingerprints (ECFP4)
5. **Train** — Random Forest regressor + Torch MLP (GPU-accelerated); ensemble prediction
6. **Virtual Screen** — score all crawled molecules with the trained model
7. **ADMET Filter** — Lipinski rules, TPSA, QED >= 0.6, PAINS/Brenk structural alerts, pIC50 >= 6.0

**Output**: `data/<disease>/final_hit_candidates.csv`

### Stage 3: Hit to Lead

**Input**: hit compounds CSV from Stage 2 (or any CSV with `canonical_smiles` and `pred_pIC50_ens` columns).

**What happens**:
1. **Scaffold Analysis** — extract Murcko scaffolds, group by generic scaffold, SAR summary
2. **Diversity Clustering** — Butina clustering on Tanimoto distance
3. **Analog Generation** — CReM fragment mutation to expand chemical space (requires CReM + fragment DB)
4. **MPO Scoring** — weighted: 40% potency + 25% QED + 20% ADMET + 15% novelty
5. **REINVENT4 RL** — optional reinforcement-learning molecular generation with potency predictor as reward
6. **Lead Ranking** — final sorted output, top-N

**Output**: `data/<disease>/final_lead_candidates.csv`

### Stage 4: Lead Optimization

**Input**: lead candidates CSV from Stage 3 + target protein PDB ID.

**What happens**:
1. **Protein Preparation** — fetch PDB from RCSB, fix with PDBFixer (add missing atoms, hydrogens), auto-detect binding site from co-crystallized ligand, convert to PDBQT
2. **Molecular Docking** — dock all leads with AutoDock Vina (CPU multi-threaded), score binding affinity (kcal/mol)
3. **Enhanced ADMET** — Synthetic Accessibility score, hERG cardiotoxicity SMARTS, CYP3A4/2D6 inhibition risk, Veber oral bioavailability rules, composite risk score
4. **MD Simulation** — (top-N only) OpenMM implicit-solvent energy minimization + MM-GBSA binding energy (GPU-accelerated via CUDA)
5. **Composite Scoring** — weighted combination of docking score, ADMET risk, MD binding energy → final ranking

**Output**: `data/<disease>/optimized_leads.csv`

## Output Files

All outputs go to `data/<disease>/` (configurable via `pipeline.out_dir`):

| File | Description | Produced by |
|---|---|---|
| `fp_cache/morgan_*.npy` | Cached Morgan fingerprints (auto-generated) | Stage 2 |
| `molecules_chemblid_smiles.csv` | Crawled molecules | Stage 2 |
| `activities_ic50.csv` | Crawled IC50 activity data | Stage 2 |
| `dataset_target_ic50.csv` | Training dataset for the selected target | Stage 2 |
| `scored_candidates.csv` | All molecules with predicted pIC50 | Stage 2 |
| `final_hit_candidates.csv` | Filtered hit compounds | Stage 2 |
| `h2l_scaffold_summary.csv` | Scaffold SAR summary | Stage 3 |
| `h2l_cluster_summary.csv` | Diversity cluster summary | Stage 3 |
| `final_lead_candidates.csv` | Final ranked lead candidates | Stage 3 |
| `<PDB_ID>.pdb` | Downloaded protein structure | Stage 4 |
| `<PDB_ID>_receptor.pdbqt` | Prepared receptor for docking | Stage 4 |
| `docking_poses/` | Best docking poses per molecule | Stage 4 |
| `optimized_leads.csv` | Final optimized candidates with all scores | Stage 4 |

## External Dependencies

| Component | Required? | Install | Used in |
|---|---|---|---|
| RDKit | Yes | `conda install -c conda-forge rdkit` | All stages |
| PyTorch | Optional | `pip install torch` | Stage 2 (MLP), Stage 3 (REINVENT4) |
| CReM | Optional | `pip install crem` + [fragment DB](https://github.com/DrrDom/crem#databases) | Stage 3 (analog gen) |
| REINVENT4 | Optional | `git clone` + `pip install -e .` from [REINVENT4](https://github.com/MolecularAI/REINVENT4) | Stage 3 (RL optimization) |
| AutoDock Vina | Optional | `pip install vina meeko gemmi` | Stage 4 (docking) |
| OpenMM | Optional | `conda install -c conda-forge openmm pdbfixer mdtraj` | Stage 4 (MD simulation) |
| openmmforcefields | Optional | `conda install -c conda-forge openmmforcefields openff-toolkit` | Stage 4 (MD ligand parameterization) |
| OriGene | Optional | See [OriGene](https://github.com/GENTEL-lab/OriGene) | Stage 1 (AI target rec.) |

All optional dependencies degrade gracefully — if not installed, the corresponding step is skipped with a warning.

## Somatic Variant Pipeline (nf-core/sarek)

T2Lead can consume tumor/normal somatic variant calls and route mutant protein structures directly into Stage 4 docking and MD. The upstream variant calling is handled by **nf-core/sarek 3.8.1** (Nextflow + Docker), producing VEP-annotated VCFs that the `variant_analysis` module parses.

### Overview

```
Tumor FASTQ ──┐
              ├─→ [nf-core/sarek] ──→ VEP-annotated VCF
Normal FASTQ ─┘  (BWA-MEM2 → GATK Mutect2 → VEP 115)
                                             │
                        T2Lead variant_analysis module
                        (vcf_parser → mutant_sequence → structure_bridge)
                                             │
                                    Stage 4: Docking + MD
                                 (mutant structure from ESMFold / FoldX)
```

### Prerequisites

- **Docker** (or Singularity/Apptainer for HPC) — nf-core/sarek pulls containers automatically
- **Java 21+** — required by Nextflow (`java -version`)
- **16 cores / 32 GB RAM** minimum for WES; 32 cores / 64 GB for WGS 30×

### Install Nextflow + Sarek

```bash
# Install Nextflow to ~/.local/bin (requires Java 21+)
make install-nextflow

# Pull nf-core/sarek 3.8.1 containers (Docker must be running)
make install-sarek

# Or manually:
curl -fsSL https://get.nextflow.io | bash
mv nextflow ~/.local/bin/
nextflow pull nf-core/sarek -r 3.8.1
```

### Run Somatic Variant Calling

Prepare a **samplesheet.csv** (tumor/normal paired):

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
patient1,XX,1,tumor_S1,lane_1,/data/tumor_R1.fastq.gz,/data/tumor_R2.fastq.gz
patient1,XX,0,normal_S1,lane_1,/data/normal_R1.fastq.gz,/data/normal_R2.fastq.gz
```

Then run:

```bash
# Via Makefile (Docker profile, GRCh38, Mutect2 + VEP annotation)
make run-sarek INPUT=samplesheet.csv

# Or directly:
nextflow run nf-core/sarek \
  -r 3.8.1 \
  -profile docker \
  --input samplesheet.csv \
  --genome GRCh38 \
  --tools mutect2,vep \
  --outdir results/sarek
```

Override defaults: `make run-sarek INPUT=samplesheet.csv SAREK_PROFILE=singularity SAREK_GENOME=GRCh37`

### Handoff to T2Lead

After Sarek completes, pass the annotated VCF to T2Lead's `variant_analysis` module (Stage 0):

```python
from drugpipe.variant_analysis.vcf_parser import VCFParser
from drugpipe.variant_analysis.mutant_sequence import MutantSequenceBuilder

# Parse VEP-annotated VCF → driver missense mutations
variants = VCFParser(driver_genes=["EGFR", "PIK3CA", "BRAF"]).parse(
    "results/sarek/annotation/sample/sample.ann.vcf"
)

# Build mutant protein sequences → FASTA for ESMFold
proteins = MutantSequenceBuilder().build(variants)
```

The mutant protein structure (via ESMFold or FoldX) then feeds directly into Stage 4 docking and MD. See [`docs/research/somatic_variant_pipeline_feasibility.md`](docs/research/somatic_variant_pipeline_feasibility.md) ([中文版](docs/research/somatic_variant_pipeline_feasibility_zh.md)) for a full feasibility analysis and tool landscape.

### Runtime Estimates (WES ~100×, 16 cores)

| Step | Tool | Wall Time |
|------|------|-----------|
| Alignment | BWA-MEM2 | ~30–45 min per sample |
| Sort + Dedup | SAMtools | ~15–20 min per sample |
| Somatic calling | GATK Mutect2 | ~2–4 hours |
| Annotation | VEP 115 | ~10–20 min |
| **Total upstream** | | **~4–6 hours** |
| ESMFold (per protein) | ESM-2 | ~30 s–2 min (GPU) |
| Docking + MD | Vina + OpenMM | as per Stage 4 |

---

## Testing

Run the test suite (no GPU, no conda required — works with a minimal pip install):

```bash
# Via Makefile
make test

# Or directly
python -m pytest tests/ -v
```

### Test Files

| File | What it tests | Key dependencies |
|------|---------------|------------------|
| `tests/test_explicit_md.py` | `ExplicitSolventRefiner.should_trigger()` — decision logic for activating TIP3P explicit-solvent MD refinement based on opt\_score spread among top-N candidates | `pandas`, `numpy` (no RDKit, no OpenMM) |
| `tests/test_md_energy_compat.py` | `_energy_as_kcal()` — converts OpenMM `Quantity` objects **or** plain floats to kcal/mol; regression test for the `float has no attribute value_in_unit` bug | `pytest`; OpenMM tests auto-skipped when not installed |
| `tests/test_variant_analysis.py` | VCF parser (VEP-annotated VCF → `SomaticVariant`), amino acid 3-letter→1-letter conversion, mutant sequence building with mocked UniProt fetch, FASTA output | `pytest` only; no chemistry libraries needed |

All tests pass without the conda/RDKit/OpenMM stack:

```
25 passed, 3 skipped   ← 3 OpenMM-only tests auto-skipped when OpenMM is absent
```

---

## Acknowledgements & References

| Stage | Reference project | What we learned / borrowed |
|---|---|---|
| **Stage 1** | [OriGene](https://github.com/GENTEL-lab/OriGene) (GENTEL-lab, CC-BY-NC-SA 4.0) | Multi-agent AI approach to therapeutic target discovery; our `origene_client.py` provides an integration wrapper for their service. |
| **Stage 1** | [OpenTargets Platform](https://platform.opentargets.org/) | Public GraphQL API for disease-target associations, genetic evidence, and tractability data. |
| **Stage 2** | [ChEMBL REST API](https://www.ebi.ac.uk/chembl/) | Public compound and bioactivity database maintained by EMBL-EBI. |
| **Stage 3** | [CReM](https://github.com/DrrDom/crem) (Polishchuk, 2020) | Context-aware fragment replacement for analog enumeration. |
| **Stage 3** | [REINVENT4](https://github.com/MolecularAI/REINVENT4) (MolecularAI, Apache 2.0) | RL-based generative molecular design; our bridge module generates REINVENT4 configs and parses its output. |
| **Stage 4** | [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) (Scripps Research) | Molecular docking engine with Python bindings. |
| **Stage 4** | [OpenMM](https://openmm.org/) (Stanford, MIT License) | GPU-accelerated molecular dynamics toolkit for MM-GBSA binding energy calculations. |
| **All** | [RDKit](https://www.rdkit.org/) | Murcko scaffolds, Butina clustering, PAINS/Brenk filters, QED, SA score, molecular descriptors. |

## Contributors

- **Lu Zhentao** — Lead Developer, responsible for overall architecture and core implementation (during internship at hupper)
- **Jiang Keying** — Developer, responsible for overall architecture and core implementation (during internship at hupper)

## License

Apache-2.0 © 2026 hupper
