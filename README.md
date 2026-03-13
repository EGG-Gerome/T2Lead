# T2Lead

[中文说明 (README_zh.md)](README_zh.md)

**Target → Hit → Lead**: an end-to-end, modular drug discovery pipeline that takes a disease name as input and outputs ranked lead candidate molecules.

```
Disease name ──► Target Discovery ──► Target-to-Hit ──► Hit-to-Lead ──► Lead Candidates CSV
 (string)        (OpenTargets /       (ChEMBL crawl /   (Scaffold /
                  OriGene)             ML + VS + ADMET)   CReM / REINVENT4)
```

## Quick I/O Summary

| Stage | Input | Output |
|---|---|---|
| **1 — Target Discovery** | Disease name `str` (e.g. `"lung cancer"`) or EFO ID (e.g. `EFO_0001378`) | Ranked target list: `chembl_id, symbol, name, score, source` |
| **2 — Target to Hit** | ChEMBL target ID `str` (e.g. `"CHEMBL204"`) | `data/final_hit_candidates.csv` — columns: `molecule_chembl_id, canonical_smiles, pred_pIC50_ens, pred_IC50_nM_ens, QED, MW, cLogP, TPSA, HBD, HBA, RotB, Rings, HeavyAtoms` |
| **3 — Hit to Lead** | Hit CSV from Stage 2 (must contain `canonical_smiles`, `pred_pIC50_ens`) | `data/final_lead_candidates.csv` — columns: `canonical_smiles, source_smiles, origin, pred_pIC50_ens, QED, mpo_score, admet_pass, HasAlert, scaffold_smi, cluster_id` |

## Important: What Happens Under the Hood

This is **NOT** a pre-trained model you just run. Here's what each stage actually does:

- **Stage 1** queries the OpenTargets public API online (seconds).
- **Stage 2** crawls ChEMBL to download molecule + IC50 data first (can take hours depending on `MAX_MOLECULES` / `MAX_ACTIVITIES` settings), then trains Random Forest + MLP models **on the fly** using the crawled data for your specific target. No pre-built model — it learns from real bioactivity data each run.
- **Stage 3** analyzes scaffolds and clusters (instant), generates analogs via CReM if a fragment DB is provided, scores everything with MPO.

## Features

- **Stage 1 — Target Discovery**: Query OpenTargets GraphQL API by disease name or EFO ID; optionally merge results with an OriGene AI agent for LLM-assisted target recommendation.
- **Stage 2 — Target to Hit**: Crawl ChEMBL molecules + IC50 activities (with checkpoint resume), build pIC50 datasets, train Random Forest + Torch MLP regressors, run virtual screening on the full molecule library, and filter by ADMET rules / QED / PAINS alerts.
- **Stage 3 — Hit to Lead**: Murcko scaffold SAR analysis, Butina diversity clustering, CReM-based analog enumeration, multi-parameter optimization (MPO) scoring, and optional REINVENT4 RL-based molecular generation.
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
│   └── utils/                     # Shared utilities
│       ├── chem.py
│       ├── io.py
│       └── http.py
├── scripts/
│   ├── run_pipeline.py            # Run full pipeline
│   └── run_stage.py               # Run a single stage
├── data/                          # Output directory (gitignored)
└── notebooks/
    └── analysis.ipynb
```

## Installation

### Prerequisites

- Python >= 3.9
- RDKit (install via conda or pip)

### Steps

```bash
# 1. Clone / enter the project
cd T2Lead

# 2. Create a conda environment (recommended for RDKit)
conda create -n t2lead python=3.11 -y
conda init && source ~/.bashrc	# # Run only after first Conda installation
conda activate t2lead
conda install -c conda-forge rdkit -y

# 3. Install core dependencies
pip install -r requirements.txt

# 4. Install the package in editable mode
pip install -e .

# 5. (Optional) Deep learning support
pip install torch

# 6. (Optional) CReM for analog generation in Stage 3
pip install crem
# You also need a CReM fragment database (~2 GB) — see:
# https://github.com/DrrDom/crem#databases
```

### Environment Variables

Copy `.env.example` to `.env` and fill in any needed keys:

```bash
cp .env.example .env
```

## Usage

### Full Pipeline

```bash
# Run all three stages for a disease
python scripts/run_pipeline.py --disease "lung cancer"

# Run with a specific target (skip Stage 1)
python scripts/run_pipeline.py --target CHEMBL204 --stages target_to_hit hit_to_lead

# Custom config + verbose logging
python scripts/run_pipeline.py -c my_config.yaml -v
```

### Single Stage

```bash
# Stage 1: Target Discovery
python scripts/run_stage.py target_discovery --disease "breast cancer"

# Stage 2: Target to Hit
python scripts/run_stage.py target_to_hit --target CHEMBL1862

# Stage 3: Hit to Lead (reads hits from data/final_hit_candidates.csv)
python scripts/run_stage.py hit_to_lead

# Stage 3 with a custom hits file
python scripts/run_stage.py hit_to_lead --hits-csv data/final_hit_candidates.csv
```

### As a Python Library

```python
from drugpipe.config import load_config
from drugpipe.pipeline import run_target_discovery, run_target_to_hit, run_hit_to_lead

cfg = load_config(overrides={
    "target_discovery": {"disease": "liver cancer"},
    "pipeline": {"stages": ["target_discovery", "target_to_hit", "hit_to_lead"]},
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
```

## Configuration

All parameters live in `configs/default_config.yaml`. Key sections:

| Section | What it controls |
|---|---|
| `pipeline.stages` | Which stages to run (`target_discovery`, `target_to_hit`, `hit_to_lead`) |
| `pipeline.seed` | Global random seed |
| `target_discovery.disease` | Input disease name or EFO ID |
| `target_discovery.origene.enabled` | Whether to call a local OriGene service |
| `target_to_hit.chembl.*` | ChEMBL crawl parameters (page_limit, max records, etc.) |
| `target_to_hit.model.*` | ML model hyperparameters (RF, MLP) |
| `target_to_hit.filter.*` | ADMET rule thresholds, pIC50 cutoff, QED minimum |
| `hit_to_lead.analog_gen.*` | CReM analog generation settings |
| `hit_to_lead.mpo.*` | MPO weight distribution |
| `hit_to_lead.reinvent4.*` | Optional REINVENT4 RL optimization |

Override any value via environment variables with the `DP_` prefix:

```bash
export DP_TARGET_TO_HIT__CHEMBL__MAX_ACTIVITIES=50000
export DP_HIT_TO_LEAD__MPO__W_POTENCY=0.5
```

## Stage Details

### Stage 1: Target Discovery

**Input**: disease name (e.g. `"lung cancer"`) or EFO ID (e.g. `EFO_0001378`) — a plain string.

**What happens**:
- Queries OpenTargets Platform GraphQL API for disease-target associations
- Maps Ensembl gene IDs to ChEMBL target IDs
- Optionally merges OriGene AI recommendations (requires separate OriGene deployment)

**Output**: ranked list of ChEMBL target IDs with association scores (printed + passed to Stage 2).

**Prerequisites**: internet connection only. No local data or models required.

### Stage 2: Target to Hit

**Input**: a ChEMBL target ID string (from Stage 1 or user-specified, e.g. `"CHEMBL204"`).

**What happens**:
1. **ChEMBL Crawl** — download molecules (SMILES) and IC50 activity data via REST API, with checkpoint resume (hours depending on scale)
2. **Dataset Build** — join molecules + activities, convert IC50 to pIC50, aggregate per-molecule duplicates
3. **Featurize** — compute 2048-bit Morgan fingerprints (ECFP4)
4. **Train** — Random Forest regressor + optional Torch MLP, trained on crawled data; ensemble prediction
5. **Virtual Screen** — score all crawled molecules with the trained model
6. **ADMET Filter** — Lipinski rules, TPSA, QED >= 0.6, PAINS/Brenk structural alerts, pIC50 >= 6.0

**Output**: `data/final_hit_candidates.csv`

**Prerequisites**: internet connection for ChEMBL API. Models are trained on the fly — no pre-trained weights needed.

### Stage 3: Hit to Lead

**Input**: hit compounds CSV from Stage 2 (or any CSV with `canonical_smiles` and `pred_pIC50_ens` columns).

**What happens**:
1. **Scaffold Analysis** — extract Murcko scaffolds, group by generic scaffold, SAR summary
2. **Diversity Clustering** — Butina clustering on Tanimoto distance
3. **Analog Generation** — CReM fragment mutation to expand chemical space (requires CReM + fragment DB)
4. **MPO Scoring** — weighted: 40% potency + 25% QED + 20% ADMET + 15% novelty
5. **Lead Ranking** — final sorted output, top-N

**Output**: `data/final_lead_candidates.csv`

**Prerequisites**: RDKit only for scaffold/clustering/MPO. CReM + fragment DB for analog generation (optional). REINVENT4 for RL optimization (optional).

## Output Files

All outputs go to `data/` (configurable via `pipeline.out_dir`):

| File | Description | Produced by |
|---|---|---|
| `molecules_chemblid_smiles.csv` | Crawled molecules | Stage 2 |
| `activities_ic50.csv` | Crawled IC50 activity data | Stage 2 |
| `dataset_target_ic50.csv` | Training dataset for the selected target | Stage 2 |
| `scored_candidates.csv` | All molecules with predicted pIC50 | Stage 2 |
| `final_hit_candidates.csv` | Filtered hit compounds | Stage 2 |
| `h2l_scaffold_summary.csv` | Scaffold SAR summary | Stage 3 |
| `h2l_cluster_summary.csv` | Diversity cluster summary | Stage 3 |
| `final_lead_candidates.csv` | Final ranked lead candidates | Stage 3 |

## External Dependencies

| Component | Required? | Install | Used in |
|---|---|---|---|
| RDKit | Yes | `conda install -c conda-forge rdkit` | All stages |
| PyTorch | Optional | `pip install torch` | Stage 2 (MLP model) |
| CReM | Optional | `pip install crem` + [fragment DB](https://github.com/DrrDom/crem#databases) | Stage 3 (analog gen) |
| REINVENT4 | Optional | See [REINVENT4](https://github.com/MolecularAI/REINVENT4) | Stage 3 (RL optimization) |
| OriGene | Optional | See [OriGene](https://github.com/GENTEL-lab/OriGene) | Stage 1 (AI target rec.) |

## Acknowledgements & References

This project was built by studying and integrating ideas from the following projects. We gratefully acknowledge them:

| Stage | Reference project | What we learned / borrowed |
|---|---|---|
| **Stage 1** | [OriGene](https://github.com/GENTEL-lab/OriGene) (GENTEL-lab, CC-BY-NC-SA 4.0) | Multi-agent AI approach to therapeutic target discovery; our `origene_client.py` provides an integration wrapper for their service. The idea of combining data-driven (OpenTargets) + AI-driven (LLM) target recommendations comes from their architecture. |
| **Stage 1** | [OpenTargets Platform](https://platform.opentargets.org/) | Public GraphQL API for disease-target associations, genetic evidence, and tractability data. |
| **Stage 2** | previous ChEMBL-based workflows and tutorials | The original single-file ChEMBL workflow (API crawl → IC50 dataset → Morgan FP → RF/MLP → virtual screening → ADMET/QED filter) was refactored into 6 modular classes. Core algorithm logic (pIC50 conversion, fingerprint featurization, RF+MLP ensemble, Lipinski/PAINS filtering) originates from this script. |
| **Stage 2** | [ChEMBL REST API](https://www.ebi.ac.uk/chembl/) | Public compound and bioactivity database maintained by EMBL-EBI. |
| **Stage 3** | [CReM](https://github.com/DrrDom/crem) (Polishchuk, 2020) | Context-aware fragment replacement for analog enumeration. |
| **Stage 3** | [REINVENT4](https://github.com/MolecularAI/REINVENT4) (MolecularAI, Apache 2.0) | RL-based generative molecular design; our bridge module generates REINVENT4 configs and parses its output. |
| **Stage 3** | [RDKit](https://www.rdkit.org/) | Murcko scaffold decomposition, Butina clustering, PAINS/Brenk filter catalogs, QED, molecular descriptors. |

## Contributors

- **Lu Zhentao** — Lead Developer, responsible for overall architecture and core implementation (during internship at hupper)
- **Jiang Keying** — Developer, responsible for overall architecture and core implementation (during internship at hupper)

## License

Apache-2.0 © 2026 hupper
