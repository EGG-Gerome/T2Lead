# Output files

[中文版 (output_reference_zh.md)](output_reference_zh.md)

Assume run root: `<pipeline.out_dir>/<disease_slug>/` when a disease is set, else `<pipeline.out_dir>/`. With `use_stage_subdirs: false`, all stage keys collapse to the same directory.

## Shared ChEMBL library (optional)

When `target_to_hit.shared_library_dir` is set (default in `default_config.yaml`), these files are stored under `<pipeline.out_dir>/<shared_library_dir>/` (not under each disease folder) so multiple diseases reuse one crawl and one Morgan `fp_cache/`:

| File / dir | Description |
|------------|-------------|
| `molecules_chemblid_smiles.csv` | Crawled molecules |
| `activities_ic50.csv` | Crawled IC50 rows |
| `crawl_state.json` | Resume checkpoint |
| `fp_cache/` | Morgan fingerprint caches |

Per-disease Stage 2 still holds `dataset_target_ic50.csv`, `model_cache/`, `scored_candidates.csv`, `final_hit_candidates.csv`, etc.

If `shared_library_dir` is empty, legacy layout keeps all of the above inside each run’s `stage2_hits/`.

## Stage 1 — `stage1_targets/`

| File | Description |
|------|-------------|
| `ranked_targets.csv` | Ranked targets from OpenTargets |

## Stage 2 — `stage2_hits/`

| File | Description |
|------|-------------|
| `molecules_chemblid_smiles.csv` | Crawled molecules (ChEMBL mode); see **Shared ChEMBL library** if `shared_library_dir` is set |
| `activities_ic50.csv` | Crawled IC50 rows (same note) |
| `crawl_state.json` | Resume checkpoint (same note) |
| `dataset_target_ic50.csv` | Training table |
| `model_cache/` | Serialized RF / MLP |
| `fp_cache/` | Morgan caches (same note; otherwise under `stage2_hits/`) |
| `scored_candidates.csv` | VS output |
| `final_hit_candidates.csv` | After ADMET / potency filters |

## Stage 3 — `stage3_leads/`

| File | Description |
|------|-------------|
| `h2l_scaffold_summary.csv` | Scaffold SAR |
| `h2l_cluster_summary.csv` | Clustering |
| `final_lead_candidates.csv` | MPO-ranked leads |
| `reinvent4_*` | Optional RL artifacts |

## Stage 4 — `stage4_optimization/`

| File | Description |
|------|-------------|
| `<PDB>.pdb`, `*_fixed.pdb`, `*_receptor.pdbqt` | Receptor |
| `docking_poses/pose_*.pdbqt` | Poses |
| `md_trajectories/` | MD logs (if enabled) |
| `optimized_leads.csv` | Final table with `docking_score`, `md_binding_energy`, `md_binding_energy_std`, `md_rmsd_mean`, `opt_score`, … |

Variant-only runs may use nested folders: `stage4_optimization/EGFR_L858R/optimized_leads.csv`.

## Other

| Path | Description |
|------|-------------|
| `logs/` | Timestamped `*_full.log`, `*_summary.log` |
| `variant_analysis/` | `mutant_fastas/`, `structures/` (under run root) |
| `variant_calling/` | Sarek samplesheet + outputs when FASTQs drive calling |
