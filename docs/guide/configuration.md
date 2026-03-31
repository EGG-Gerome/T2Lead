# Configuration reference

Primary file: **`configs/default_config.yaml`**. Overrides merge in this order (lowest first):

1. Default YAML  
2. User `-c my.yaml`  
3. Environment variables `DP_*`  
4. CLI / Python `overrides`

## `DP_` mapping

Double underscores nest: `DP_TARGET_TO_HIT__CHEMBL__PAGE_LIMIT=500` → `target_to_hit.chembl.page_limit`.

## Notable keys

| Section | Keys |
|---------|------|
| `pipeline` | `stages`, `out_dir`, `device`, `output_layout.use_stage_subdirs` |
| `variant_analysis` | `enabled`, `vcf_path`, FASTQ paths, `driver_genes`, `min_impact`, `sarek.*` |
| `target_to_hit` | `external_activities_csv`, `screening_library_csv`, `docking_only`, `chembl.*`, `model.*`, `filter.*` |
| `hit_to_lead` | `analog_gen.*`, `mpo.*`, `reinvent4.*` |
| `lead_optimization` | `pdb_id`, `protein_sequence`, `docking.*`, `admet_deep.*` (incl. `hard_filter.*`), `md_simulation.*` (incl. `ensemble.*`), `scoring.*`, `explicit_md.*` |

## Stage 4 scoring defaults

- `scoring.w_docking`: 0.25  
- `scoring.w_md_energy`: 0.40  
- `scoring.w_stability`: 0.35  
- `scoring.w_cyp_soft`: 0.08 (set `0` to disable)

ADMET safety gates: `lead_optimization.admet_deep.hard_filter` (`drop_herg`, `drop_veber_fail`, `drop_high_sa`).

## Ensemble MM-GBSA

Under `lead_optimization.md_simulation.ensemble`:

- `enabled`, `n_runs`, `equilibration_ps`, `production_ps`, `sample_interval_ps`

Binding statistics are **approximate** (fixed receptor/ligand single-structure energies subtracted from time-varying complex energy).

## Autodl / disk

If `./data` is used and `/autodl-fs/data` exists, `get_out_dir()` may redirect defaults to `/autodl-fs/data/T2Lead`. Set `DP_PIPELINE__OUT_DIR` explicitly for full control.
