# Disease-driven workflow (Stages 1–4)

[中文版 (disease_pipeline_zh.md)](disease_pipeline_zh.md)

This path demonstrates end-to-end discovery from a **disease label** through public data and generative tools, then **Stage 4** rescoring.

## Run

```bash
python scripts/run_pipeline.py --disease "breast cancer" -v
```

Skip Stage 1 with a known ChEMBL target:

```bash
python scripts/run_pipeline.py --target CHEMBL4005 \
  --stages target_to_hit hit_to_lead lead_optimization -v
```

## Custom data (novel target)

IC50 CSV columns: `molecule_chembl_id`, `target_chembl_id`, `standard_value`, `standard_type`, `standard_units`.

```bash
python scripts/run_pipeline.py \
  --target MY_TARGET_ID \
  --activities-csv data/user_inputs/activities/my_ic50.csv \
  --screening-library data/user_inputs/screening_library/my_mols.csv \
  -v
```

Docking-only (no IC50 to train on):

```bash
python scripts/run_pipeline.py --target MY_TARGET --docking-only -v
```

## Stage outputs

With `pipeline.output_layout.use_stage_subdirs: true`, see [output_reference.md](output_reference.md).

## Single-stage CLI

```bash
python scripts/run_stage.py target_discovery --disease "lung cancer"
python scripts/run_stage.py target_to_hit --target CHEMBL204
python scripts/run_stage.py hit_to_lead   # reads stage2_hits/final_hit_candidates.csv by default
```

`hit_to_lead` resolves the hit CSV from `run_root_for_config` + `stage2_hits/` (same logic as the full pipeline).
