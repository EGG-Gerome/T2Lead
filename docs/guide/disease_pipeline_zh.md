# 疾病驱动工作流（阶段一至四）

[English (disease_pipeline.md)](disease_pipeline.md)

该路径演示从**疾病标签**经公开数据与生成式工具到**阶段四**再评分的端到端发现流程。

## 运行

```bash
python scripts/run_pipeline.py --disease "breast cancer" -v
```

若已知 ChEMBL 靶点，可跳过阶段一：

```bash
python scripts/run_pipeline.py --target CHEMBL4005 \
  --stages target_to_hit hit_to_lead lead_optimization -v
```

## 自定义数据（新靶点）

IC50 CSV 列：`molecule_chembl_id`、`target_chembl_id`、`standard_value`、`standard_type`、`standard_units`。

```bash
python scripts/run_pipeline.py \
  --target MY_TARGET_ID \
  --activities-csv data/user_inputs/activities/my_ic50.csv \
  --screening-library data/user_inputs/screening_library/my_mols.csv \
  -v
```

纯对接（无 IC50 训练）：

```bash
python scripts/run_pipeline.py --target MY_TARGET --docking-only -v
```

## 各阶段输出

当 `pipeline.output_layout.use_stage_subdirs: true` 时，详见 [output_reference_zh.md](output_reference_zh.md)。

## 单阶段 CLI

```bash
python scripts/run_stage.py target_discovery --disease "lung cancer"
python scripts/run_stage.py target_to_hit --target CHEMBL204
python scripts/run_stage.py hit_to_lead   # 默认读取 stage2_hits/final_hit_candidates.csv
```

`hit_to_lead` 从 `run_root_for_config` + `stage2_hits/` 解析 Hit CSV，逻辑与完整流水线一致。
