# 配置参考

[English (configuration.md)](configuration.md)

主配置文件：**`configs/default_config.yaml`**。覆盖项按以下顺序合并（后者优先）：

1. 默认 YAML  
2. 用户 `-c my.yaml`  
3. 环境变量 `DP_*`  
4. CLI / Python `overrides`

## `DP_` 映射规则

双下划线表示嵌套：`DP_TARGET_TO_HIT__CHEMBL__PAGE_LIMIT=500` → `target_to_hit.chembl.page_limit`。

## 常用配置项

| 区块 | 键 |
|------|-----|
| `pipeline` | `stages`、`out_dir`、`device`、`output_layout.use_stage_subdirs`、`output_layout.variant_isolated_runs`、`output_layout.variant_run_id` |
| `variant_analysis` | `enabled`、`vcf_path`、FASTQ 路径、`sample_id`、`driver_genes`、`min_impact`、`auto_stage23`、`sarek.*` |
| `target_to_hit` | `external_activities_csv`、`screening_library_csv`、`docking_only`、`chembl.*`、`model.*`、`filter.*` |
| `hit_to_lead` | `analog_gen.*`、`mpo.*`、`reinvent4.*` |
| `lead_optimization` | `pdb_id`、`protein_sequence`、`docking.*`、`admet_deep.*`（含 `hard_filter.*`）、`md_simulation.*`（含 `ensemble.*`）、`scoring.*`、`explicit_md.*` |

GPU 相关示例见：[gpu_execution_zh.md](gpu_execution_zh.md)。

## 阶段四评分默认权重

- `scoring.w_docking`：0.25  
- `scoring.w_md_energy`：0.40  
- `scoring.w_stability`：0.35  
- `scoring.w_cyp_soft`：0.08（设为 `0` 可关闭）

ADMET 安全闸门：`lead_optimization.admet_deep.hard_filter`（`drop_herg`、`drop_veber_fail`、`drop_high_sa`）。

## 集成 MM-GBSA

在 `lead_optimization.md_simulation.ensemble` 下：

- `enabled`、`n_runs`、`equilibration_ps`、`production_ps`、`sample_interval_ps`

结合能统计为**近似值**（固定受体/配体单结构能量与时间演化的复合物能量相减）。

## AutoDL / 磁盘

若使用 `./data` 且存在 `/autodl-fs/data`，`get_out_dir()` 可能将默认输出重定向到 `/autodl-fs/data/T2Lead`。需要完全可控时请显式设置 `DP_PIPELINE__OUT_DIR`。
