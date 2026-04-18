# 输出文件说明

[English (output_reference.md)](output_reference.md)

假定运行根目录为：
- 标准路径：设置疾病时为 `<pipeline.out_dir>/<disease_slug>/`，否则为 `<pipeline.out_dir>/`。
- 变异路径（`variant_analysis.enabled: true` 且 `pipeline.output_layout.variant_isolated_runs: true`，默认开启）：`<pipeline.out_dir>/<disease_slug>/variant_runs/<sample_id>/<run_id>/`。

其中 `sample_id` 优先取 `variant_analysis.sample_id`（为空时从 VCF/FASTQ 文件名推断），`run_id` 优先取 `pipeline.output_layout.variant_run_id`（为空时自动使用启动时间戳）。若 `use_stage_subdirs: false`，各阶段键将落在同一目录。

## 共享 ChEMBL 库（可选）

当配置了 `target_to_hit.shared_library_dir`（默认见 `default_config.yaml`）时，下列内容位于 `<pipeline.out_dir>/<shared_library_dir>/`，供多种疾病共用同一次爬取与同一份 Morgan `fp_cache/`：

| 文件 / 目录 | 说明 |
|-------------|------|
| `molecules_chemblid_smiles.csv` | 爬取的分子 |
| `activities_ic50.csv` | 爬取的 IC50 |
| `crawl_state.json` | 断点续跑 |
| `fp_cache/` | Morgan 指纹缓存 |

各疾病的 `stage2_hits/` 仍保存 `dataset_target_ic50.csv`、`model_cache/`、`scored_candidates.csv`、`final_hit_candidates.csv` 等。若 `shared_library_dir` 为空，则恢复旧行为：上述文件全部在各自的 `stage2_hits/` 下。

## 阶段一 — `stage1_targets/`

| 文件 | 说明 |
|------|------|
| `ranked_targets.csv` | OpenTargets 排序后的靶点 |

## 阶段二 — `stage2_hits/`

| 文件 | 说明 |
|------|------|
| `molecules_chemblid_smiles.csv` | ChEMBL 模式爬取的分子（若启用共享库则见上文 **共享 ChEMBL 库**） |
| `activities_ic50.csv` | 爬取的 IC50 行（同上） |
| `crawl_state.json` | 断点续跑检查点（同上） |
| `dataset_target_ic50.csv` | 训练用表 |
| `model_cache/` | 序列化的 RF / MLP |
| `fp_cache/` | Morgan 指纹缓存（同上；否则在本目录下） |
| `scored_candidates.csv` | 虚拟筛选输出 |
| `final_hit_candidates.csv` | ADMET / 效力过滤后 |

## 阶段三 — `stage3_leads/`

| 文件 | 说明 |
|------|------|
| `h2l_scaffold_summary.csv` | 骨架 SAR |
| `h2l_cluster_summary.csv` | 聚类摘要 |
| `final_lead_candidates.csv` | MPO 排序后的先导 |
| `reinvent4_*` | 可选 RL 产物 |

## 阶段四 — `stage4_optimization/`

| 文件 | 说明 |
|------|------|
| `<PDB>.pdb`、`*_fixed.pdb`、`*_receptor.pdbqt` | 受体 |
| `docking_poses/pose_*.pdbqt` | 对接构象 |
| `md_trajectories/` | MD 日志（若启用） |
| `visual_assets/` | 可选的二维结构导出（lead/benchmark 网格图与单分子 SVG/PNG） |
| `optimized_leads.csv` | 最终表，含 `docking_score`、`md_binding_energy`、`md_binding_energy_std`、`md_rmsd_mean`、`opt_score` 等 |

仅变异运行时可能使用嵌套目录：`stage4_optimization/EGFR_L858R/optimized_leads.csv`。

## 其他

| 路径 | 说明 |
|------|------|
| `logs/` | 带时间戳的 `*_full.log`、`*_summary.log` |
| `variant_analysis/` | 运行根下 `mutant_fastas/`、`structures/` |
| `variant_analysis/stage4_resume_checkpoint.csv` | 变异级 Stage4 断点续跑状态文件 |
| `patient_aggregation/` | 跨变异病人级汇总（`patient_recommendations.csv`、`recommendation_notes_zh.md`、`dashboard.html`） |
| `variant_calling/` | FASTQ 驱动检出时的 sarek 样本表与输出 |

首次使用者可先看目录图解：
[variant_runs_structure_zh.md](variant_runs_structure_zh.md)。

实用脚本：
- `python scripts/export_stage4_2d.py --disease "breast cancer" --out-dir /root/autodl-fs/T2Lead`
- `python scripts/organize_stage4_assets.py --disease "breast cancer" --out-dir /root/autodl-fs/T2Lead`
