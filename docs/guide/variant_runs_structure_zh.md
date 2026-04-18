# `variant_runs` 目录结构说明（中文）

适用路径：

`/root/autodl-fs/T2Lead_mainpath/<disease>/variant_runs/<sample_id>/<run_id>/`

示例：

`/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/`

## `variant_runs/` 一级子目录含义

`variant_runs/` 下的一级目录通常就是 `sample_id`（也可能是测试标签），例如：

- `xena_full_max_variants`
- `single_min_vcf`
- `patient_demo`

每个一级目录下会有多个 `run_id` 时间戳目录，分别对应不同批次运行。

## 顶层目录含义

- `logs/`
  - 运行日志：`*_full.log`（详细）与 `*_summary.log`（摘要）。
- `variant_analysis/`
  - 变异路径核心中间产物：
  - `converted_inputs/`：Xena TSV 转换后的 VEP 风格 VCF。
  - `mutant_fastas/`：每个变异位点对应一个突变蛋白 FASTA。
  - `structures/`：结构文件（实验 PDB、`*_mut_localopt.pdb`、`*_esmfold.pdb`）。
  - `stage4_resume_checkpoint.csv`：变异级 Stage4 断点续跑状态。
- `stage1_targets/`
  - 靶点发现结果。
- `stage2_hits/`
  - 命中筛选结果。
- `stage3_leads/`
  - 先导分子结果。
- `stage4_optimization/`
  - 对接/MD 优化结果（`optimized_leads.csv`、姿势、轨迹等）。
- `patient_aggregation/`
  - 病人级聚合结果（跨变异汇总后的推荐清单与 dashboard）。

## 常见疑问

### 为什么会有 FASTA？

不是输入 FASTQ。这里的 FASTA 是由 VCF/TSV 解析出的蛋白突变序列中间文件，用于后续结构建模。

### `*_mut_localopt.pdb` 是什么？

表示在实验模板结构上做点突变后，再对突变附近做局部约束优化得到的突变结构。

### 只看到 `variant_analysis`，后续阶段为空？

通常表示任务中断或仍在前置处理中。先看 `logs/*_full.log` 最后时间与最后成功步骤。
