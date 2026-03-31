# 变异驱动工作流（核心路径）

[English (variant_pipeline.md)](variant_pipeline.md)

## 概览

1. 获得**肿瘤/正常**比对与体细胞变异（例如 **nf-core/sarek** + Mutect2 + VEP），或直接提供**已带 VEP 注释的 VCF**。
2. 在配置中启用 **`variant_analysis`**，并设置 `vcf_path`**或**四条 FASTQ 路径。
3. 准备（或复用）阶段二、三的 **`stage3_leads/final_lead_candidates.csv`**，以便阶段四有待对接分子。
4. 运行流水线：解析 VCF → 构建突变序列 → 解析结构（PDB / ESMFold）→ 按变异在 `stage4_optimization/<GENE_MUTATION>/` 下运行**阶段四**。

## 配置示例

在 `configs/default_config.yaml`（或自定义 YAML）中：

```yaml
pipeline:
  stages:
    - lead_optimization    # 常见：仅阶段四，上游基因组学已完成

variant_analysis:
  enabled: true
  vcf_path: "data/user_inputs/vcf/sample.ann.vcf"   # 或留空并配置 FASTQ
  # tumor_fastq_r1, tumor_fastq_r2, normal_fastq_r1, normal_fastq_r2
  driver_genes: []        # 空表示全部；或 ["EGFR", "PIK3CA", ...]
  min_impact: "MODERATE"
```

**自动衔接：** 设置 `vcf_path` 时，Python 流水线直接解析该文件。若配置 FASTQ 且 `vcf_path` 为空，`SarekRunner` 会调用 Nextflow；得到的 VCF 路径随后被解析——除非单独跑 sarek，否则不必手动拷贝到固定 `results/sarek/...` 路径。

## 独立运行 sarek（Makefile）

```bash
make install-sarek
make run-sarek INPUT=data/user_inputs/samplesheet.csv
```

将 `variant_analysis.vcf_path` 指向 Nextflow `outdir` 产出的**已注释** VCF（具体文件名见 sarek 文档）。

## 仅跑阶段四的前置条件

磁盘上**必须**已有先导分子：

`<out_dir>/<disease_slug>/stage3_leads/final_lead_candidates.csv`

（或对同一疾病/配置先完整跑过阶段二、三。）

## 按突变的输出

在 `stage4_optimization/<GENE_MUTATION>/` 下：

- `optimized_leads.csv`
- 受体 PDB/PDBQT、`docking_poses/`、可选 `md_trajectories/`

突变结构由编排器写入 `lead_optimization._mutant_pdb_path` 并复制到阶段四工作目录。

## 参考

- [somatic_variant_pipeline_feasibility_zh.md](../research/somatic_variant_pipeline_feasibility_zh.md)
- [reproduction_steps_zh.md](../reproducibility/reproduction_steps_zh.md)
