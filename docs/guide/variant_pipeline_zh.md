# 变异驱动工作流（核心路径）

[English (variant_pipeline.md)](variant_pipeline.md)

## 概览

1. 获得**肿瘤/正常**比对与体细胞变异（例如 **nf-core/sarek** + Mutect2 + VEP），或直接提供**已带 VEP 注释的 VCF**。
2. 在配置中启用 **`variant_analysis`**，并设置 `vcf_path`**或**四条 FASTQ 路径。
3. 运行流水线：解析 VCF → 构建突变序列 → 解析结构（PDB / ESMFold）→ 按变异运行**阶段四**。

## 先导分子来源（二选一）

| 模式 | 配置 | 说明 |
|------|------|------|
| **自动 Stage 2-3**（推荐） | `auto_stage23: true` | 每个突变基因自动查 ChEMBL 靶点 → 筛选化合物 → 生成先导 → 对接突变结构。**无需预先准备任何 CSV。** |
| **复用已有先导** | `auto_stage23: false` | 使用之前疾病路径生成的 `stage3_leads/final_lead_candidates.csv`。快，但分子不一定针对该突变靶点。 |

## 配置示例

在 `configs/default_config.yaml`（或自定义 YAML，如 `configs/variant_breast.yaml`）中：

```yaml
pipeline:
  stages:
    - lead_optimization

variant_analysis:
  enabled: true
  vcf_path: "data/user_inputs/vcf/sample.ann.vcf"   # 或留空并配置 FASTQ
  sample_id: "patient_001"   # 可选但推荐，多病人时用于输出隔离
  driver_genes: []        # 空表示全部；或 ["EGFR", "PIK3CA", ...]
  min_impact: "MODERATE"
  auto_stage23: true      # 每个突变基因自动运行 Stage 2-3
```

## 命令行方式（无需改文件）

```bash
# 提供 VCF + 自动 Stage 2-3
python scripts/run_pipeline.py \
  --disease "breast cancer" \
  --vcf-path /path/to/your.vcf.gz \
  --auto-stage23 \
  --stages lead_optimization \
  -v

# 仅指定配置文件
python scripts/run_pipeline.py -c configs/variant_breast.yaml -v
```

## 主路径一键脚本（推荐）

多病人隔离运行且不想写长命令时，直接用：

```bash
python scripts/run_mainpath.py \
  --disease "breast cancer" \
  --sample-id patient_001 \
  --vcf-path /path/to/patient_001.ann.vcf.gz \
  -v
```

`run_mainpath.py` 默认值：
- 输出根目录：`/root/autodl-fs/T2Lead_mainpath`
- `auto_stage23: true`
- Sarek profile：`singularity`
- 默认关闭 ChEMBL 共享缓存（优先隔离）

可用 CLI 参数：

| 参数 | 说明 |
|------|------|
| `--disease "名称"` | 设置疾病名（决定输出子目录） |
| `--vcf-path /路径` | 自动启用 variant_analysis 并设置 VCF 路径 |
| `--auto-stage23` | 启用每个突变基因的自动 Stage 2-3 |
| `--stages s1 s2 ...` | 覆盖运行阶段 |

## 自动衔接

设置 `vcf_path` 时，Python 流水线直接解析该文件。若配置 FASTQ 且 `vcf_path` 为空，`SarekRunner` 会调用 Nextflow；得到的 VCF 路径随后被解析——除非单独跑 sarek，否则不必手动拷贝到固定 `results/sarek/...` 路径。

重点：
- **samplesheet 只用于 FASTQ 路线（Sarek）**。
- **VCF 路线不需要 samplesheet**，直接指向一个注释后的 VCF 文件。
- 做完 VEP 后它仍然是 VCF（`.vcf` / `.vcf.gz`），只是 INFO 里带 `CSQ` 注释字段。

## 独立运行 sarek（Makefile）

```bash
make install-sarek
make run-sarek INPUT=data/user_inputs/samplesheet.csv
```

将 `variant_analysis.vcf_path` 指向 Nextflow `outdir` 产出的**已注释** VCF（具体文件名见 sarek 文档）。

## 按突变的输出

默认（`variant_isolated_runs: true`）会先按运行隔离到：

- `<pipeline.out_dir>/<disease_slug>/variant_runs/<sample_id>/<run_id>/...`

随后在每个突变目录 `stage4_optimization/<GENE_MUTATION>/` 下：

- `optimized_leads.csv`
- 受体 PDB/PDBQT、`docking_poses/`、可选 `md_trajectories/`

当 `auto_stage23: true` 时，每个突变目录下还会包含：

- `stage2_hits/` — 该基因靶点的 ChEMBL 筛选结果
- `stage3_leads/` — 该基因靶点的先导分子

突变结构由编排器写入 `lead_optimization._mutant_pdb_path` 并复制到阶段四工作目录。

## 参考

- [somatic_variant_pipeline_feasibility_zh.md](../research/somatic_variant_pipeline_feasibility_zh.md)
- [reproduction_steps_zh.md](../reproducibility/reproduction_steps_zh.md)
