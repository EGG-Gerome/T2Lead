# 复现检查清单

[English (reproduction_steps.md)](reproduction_steps.md)

用于在干净机器上验证 T2Lead，或作为论文/补充材料中的复现步骤。

## 1. 环境

- [ ] 安装 Conda；在仓库根目录执行 `make install`（或按 `Dockerfile` 构建）。
- [ ] `make test` 全部通过（未安装 OpenMM 时，4 个 OpenMM 相关测试可能跳过）。

## 2. 最小疾病路径

```bash
export DP_PIPELINE__OUT_DIR=/tmp/t2lead_smoke
python scripts/run_pipeline.py --disease "breast cancer" \
  --stages target_discovery target_to_hit \
  -v
```

预期生成 `stage1_targets/ranked_targets.csv` 并填充 `stage2_hits/`（ChEMBL 爬取可能较久；冒烟测试可在 YAML 中降低分页/条数限制）。

## 3. 变异路径（集成）

- [ ] 安装 Nextflow + Java 21；Docker 已启动以运行 sarek。
- [ ] 三选一输入路径：
  - 肿瘤/正常 FASTQ + sarek 样本表（列定义见 nf-core/sarek 文档），
  - VEP 注释 VCF，
  - Xena somatic mutation TSV（流程自动转换）。
- [ ] 执行 `make run-sarek INPUT=...`（FASTQ 路线）**或**通过 `variant_analysis.vcf_path` 提供 VCF/TSV。
- [ ] 确保存在 `stage3_leads/final_lead_candidates.csv`（先跑阶段二、三，或复制含 `canonical_smiles` 的小型测试 CSV）。
- [ ] 设置 `variant_analysis.enabled: true` 并执行 `python scripts/run_pipeline.py -c variant.yaml -v`。
- [ ] 确认出现 `stage4_optimization/<GENE_MUT>/optimized_leads.csv`。

## 4. 合成 VCF / 单元级

- [ ] `python -m pytest tests/test_variant_analysis.py -v`（无需网络）。

## 5. 可选 GPU

- [ ] 使用 CUDA 版 PyTorch 完成 `make install` 后，以 `pipeline.device: cuda` 跑阶段四，并在日志中确认 OpenMM 选用 CUDA。

## 故障排除

- **optimized_leads 为空：** 检查 ADMET 硬过滤是否删光所有行；可暂时放宽 `sa_score_max` 或关闭部分 `hard_filter` 开关做排查。
- **Sarek / VCF 路径：** Makefile 默认将结果写入 `variant_calling/sarek_results/`；请把 `vcf_path` 指向该次运行产出的真实已注释 VCF。
- **Xena TSV 路径：** 若 `vcf_path` 指向 Xena TSV，流程会先转换为 VEP 风格 VCF；请在 `variant_analysis/converted_inputs/` 查转换结果。
- **双输出根目录：** 始终设置 `DP_PIPELINE__OUT_DIR`，使中间产物与最终结果在同一目录树。

全外显子数据准备（公开数据、BQSR、panel 与参考基因组一致等）请遵循机构规范；T2Lead 消费的是 **FASTQ / 注释 VCF / Xena TSV** 输入，但用于变异解析的核心载体是 VCF 语义而非原始读段。
