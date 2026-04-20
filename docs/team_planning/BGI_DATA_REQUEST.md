# 向华大基因（BGI）的统一数据需求

> **单一数据入口**：所有下游模块**不单独与 BGI 对接**。此文档由 BioDriver 负责人 + Lead 统一向 BGI 提出。

## 1) 必须交付（Must Have）

### A. 单细胞转录组（scRNA-seq）
- **计数矩阵**：`.h5ad`（首选）或 `.rds`（Seurat）或 10x `.h5`
- **细胞级元数据**（cell-level metadata），至少包含：
  - `cell_barcode`, `sample_id`, `patient_id`
  - `condition`（肿瘤 / 正常对照）
  - `batch`, `tissue/source`
- **基因注释**：`gene_id`（Ensembl）、`gene_symbol`（HGNC 格式优先）
- **质控结果**：每样本测序深度、有效细胞数、线粒体比例分布、双细胞处理说明

### B. 配对全外显子组测序（WES，肿瘤/正常配对）
- **格式**：FASTQ（或已比对好的 BAM + 索引）
- **样本要求**：每个病人 **肿瘤 + 正常对照** 双样本（必须配对，否则没法做体细胞变异）
- **参考基因组版本**：建议 GRCh38（hg38），带 ALT-aware contigs

### C. 样本与临床最小元数据（脱敏）
- `patient_id`, `age_group`, `sex`, `cancer_type`, `stage`（若有）

## 2) 强烈建议（Strongly Recommended）

- 上游流程说明：软件版本、参考基因组、比对/定量参数
- 初步细胞类型注释与聚类结果（作为参考，不直接当最终结论）
- BGI 内部 DEA 结果（如有）：`gene, log2FC, p-value, adj_p`（作为参考，BioDriver 仍会复算）

## 3) 可选（Nice to Have）

- WGS（全基因组）：如果预算允许。对 ADC / 疫苗 / RNAi 这几条路**WES 通常够用**；WGS 主要在"非编码区驱动突变"场景有优势。
- 甲基化 / ATAC-seq 数据（作为二级证据，不是第一阶段必需）
- 治疗史、响应状态（用于后续分层分析）

## 4) 不需要从 BGI 拿（内部可算）

- **HLA 分型**：可用 OptiType / HLA-HD / arcasHLA 从 WES 直接算，由 BioDriver 负责。
- **新抗原候选**：BioDriver 从变异 + HLA 分型内部产出，不要求 BGI 提供。
- **药物-靶点打分**：纯内部建模。

## 5) 关键澄清问题（给 BGI）

1. scRNA-seq 的 FASTQ 与计数矩阵是否同一时间可获取？
2. 肿瘤/正常配对的 WES 是否严格配对？（同一患者）
3. 数据脱敏合规性：能否提供脱敏后的临床最小元数据？
4. 数据交付方式：S3 / FTP / 物理盘？
5. 数据规模：样本数、每样本原始数据体积？

## 6) 对团队内部的使用约束

- 任何模块**不得**绕过 BioDriver 单独对接 BGI。
- 任何需要新数据的请求，统一提给 BioDriver 负责人汇总。
- BioDriver 负责维护"**BGI Delivery Log**"（数据到货时间、版本、QC 结果）。
