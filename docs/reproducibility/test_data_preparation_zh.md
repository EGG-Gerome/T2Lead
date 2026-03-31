# 测试数据准备（WES / WGS）

[English (test_data_preparation.md)](test_data_preparation.md)

T2Lead 的 `variant_analysis` 模块支持两种输入：

1. **成对 FASTQ**（肿瘤与正常样本；当 `vcf_path` 为空时交给 nf-core/sarek），或  
2. **带 VEP 注释的 VCF**，其 INFO 字段需与 `vcf_parser.py` 中的解析逻辑兼容。

## 公开参考数据

- **1000 Genomes / gnomAD / TCGA** — 按期刊与 dbGaP 政策下载 BAM/FASTQ；若仅有 BAM 需转为 FASTQ。  
- **合成 / 最小 VCF** — CI 式检查可手写小 VCF（一条错义突变、PASS 过滤）；VEP 风格 INFO 模式见 `tests/test_variant_analysis.py`。

## 样本表（sarek）

nf-core/sarek 3.x 期望 CSV 列类似：`patient,sex,status,sample,lane,fastq_1,fastq_2`（以官方文档为准）。将 CSV 路径传给 `make run-sarek INPUT=...`。

## 同意书与范围

使用去标识化研究数据；外显子 panel 应与 sarek 配置的**基因组版本**（`GRCh38` vs `GRCh37`）及 VEP 注释一致。
