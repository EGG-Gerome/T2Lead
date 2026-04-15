# 产物索引：20260413_191441

## 运行根目录

- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/`
  - 本次运行的完整输出根目录。

## 日志

- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/logs/20260413_191442_breast_cancer_full.log`
  - 详细日志（debug/info/warn）。
- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/logs/20260413_191442_breast_cancer_summary.log`
  - 摘要日志。

## 输入转换

- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/variant_analysis/converted_inputs/xena_converted.TCGA-AN-A046-01A.vep.vcf`
  - 由 Xena TSV 转换得到的 VEP 风格 VCF。

## 突变 FASTA

- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/variant_analysis/mutant_fastas/`
  - 每个变异对应一个 FASTA。
  - 本次快照观察数量：`4638` 个文件。

## 结构输出

- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/variant_analysis/structures/`
  - 实验模板结构（`*.pdb`，如 `6e1f.pdb`）。
  - 点突变局部优化结构（`*_mut_localopt.pdb`，如 `6e1f_R251Q_mut_localopt.pdb`）。
  - ESMFold 结构（`*_esmfold.pdb`，如 `AGMAT_E320K_esmfold.pdb`）。

## 主路径共享缓存

- `/root/autodl-fs/T2Lead_mainpath/shared/chembl_stage2` -> 软链接到 `/root/autodl-fs/T2Lead/shared/chembl_stage2`
- `/root/autodl-fs/T2Lead_mainpath/shared/report_cache` -> 软链接到 `/root/autodl-fs/T2Lead/shared/report_cache`

## 完成度说明

本次状态为 `interrupted`（未完成）。当前产物集属于中断前的部分结果。
