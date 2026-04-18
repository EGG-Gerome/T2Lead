# BRCA 全突变 Xena 运行状态与产物复盘

状态：已中断（未完成）  
日期：2026-04-14  
Run ID：`20260413_191441`  
主日志：`/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/logs/20260413_191442_breast_cancer_full.log`

## 执行摘要

本次任务采用全突变 Xena TSV 输入，单卡设置（`parallel_workers=1`），并使用“实验结构优先 + 点突变建模 + 局部优化，失败再回退 ESMFold”的策略。任务启动正常，且生成了较多变异分析产物，但在进入完整流水线收敛前停止。后续检查时未发现活跃的 `run_mainpath` 进程。

置信度说明：关于“进程为何停止”属于高概率推断，不是已证实根因。原因是可见日志中没有明确的终止信号（如 kill 信息）或 Python Traceback。

## 成功点

- Xena TSV 自动转换正常，生成了 VEP 风格输入 VCF。
- 大规模突变 FASTA 生成成功（数量达到数千级）。
- 多个位点通过实验 PDB 下载成功进入结构解析流程。
- 多个位点成功完成点突变建模 + 局部约束优化（产出 `*_mut_localopt.pdb`）。
- 共享缓存策略生效（ChEMBL/report cache 指向 shared 路径）。

## 主要问题与可能原因

1. **任务未完成 / 进程不在运行**
   - 现象：日志在 `2026-04-13 22:13:35 +0800` 后停止更新，且没有运行中的主进程。
   - 可能原因：前台会话中断（SSH 断连、终端关闭、人工中断、系统级终止）。

2. **RCSB 查询偶发失败（HTTP 400）**
   - 现象：`POST /rcsbsearch/v2/query` 间歇性返回 `400`。
   - 可能原因：特定查询在服务端校验失败或 API 侧短暂异常。
   - 影响：会触发回退；多数位点仍可通过候选替代条目继续。

3. **长序列 ESMFold 被拒（HTTP 413）**
   - 现象：`POST /foldSequence/v1/pdb/` 返回 `413 Request Entity Too Large`。
   - 可能原因：外部 ESMFold API 对请求体/序列长度存在限制。
   - 影响：此类位点需回退其他来源；若都不可用则结构无法解析。

4. **部分 UniProt 在 AlphaFold DB 不存在（HTTP 404）**
   - 现象：`AF-<UniProt>-F1-model_v4.pdb` 下载不到。
   - 可能原因：该 accession/version 无公开 AlphaFold DB 条目。
   - 影响：对应位点可能出现 `No structure resolved`。

5. **部分候选 PDB ID 不可下载（HTTP 404）**
   - 现象：如 `8C6J.pdb` 不存在，后切换 `6ID1.pdb` 成功。
   - 可能原因：条目废弃、更新或下载路径不可用。
   - 影响：主要是重试和切换导致延迟，通常不致命。

## 本次尝试效果评估

- 新的结构桥接策略在真实负载下可工作（不仅是单测/smoke），已有多例 `*_mut_localopt.pdb` 产物验证。
- 主要失败模式来自外部数据覆盖/API 限制，不是核心流程逻辑整体失效。
- 但本次未跑完整条下游流程，因此 4k+ 变异规模的端到端 Stage2/3/4 最终结果尚未形成。

## 关键生成文件及含义

### 运行根目录
- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/`
  - 本次尝试的专属输出根目录。

### 日志
- `.../logs/20260413_191442_breast_cancer_full.log`
  - 详细运行日志（debug/info/warning）。
- `.../logs/20260413_191442_breast_cancer_summary.log`
  - 摘要日志。

### 输入转换产物
- `.../variant_analysis/converted_inputs/xena_converted.TCGA-AN-A046-01A.vep.vcf`
  - Xena TSV 转换后的 VEP 风格 VCF，供后续变异流程使用。

### 突变序列产物
- `.../variant_analysis/mutant_fastas/`
  - 每个变异对应的突变 FASTA 文件（本次快照共 4638 个）。

### 结构产物
- `.../variant_analysis/structures/`
  - 实验模板结构（如 `6e1f.pdb`, `9ejy.pdb`）。
  - 点突变局部优化结构（如 `6e1f_R251Q_mut_localopt.pdb`）。
  - 部分 ESMFold 结构（如 `AGMAT_E320K_esmfold.pdb`）。

## 下次运行建议

1. 用 `tmux` / `screen` / `nohup` 启动长任务，避免 SSH 断开导致中断。
2. 当前单卡机器维持 `parallel_workers=1`。
3. 对超大变异集建议分批执行（按驱动基因或分块），提高可控性和可恢复性。
4. 每次运行都在 `docs/experiments/runs/<date>_<run-name>_<run-id>/` 归档。
