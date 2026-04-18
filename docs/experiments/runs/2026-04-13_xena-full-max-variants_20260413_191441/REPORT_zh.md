# 运行报告：BRCA Xena 全突变（max_variants）

- status: `interrupted`
- run_id: `20260413_191441`
- sample_id: `TCGA-AN-A046-01A`
- strategy: `xena-sample-strategy=max_variants`
- runtime mode: `parallel_workers=1`, `device=cuda`

## 目标

对选定 Xena 样本执行突变驱动流程，并采用最佳实践结构策略：

1. 优先实验 PDB；
2. 执行点突变建模；
3. 进行局部约束优化；
4. 无模板时回退 ESMFold。

## 结果快照

- 输入转换成功。
- 变异分析阶段产出了大规模 FASTA 与结构文件。
- 未达到流水线最终完成态；复查时主进程已不在运行。

## 主要成功点

- Xena TSV -> VEP-VCF 转换成功。
- 多个位点成功产出 `*_mut_localopt.pdb`。
- 共享缓存策略已生效，避免了关键缓存目录重复下载。

## 主要问题与原因

1. **完成前中断**
   - 最后日志时间：`2026-04-13 22:13:35 +0800`。
   - 后续未发现活跃 `run_mainpath` 进程。
   - 高概率为会话/进程中断，而非稳定可复现代码崩溃。

2. **外部结构服务覆盖/限制**
   - RCSB 查询偶发 `400`。
   - 长序列触发 ESMFold `413`。
   - 部分 accession 在 AlphaFold DB 返回 `404`。
   - 影响：部分变异结构解析失败或回退后仍不可得。

## 最后确认成功步骤

`SYF2 K71T` 在 `8C6J` 下载失败后，已成功下载替代模板 `6ID1.pdb`。

## 下一步建议

1. 使用 `tmux` / `screen` / `nohup` 重新启动。
2. 单卡环境保持 `parallel_workers=1`。
3. 全突变建议分批执行，提升可控性与可恢复性。
