# 本次使用命令说明

## 启动意图（单 GPU）

环境变量：

- `DP_PIPELINE__DEVICE=cuda`
- `PARALLEL_WORKERS=1`

核心参数：

- `--vcf-path /root/autodl-fs/T2Lead_mainpath/data/test_inputs/xena/TCGA-BRCA.somaticmutation_wxs.tsv`
- `--xena-sample-strategy max_variants`
- `--sample-id xena_full_max_variants`
- 启用共享缓存模式（当前脚本逻辑默认共享 ChEMBL 缓存）

## 监控命令

- 进程检查：
  - `pgrep -af "run_mainpath|python.*drugpipe"`
- 日志检查：
  - 查看 `.../logs/*_full.log`
- 新鲜度检查：
  - `stat .../logs/*_full.log`

## 下次建议启动方式（长任务）

使用 `tmux` 或 `nohup`，避免 SSH 断开导致前台任务终止。
