# 实验记录规范

该目录是每次实验运行记录的统一归档位置。

## 目录结构

每次运行都使用如下路径：

`docs/experiments/runs/<date>_<run-name>_<run-id>/`

示例：

`docs/experiments/runs/2026-04-13_xena-full-max-variants_20260413_191441/`

## 每次运行必须包含的文件

1. `REPORT.md`
   - 面向人的最终实验报告。
   - 必须包含：目标、配置、运行状态、成功点、失败点、根因、产物、下一步。

2. `ARTIFACTS.md`
   - 面向机器/运维的产物索引。
   - 必须包含：绝对路径、文件用途、完整/部分状态。

3. `COMMANDS.md`
   - 实际使用的启动命令、环境变量、监控命令。

4. `METADATA.yaml`
   - 结构化元数据（run_id、病种、sample_id、输入来源、起止时间、状态、备注等）。

> 建议由 Agent 在每次运行后自动创建/更新上述文件，不需要手工维护。

## 语言建议

- 团队内部默认中文：补充 `REPORT_zh.md`、`ARTIFACTS_zh.md`、`COMMANDS_zh.md`。
- 对外协作再补英文版本。

## 推荐状态值

- `completed`（完成）
- `interrupted`（中断）
- `failed`（失败）
- `partial-success`（部分成功）

## 编写规则

- 所有关键产物都使用绝对路径。
- 严格区分“观察事实”和“原因推断”。
- 明确列出尚未解决的阻塞项。
- 若为中断任务，必须写明最后日志时间和最后确认成功步骤。
