# 团队任务分工（Team Planning）

> 本目录是“**团队项目管理文档**”，不是研究笔记。面向 6 大模块负责人与新人。
> 所有文件**不进入 git 主线**，按需自行同步。

## 目录结构

- [`MODULES_OVERVIEW.md`](./MODULES_OVERVIEW.md) — 项目总览与 6 大模块职责总表
- [`ARCHITECTURE_AND_DATAFLOW.md`](./ARCHITECTURE_AND_DATAFLOW.md) — 架构图 + 数据流图（Mermaid）
- [`SIMHUB_CONTRACT.md`](./SIMHUB_CONTRACT.md) — **SimHub 输入契约单一真相源**（下游模块必读）
- [`BGI_DATA_REQUEST.md`](./BGI_DATA_REQUEST.md) — 对华大基因（BGI）的统一数据需求（单一来源）
- [`GIT_TEAM_STRATEGY.md`](./GIT_TEAM_STRATEGY.md) — Git 团队协作策略（单仓 vs 多仓 + 权限方案）
- [`modules/`](./modules) — 各模块任务书（抛砖引玉版，细节由负责人自主探索）
  - `BioDriver.md`
  - `DrugReflector_姜可盈.md`
  - `RNAi_Therapy.md`
  - `PrecisionDelivery_ADC.md`
  - `ImmunoGen_mRNA_Vaccine.md`
  - `SimulationHub_Lead.md`

## 使用方式

1. 每个负责人先读 `MODULES_OVERVIEW.md` 建立全局观。
2. 再读 `ARCHITECTURE_AND_DATAFLOW.md` 理解自己在数据流中的位置。
3. 然后读自己的 `modules/<module>.md` 任务书。
4. 向 BGI 提数据需求时，**统一参照** `BGI_DATA_REQUEST.md`（由 BioDriver 负责人与 Lead 一起对接）。

## 一阶段目标（首要）

**把 6 大模块联通并都能跑通（端到端 smoke run）**，不追求每个模块都产业级完美。
