# 项目总览与模块职责总表

## 项目愿景（一句话）

> 从病人的组学数据出发，走四条互补的药物发现路径（小分子老药新用 / RNA 干扰 / ADC / mRNA 疫苗），
> 最后统一用物理仿真做质量把关，输出可递交临床合作方（广科医）做湿实验验证的候选清单。

## 六大模块一览

| #    | 模块                                    | 负责人                  | 定位                                                         | 核心交付                                         |
| ---- | --------------------------------------- | ----------------------- | ------------------------------------------------------------ | ------------------------------------------------ |
| 1    | **BioDriver**（生物驱动/发现层）        | 待定（计算生物学/生信） | 唯一对接 BGI 原始数据的**上游中台**，产出标准化"知识产品"分发给 4 下游 | 疾病签名、致病靶点、新抗原、HLA 分型、胞外域结构 |
| 2    | **DrugReflector**（表型药物/老药新用）  | 姜可盈（暂定）          | **L1** 用单细胞签名做小分子重定位；**L2**（阶段二）对弱 hit 触发类似物生成（沿用 T2Lead 资产） | Top-N 小分子候选 + 打分；（阶段二）类似物批次    |
| 3    | **RNAi-Therapy**（RNA 干扰）            | 待定（核酸药物方向）    | 针对致病基因设计 siRNA                                       | siRNA 候选序列 + 脱靶评估                        |
| 4    | **PrecisionDelivery**（ADC 精密递送）   | 待定（蛋白质工程方向）  | 针对肿瘤表面蛋白设计 ADC                                     | 抗体序列 + linker + payload 方案                 |
| 5    | **ImmunoGen**（mRNA 疫苗）              | 待定（免疫/AI 序列）    | 基于新抗原设计 mRNA 多价疫苗                                 | 候选 mRNA 全序列 + 免疫原性评分                  |
| 6    | **Simulation Hub**（物理仿真/验证中台 + L2 触发器） | 卢振涛（暂定）          | 对 4 下游的分子做统一 MD/能量/稳定性验证；对小分子弱 hit 触发 L2 类似物生成 | ΔG、RMSD 稳定性报告、通过/拒绝标签、L2 触发事件  |

## 上下游关系（总线图）

```
BGI ──► BioDriver ──┬──► DrugReflector (L1) ──┐
                    ├──► RNAi-Therapy ─────────┤
                    ├──► PrecisionDelivery ────┼──► Simulation Hub ──► 广科医（湿实验）
                    └──► ImmunoGen ────────────┘              ▲ │
                                                              │ │ (仅小分子弱 hit)
                                                              │ ▼
                                                    Analog Generation (L2，阶段二)
```

**核心原则**：

- BGI 的数据**只对接 BioDriver 一个出口**，避免 4 个下游各自对接华大造成混乱与重复工作。
- 4 个下游**只消费 BioDriver 输出**（除非特殊情况见各自任务书）。
- Simulation Hub 对外按 **`molecule_type` 分派契约**：小分子走 `receptor.pdb + ligand.sdf`；多肽/抗体/siRNA 走 `complex.pdb` 多链方案（**不能用 SDF**，否则力场会崩）。**字段级唯一真相源**：[`SIMHUB_CONTRACT.md`](./SIMHUB_CONTRACT.md)。

详细数据流、字段级接口与 Mermaid 图：见 [`ARCHITECTURE_AND_DATAFLOW.md`](./ARCHITECTURE_AND_DATAFLOW.md)。

## 第一阶段目标（必达，聚焦）

**不是 6 个模块同时启动**。第一阶段只启动：

- BioDriver（必须，作为上游中台）
- DrugReflector（姜可盈暂定）
- **1 个新人二选一的下游**（RNAi / ADC / ImmunoGen 任选一个）
- Simulation Hub（Lead）

目标：

1. 选中的这 4 个模块都产出至少 **1 个真实可跑的 smoke sample**。
2. 模块间用"文件 + 约定字段"交接（不要强行上 RPC / 消息队列）。
3. 联通性测试：`BioDriver → DrugReflector + 选中下游 → Simulation Hub` 端到端跑通一次。
4. 每个模块交一份 `REPORT.md`：输入、工具链、输出、风险与限制。

## 第二阶段目标（扩张）

- 启动剩下的 2 个下游模块（补齐 4 条药物发现线）。
- **激活 DrugReflector 的 Level 2（类似物生成）**：沿用 T2Lead 现有资产（骨架感知生成 + pIC50 ML + ADMET + MPO），由 Simulation Hub 对弱 hit 自动触发，形成 L1↔L2 闭环。
- 多病人批量跑。
- 广科医"决策级候选档案"统一 dashboard。


## 关于 BGI 数据

**统一需求在 `BGI_DATA_REQUEST.md`，由 BioDriver 负责人 + Lead 一起发送。**
各下游模块在自己任务书里**不重复提**，只说明"从 BioDriver 接什么"。

## 关于 Git 协作

见 [`GIT_TEAM_STRATEGY.md`](./GIT_TEAM_STRATEGY.md)。  
推荐：**GitHub Organization + 多仓库**（每模块一仓库），成员只能写自己仓库，修改其他的要 PR。