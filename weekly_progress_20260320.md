# T2Lead 周报（2026-03-14 ~ 2026-03-20）

## 上周回顾

上周（Week 06）完成了 **IsoDDE-BC-Dify MVP**：基于 Dify 平台构建了乳腺癌个体化用药推荐原型，实现了「基因突变 → 蛋白序列 → IsoDDE 亲和力预测 → 临床权重排序」的概念验证流程。该方案依赖 IsoDDE 单点预测 + 手动临床权重表，覆盖约 50-100 个已上市药物，主要面向"已有突变患者"的选药场景。

## 本周工作概述

本周将项目从"患者选药"的下游应用切换到"药物发现"的上游引擎——**T2Lead 端到端自动化流水线**的开发，从 Version 1 基础框架推进到四阶段全链路可跑通状态。相比上周 IsoDDE MVP 的手动脚本串联，本周实现了完全自动化的 Target → Hit → Lead → Optimization 全流程，从输入一个疾病名称到输出可优化先导化合物，只需一条命令。

---

## 一、核心进展

### 1. 四阶段流水线全链路跑通

上周的 IsoDDE MVP 需要手动执行 12 个步骤、在本地和 AutoDL 之间反复传输文件。本周 T2Lead 实现了**单命令全自动运行**：`make run DISEASE="breast cancer"` 即可完成从靶点发现到先导优化的全部四个阶段，总耗时约 35 分钟。

**Stage 1 — 靶点发现**：接入 OpenTargets GraphQL API，自动解析疾病关联靶点并按综合评分排序。修复了 API 字段映射（`proteinAnnotations` → `dbXrefs`）和断点续跑逻辑，确保 Stage 标记与缓存状态正确判断已完成步骤。

**Stage 2 — 虚拟筛选（Target to Hit）**：改进了靶点选择策略——不再盲选评分最高的靶点，而是先查询各候选靶点在 ChEMBL 的实际 IC50 数据量，选择满足最低训练阈值（≥ 200 条）的靶点。本次乳腺癌场景中，BRCA1 虽然评分最高但仅有 20 条 IC50 记录，流水线自动选择了拥有 9728 条记录的 **PIK3CA**，保证了模型训练的数据充分性。虚拟筛选环节从 CREM 数据库的 285 万分子中经过 pIC50 预测、QED、结构警报和 ADMET 规则层层过滤，最终筛出 200 个 Hit。

**Stage 3 — Hit to Lead**：对 100 个 Hit 进行骨架分析（28 个独立骨架）、Butina 多样性聚类（31 个簇），并通过 CREM 生成 4672 个类似物。经 MPO 多目标评分排序后输出 50 个 Lead 候选。同时集成了 REINVENT4 强化学习分子生成桥接模块（详见第 4 点）。

**Stage 4 — 先导优化**：这是本周从零完成的全新阶段。包含蛋白准备（PDB 下载 → PDBFixer 修复 → PDBQT 转换）、分子对接（AutoDock Vina）、增强 ADMET 评估和 MM-GBSA 结合能计算，最终通过 composite score 综合排序，输出 10 个最优先导化合物。

### 2. 大规模计算性能优化

流水线处理 285 万分子的计算瓶颈本周得到显著改善：

- **分子属性计算并行化**：将 QED、描述符、PAINS 检测改为多进程并行（自动检测 CPU 核数，上限 16 workers），加速约 **10-15 倍**。
- **分子对接并行化**：从逐个分子串行对接改为 `ProcessPoolExecutor` 并行，208 核机器同时对接 6 个分子（每分子分配 32 核 exhaustiveness），CPU 利用率从 ~15% 提升至 **~95%**。50 个分子的对接在 8 分钟内完成。
- **指纹计算优化**：替换已废弃的 `AllChem` API 为 `rdFingerprintGenerator`，新增单例缓存避免重复初始化。
- **PDB 修复加速**：跳过缺失残基环重建步骤，PDBFixer 耗时从 15 分钟降至 2 分钟。

### 3. 全链路缓存系统

为解决迭代开发中反复重跑的时间成本，实现了三级缓存体系：

| 缓存对象 | 存储格式 | 效果 |
|---|---|---|
| Morgan 分子指纹（285 万） | fp_cache/*.npy | 11 分钟 → <1 秒 |
| RF + MLP 训练模型 | model_cache/ (joblib + .pt) | 训练 → 直接加载 |
| 虚拟筛选打分结果 | scored_candidates.csv | 完整筛选 → 跳过 |

首次完整运行约 30 分钟，参数未变时重跑 **< 10 秒**。

### 4. REINVENT4 部署与集成

本地部署了 AstraZeneca 开源的 REINVENT4 分子生成框架（`/root/REINVENT4`），并修改其设备管理逻辑以同时支持 **CPU / CUDA / MPS** 三种后端，确保在不同硬件环境下均可运行。在 T2Lead 中编写了 `reinvent_bridge.py` 桥接模块：自动生成 REINVENT4 配置文件、序列化 pIC50 预测模型为评分脚本、以子进程启动 RL 优化并收集生成的 SMILES。目前 REINVENT4 在当前环境下启动后崩溃（见"待解决问题"），桥接模块已做了 graceful fallback，不影响主流程。

### 5. GPU 兼容性与环境适配

服务器配置 RTX 5090（Blackwell 架构，sm_120），对 CUDA 生态提出了较高的兼容性要求：

- PyTorch 升级至 `2.10.0+cu128`，确保 Blackwell GPU 的原生支持。
- 为 REINVENT4 入口脚本自动注入 `LD_LIBRARY_PATH`，解决 libstdc++ 动态库加载问题。
- README 新增 RTX 50/40/30 系列 CUDA 版本对照表和具体安装命令。

### 6. 多场景入口与 ESMFold 集成

在标准的"疾病名 → 全自动"流程之外，新增了两种使用场景：

- **场景 B**（`--activities-csv`）：导入自有 IC50 数据，跳过 ChEMBL 爬取，适用于已有文献或实验室数据的情况。
- **场景 C**（`--docking-only`）：纯对接模式，跳过 ML 训练，适用于零活性数据的全新靶点。
- 新增 `--protein-sequence` 参数 + **ESMFold API 集成**：只需氨基酸序列即可自动预测蛋白 3D 结构，无需 PDB ID 或本地 GPU，实现 `序列 → 结构 → 修复 → 对接` 的全自动闭环。

### 7. 工程化与代码质量

- **Makefile**：新增 `install` / `run` / `run-docking` / `clean` / `clean-all` / `clean-logs` 等标准化目标，加入 `--no-capture-output` 以支持运行时日志实时流式输出。
- **双日志系统**：每次运行自动生成 `*_full.log`（完整 DEBUG 级别）和 `*_summary.log`（仅关键节点、指标和错误警告），方便调试和快速回顾。
- **商业级中英双语注释**：对 `scripts/` 和 `src/` 下所有核心模块补齐面向维护的业务语义注释，统一模块级与函数级说明格式，为后续团队协作和代码审查做好准备。
- **模型评估增强**：新增 5-fold 交叉验证，报告 per-fold RMSE/R² 及 mean±std，比单一 test set 指标更可靠。

---

## 二、本周运行结果

以乳腺癌（breast cancer）为测试场景，最近一次完整运行（2026-03-20 14:05 ~ 14:40）的结果：

| Stage | 状态 | 关键指标 |
|---|---|---|
| Stage 1 — 靶点发现 | 成功 | 5 个候选靶点，自动选中 PIK3CA（CHEMBL4005，9728 条 IC50） |
| Stage 2 — 虚拟筛选 | 成功 | 7968 分子训练集；RF R²=0.771，RMSE=0.600；5-fold CV R²=0.730±0.019；200 个 Hit |
| Stage 3 — Hit to Lead | 成功 | 28 个骨架、31 个 Butina 簇；4672 个类似物；50 个 Lead（REINVENT4 模块失败，已 fallback） |
| Stage 4 — 先导优化 | 部分成功 | 50/50 对接通过（全部优于 −6.0 kcal/mol）；0/50 ADMET 高风险；MM-GBSA 因 CUDA 兼容性问题未能计算（见下文） |

**最终输出 10 个先导化合物**，Top 1 综合评分 0.969：

```
SMILES:   COc1ncc(-c2cc3c(C)nc(N)nc3n([C@H]3CC[C@H](C(C)=O)CC3)c2=O)cc1F
pIC50:    8.61（预测活性 IC50 ≈ 2.5 nM）
QED:      0.683
Docking:  −9.869 kcal/mol
MPO:      0.879
ADMET:    0.086（低风险）
```

---

## 三、与上周对比

| 维度 | Week 06 — IsoDDE-BC-Dify MVP | Week 07 — T2Lead Pipeline |
|---|---|---|
| 定位 | 患者个体化选药（已有药物重排序） | 药物发现（从疾病名到新化合物） |
| 自动化程度 | 12 个手动步骤，需本地/云端反复传文件 | **单命令全自动**，一条 `make run` 跑完 |
| 分子空间 | 50-100 个已上市药物 | **285 万** CREM 化合物库 + 4672 类似物生成 |
| 筛选能力 | IsoDDE 单点亲和力预测 | RF + MLP 模型预测 + QED/PAINS/ADMET 多层过滤 + Vina 对接验证 |
| 验证手段 | 临床指南权重打分 | 分子对接（Vina）+ MM-GBSA 结合能 + ADMET 深度评估 |
| 输出质量 | 药物排序表（重排已知药物） | 10 个全新先导化合物（含 SMILES、预测活性、对接构象） |
| 模型评估 | 无（直接使用 IsoDDE 输出） | 5-fold CV，R²=0.730±0.019 |
| 计算时间 | 分钟级（药物数少） | ~35 分钟首次运行，<10 秒缓存重跑 |

两个项目定位互补：IsoDDE MVP 面向临床端"有突变、选已有药"，T2Lead 面向研发端"有疾病、找新分子"。

---

## 四、待解决问题

### 1. MM-GBSA 计算失败

OpenMM 在 CUDA 模式下报 `CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`，原因是当前 OpenMM 版本编译的 PTX 不支持 Blackwell 架构（sm_120）。10 个先导化合物的 MD 结合能和构象稳定性均未能计算，composite score 已自动将 MD 权重重分配给其他维度，但评分完整性有所降低。

**计划**：优先尝试增加 CPU fallback 模式，同时关注 OpenMM 新版本对 CUDA 12.8 的支持。

### 2. REINVENT4 启动崩溃

REINVENT4 subprocess 返回 rc=1，traceback 在日志中被截断。目前桥接模块已做 graceful fallback（返回空 DataFrame），不影响主流程，但基于强化学习的分子生成优化功能处于不可用状态。

**计划**：检查完整 stderr 输出，排查 Python 版本或 libstdc++ 兼容性问题。

### 3. 其他待补齐项

- ESMFold 集成已实现但未在真实运行中测试（本次使用了已有的 4JPS PDB 结构）。
- 场景 B / C 的 `--activities-csv` 和 `--docking-only` 模式尚无端到端验证。
- `md_trajectories/` 目录已建立但为空，轨迹分析功能待 MM-GBSA 跑通后补全。
- 对接构象和分子结构的可视化脚本尚未编写。
- 单元测试覆盖率为零，需逐步补充。

---

## 五、下周计划

1. 解决 MM-GBSA CUDA 兼容性问题（CPU fallback 或升级 OpenMM），补全先导化合物的结合能数据。
2. 调试 REINVENT4 崩溃，使 RL 分子生成功能可用。
3. 对 `--docking-only` 和 `--activities-csv` 场景各跑一个端到端测试。
4. 验证 composite score 在 MD 数据完整时的权重分配逻辑。
