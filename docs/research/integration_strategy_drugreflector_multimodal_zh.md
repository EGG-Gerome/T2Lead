# T2Lead 多模态治疗策略整合方案 — 讨论稿 v0.1

> **状态**: 草稿 / 待多轮讨论确认  
> **日期**: 2026-04-15  
> **范围**: DrugReflector 整合评估、多组学输入升级、小分子/RNA疫苗/放疗三轨决策、外部工具链优化

---

## 目录

1. [现状总结：T2Lead 当前架构](#1-现状总结t2lead-当前架构)
2. [DrugReflector 深度评估](#2-drugreflector-深度评估)
   - 2.1 [它做了什么](#21-它做了什么)
   - 2.2 [它没做什么 — 为什么仍需 H2L/LO](#22-它没做什么--为什么仍需-h2llo)
   - 2.3 [与 T2Lead 互补性分析](#23-与-t2lead-互补性分析)
   - 2.4 [整合方案](#24-整合方案)
3. [多组学输入升级路径](#3-多组学输入升级路径)
4. [外部工具链定位与优化](#4-外部工具链定位与优化)
   - 4.1 [结构预测：ESMFold vs ColabFold vs AF2 vs AF3](#41-结构预测esmfold-vs-colabfold-vs-af2-vs-af3)
   - 4.2 [分子生成与优化：REINVENT4](#42-分子生成与优化reinvent4)
   - 4.3 [分子动力学：OpenMM](#43-分子动力学openmm)
   - 4.4 [推荐工具链组合](#44-推荐工具链组合)
5. [三轨治疗决策管线设计](#5-三轨治疗决策管线设计)
   - 5.1 [Track A：小分子配体（靶向药物）](#51-track-a小分子配体靶向药物)
   - 5.2 [Track B：RNA 疫苗 / mRNA 治疗](#52-track-brna-疫苗--mrna-治疗)
   - 5.3 [Track C：放射治疗（质子/重离子/TRT）](#53-track-c放射治疗质子重离子trt)
   - 5.4 [决策树与分流逻辑](#54-决策树与分流逻辑)
6. [整体架构蓝图（新）](#6-整体架构蓝图新)
7. [与现有代码的兼容性评估](#7-与现有代码的兼容性评估)
8. [待讨论问题清单](#8-待讨论问题清单)
9. [参考文献](#9-参考文献)

---

## 1. 现状总结：T2Lead 当前架构

### 1.1 管线主路径

```
患者数据 (VCF / Xena TSV / FASTQ+Sarek)
    │
    ├─ 变异驱动路径 (main path / variant_analysis)
    │   ├─ VCF 解析 → 体细胞突变序列构建 (MutantSequenceBuilder)
    │   ├─ 结构获取 (StructureBridge: 实验 PDB → AF-DB → ESMFold API)
    │   ├─ 突变建模 + OpenMM 局部最小化
    │   ├─ [可选] 自动 Stage 2-3 (ChEMBL → Hit → CReM/REINVENT4 类似物)
    │   └─ Stage 4: Vina 对接 → DeepADMET → OpenMM MD → 排名
    │       └─ 患者汇聚 (patient_aggregation)
    │
    └─ 疾病驱动路径 (disease path)
        ├─ Stage 1: OpenTargets 靶点发现
        ├─ Stage 2: ChEMBL/ML 虚拟筛选
        ├─ Stage 3: CReM/MPO/REINVENT4 (H2L)
        └─ Stage 4: Vina/ADMET/OpenMM (LO)
```

### 1.2 当前输入限制

| 数据类型 | 当前支持 | 备注 |
|---------|---------|------|
| WES / WGS (VCF) | ✅ 已实现 | variant_analysis 主入口 |
| FASTQ → VCF (Sarek) | ✅ 已实现 | SarekRunner |
| Xena 公共数据 (TSV) | ✅ 已实现 | xena_adapter |
| 转录组 (RNA-seq) | ❌ 未实现 | — |
| 蛋白质组 | ❌ 未实现 | — |
| 代谢组 | ❌ 未实现 | — |
| 单细胞 scRNA-seq | ❌ 未实现 | — |

### 1.3 当前外部工具集成情况

| 工具 | 集成方式 | 可靠性 |
|------|---------|--------|
| **OpenMM** | Conda 安装，Stage 4 MD + 变异 minimization | ✅ 稳定 |
| **REINVENT4** | 可选 Stage 3 bridge | ⚠️ 可选，需单独安装 |
| **ESMFold** | HTTP API (esmatlas.com) | ⚠️ 长序列 413 错误 |
| **AlphaFold DB** | EBI 静态 PDB 下载 | ⚠️ 404 当无模型 |
| **ColabFold** | 未集成（仅文档/roadmap 提及） | ❌ |
| **AlphaFold 2/3** | 未集成（无本地推理） | ❌ |

---

## 2. DrugReflector 深度评估

### 2.1 它做了什么

**DrugReflector** (Cellarity, Science 2025) 是一个基于基因表达谱的化合物排序模型：

- **输入**: 基因表达差异签名（来自 single-cell atlas 的细胞状态转换定义）
- **模型**: 3 个 MLP 分类器的 ensemble，训练于 CMap 9,597 扰动 × 52 细胞系
- **输出**: 化合物类别预测 + 排名列表（~100 个优先化合物）
- **关键优势**:
  - 表型驱动（phenotype-first），不依赖已知靶点
  - 命中率比随机筛选高 **13-17 倍**
  - 支持主动学习 (lab-in-the-loop)：实验反馈 → 模型迭代 → 命中率再翻倍

**核心工作流**:
```
单细胞图谱 → 目标细胞状态转换签名 (gene signature)
    → DrugReflector 排名 → ~100 候选化合物
        → 表型验证 → 转录组反馈 → 模型迭代
```

### 2.2 它没做什么 — 为什么仍需 H2L/LO

DrugReflector 的输出是一个 **已知化合物的排名列表**，而非新分子的设计。以下能力它 **不具备**：

| 能力 | DrugReflector | T2Lead |
|------|:---:|:---:|
| 从多组学识别靶点 | 间接（表型排名，不输出靶点） | ✅ OpenTargets + 变异分析 |
| 新分子生成 (de novo) | ❌ | ✅ CReM + REINVENT4 |
| 类似物扩展 | ❌ | ✅ CReM analog_gen |
| Hit-to-Lead 优化 | ❌ | ✅ MPO, scaffold analysis |
| Lead Optimization | ❌ | ✅ Vina + ADMET + OpenMM MD |
| 蛋白-配体结构验证 | ❌ | ✅ 对接 + MD 模拟 |
| ADMET 预测 | ❌ | ✅ DeepADMET |
| 患者级别汇聚推荐 | ❌ | ✅ patient_aggregation |

**结论**:
> DrugReflector 擅长的是"从大量已知化合物中快速找到可能有效的苗头"，但它 **不替代** H2L 和 LO 过程。它找到的命中化合物（Hit）仍然需要：
> 1. 类似物扩展和结构优化（H2L）  
> 2. 对接验证、ADMET 评估、MD 稳定性验证（LO）  
> 3. 这正是 T2Lead Stage 3-4 做的事情

### 2.3 与 T2Lead 互补性分析

```
DrugReflector 的价值定位:

传统 T2Lead 路径:
  靶点发现 → ChEMBL 虚拟筛选 → Hit → H2L → LO
  问题: 依赖已知靶点, 虚拟筛选覆盖有限

加入 DrugReflector 后:
  多组学 → 基因签名 → DrugReflector 表型驱动排名 → Hit (更高命中率)
                                                      ↓
                                          → T2Lead H2L → LO (结构验证)
```

**互补优势**:
1. **无靶点也能工作**: 当变异分析未能找到可用靶点（druggable target）时，DrugReflector 可以绕过靶点，直接用表型签名搜索化合物
2. **更高命中率**: 13-17x 优于随机筛选，可以显著提升 Stage 2 的效率
3. **主动学习回路**: 实验数据反馈让模型越用越准
4. **与靶点路径并行**: 不冲突，可以同时跑两条路

**局限性**:
1. 仅限已知化合物库（CMap ~20K 化合物），无法发现全新分子
2. 训练数据偏向造血和肿瘤细胞系，其他组织可能表现不佳
3. 需要高质量的基因表达差异签名作为输入
4. 不提供分子-蛋白相互作用的结构解释

### 2.4 整合方案

**推荐: 将 DrugReflector 作为 Stage 2 的并行通路**

```python
# 概念性伪代码
class PipelineV2:
    def run_hit_discovery(self, config):
        hits_traditional = []
        hits_drugreflector = []

        # 通路 A: 传统靶点驱动 (现有 Stage 1-2)
        if targets := self.run_target_discovery():
            hits_traditional = self.run_chembl_screen(targets)

        # 通路 B: 表型驱动 (DrugReflector)
        if gene_signature := self.extract_gene_signature(omics_data):
            hits_drugreflector = self.run_drugreflector(gene_signature)

        # 合并去重
        all_hits = self.merge_and_deduplicate(
            hits_traditional, hits_drugreflector
        )

        # 继续 H2L → LO
        leads = self.run_hit_to_lead(all_hits)
        optimized = self.run_lead_optimization(leads)
        return optimized
```

**实现步骤**:
1. 新增 `src/drugpipe/target_to_hit/drugreflector_bridge.py`
2. 需要 RNA-seq / scRNA-seq 数据预处理模块（提取基因签名）
3. 在 `pipeline.py` 的 `run_target_to_hit` 中增加分支
4. 合并排名策略（加权融合或投票）

---

## 3. 多组学输入升级路径

### 3.1 现状 → 目标

```
现状:                          目标:
WES/WGS (VCF)           →    WES/WGS (VCF)
FASTQ (Sarek)            →    FASTQ (Sarek)
Xena TSV                 →    Xena TSV
                              RNA-seq (bulk + single-cell)    ← 新增
                              蛋白质组 (proteomics)           ← 新增
                              代谢组 (metabolomics)           ← 新增
                              临床表型数据                     ← 新增
```

### 3.2 各组学数据的用途

| 组学类型 | 管线中的作用 | 对应工具/模块 |
|---------|------------|-------------|
| **WES/WGS** | 体细胞突变 → 新抗原 → 靶点/药物靶标 | 现有 variant_analysis |
| **RNA-seq (bulk)** | 差异表达基因 → DrugReflector 输入; 融合基因检测 | DESeq2/edgeR → gene signature |
| **scRNA-seq** | 细胞状态图谱 → 精确表型签名; 肿瘤微环境分析 | Scanpy/Seurat → DrugReflector |
| **蛋白质组** | 验证转录本表达; 翻译后修饰 → 药靶 | OmniNeo 框架参考 |
| **代谢组** | 代谢通路异常 → 代谢酶靶点; 药物代谢预测 | 通路富集 → 靶点补充 |
| **临床数据** | 分期/分型 → 治疗方案分流; 预后评估 | 决策树输入 |

### 3.3 数据预处理模块设计

```
src/drugpipe/omics_integration/   ← 新增顶层模块
    __init__.py
    bulk_rnaseq.py          # DESeq2 差异表达 → gene signature
    scrna_processor.py      # Scanpy pipeline → cell state signatures
    proteomics_loader.py    # MaxQuant/DIA-NN 输出 → 蛋白表达矩阵
    metabolomics_loader.py  # 代谢组数据标准化
    omics_integrator.py     # 多组学融合 (MOFA+ / mixOmics 参考)
    signature_extractor.py  # 统一输出: gene_signature for DrugReflector
```

---

## 4. 外部工具链定位与优化

### 4.1 结构预测：ESMFold vs ColabFold vs AF2 vs AF3

| 工具 | 速度 | 准确度 | 特点 | 推荐场景 |
|------|------|--------|------|---------|
| **ESMFold** | ⚡ 最快 (~秒级) | 中等 (pLDDT ~70-85) | 单序列，无 MSA | 快速初筛; 短序列 (<400aa) |
| **ColabFold/AF2** | 中等 (~分钟级) | 高 (pLDDT ~85-95) | MSA-based, 经典验证 | 单体/多聚体标准预测 |
| **AlphaFold 3** | 慢 (~分钟-小时) | 最高 | 蛋白-配体/核酸/离子 | 药物-靶标复合物; RNA 结构 |
| **AF-DB 下载** | ⚡ 即时 | 高 (预计算) | 仅已知蛋白 | 野生型已知蛋白 |

**推荐分级策略 (Tiered Structure Resolution)**:

```
需要蛋白结构
    │
    ├─ PDB 实验结构存在? → 使用实验结构 (最高置信)
    │
    ├─ AF-DB 有预测结构? → 下载 AF-DB (快速, 高质量)
    │
    ├─ 是突变体/需要配体复合物?
    │   ├─ 是 → ColabFold/AF2 预测 (MSA 保证准确度)
    │   │       → AF3 预测蛋白-配体复合物 (如可用)
    │   └─ 否 → ESMFold 快速预测 (单序列足够)
    │
    └─ 序列 > 1000aa?
        ├─ 是 → ColabFold (分块策略) 或 AF3
        └─ 否 → ESMFold (快速) → ColabFold (验证)
```

**整合建议**:
1. **本地部署 ColabFold**: 解决 ESMFold API 的 413 错误和 AF-DB 404 问题
2. **AF3 作为高精度后备**: 用于药物-蛋白复合物建模（需申请权重）
3. **ESMFold 降级为预筛工具**: 快速过滤，不做最终决策依据
4. 在 `StructureBridge` 中实现这个分级策略

### 4.2 分子生成与优化：REINVENT4

**现有集成**: Stage 3 可选 bridge

**优化建议**:
1. **从可选提升为默认**: REINVENT4 在 de novo 设计和 scaffold hopping 上显著优于 CReM
2. **多目标优化**: 利用 REINVENT4 的 RL 框架，将 Vina docking score + ADMET + 合成可及性作为奖励函数
3. **与 DrugReflector 联动**: DrugReflector 排名 → top hits 的骨架 → REINVENT4 scaffold hopping → 新类似物
4. **linker design**: 对 PPI (protein-protein interaction) 靶点，用 REINVENT4 的 linker 模式

### 4.3 分子动力学：OpenMM

**现有集成**: Stage 4 隐式溶剂 MD + MM-GBSA 评分

**优化建议**:
1. **显式溶剂 MD**: 对 top candidates 运行显式溶剂模拟 (已有 config 但标记可选)，提高评分准确度
2. **FEP (自由能微扰)**: 对最终候选分子之间做相对结合自由能比较
3. **增强采样**: 对 hard-to-dock 靶点使用 metadynamics / replica exchange
4. **GPU 批量化**: 当前已支持多 GPU 调度，可进一步优化并行效率

### 4.4 推荐工具链组合

```
                    ┌─────────────────────────────────────────┐
                    │         完整工具链 (推荐)                │
                    ├─────────────────────────────────────────┤
  输入层             │  多组学 → omics_integrator              │
                    ├─────────────────────────────────────────┤
  靶点/表型层        │  OpenTargets (靶点) ∥ DrugReflector     │
                    │  (表型驱动)                              │
                    ├─────────────────────────────────────────┤
  Hit 发现          │  ChEMBL + DrugReflector 合并             │
                    ├─────────────────────────────────────────┤
  H2L               │  CReM + REINVENT4 (骨架跳跃/de novo)    │
                    ├─────────────────────────────────────────┤
  结构预测           │  PDB > AF-DB > ColabFold/AF2 > ESMFold  │
                    │  AF3 (复合物)                            │
                    ├─────────────────────────────────────────┤
  LO                │  Vina → DeepADMET → OpenMM MD/FEP       │
                    │  → 合成可及性 → 排名                     │
                    └─────────────────────────────────────────┘
```

---

## 5. 三轨治疗决策管线设计

### 5.1 Track A：小分子配体（靶向药物）

**适用条件**: 存在可成药靶点 (druggable target)

```
多组学分析 → 识别驱动突变/过表达靶点
    → 靶点可成药性评估 (DGIdb, OpenTargets tractability)
        → 可成药 → T2Lead 完整管线 (Stage 1-4)
        → 不可成药 → 转 Track B / Track C
```

**现有管线完全覆盖，需要的改进**:
- 增加可成药性 (druggability) 评估模块
- 增加 DrugReflector 表型路径作为并行 Hit 源
- 增强结构预测工具链

### 5.2 Track B：RNA 疫苗 / mRNA 治疗

**适用条件**:
1. 存在新抗原 (neoantigen)，但没有可成药的蛋白靶点
2. 免疫系统可被激活对抗肿瘤

**两种子方向**:

#### B1: 新抗原 mRNA 疫苗 (个性化肿瘤疫苗)
```
WES/WGS → 体细胞突变
    → 新抗原预测 (HLA 结合亲和力 + 免疫原性)
        → mRNA 序列设计
            → 密码子优化 + UTR 设计
                → 候选疫苗序列
```

**需要新增的模块**:
```
src/drugpipe/rna_vaccine/            ← 新增
    __init__.py
    neoantigen_predictor.py   # pVACseq / NetMHCpan / OmniNeo
    immunogenicity_scorer.py  # T细胞识别潜力评估
    mrna_designer.py          # 密码子优化, UTR 选择, polyA
    structure_predictor.py    # RNA 二级结构 (ViennaRNA/LinearFold)
    delivery_optimizer.py     # LNP 配方优化 (可选)
```

**关键外部工具**:
- **pVACseq / pVACtools**: 新抗原预测金标准
- **NetMHCpan 4.1**: HLA-I/II 结合预测
- **OmniNeo**: 多组学新抗原优化（含蛋白质组验证）
- **LinearFold**: 快速 RNA 二级结构预测
- **AF3**: RNA 三维结构预测（支持 RNA 输入）

#### B2: mRNA 治疗 (基因替代/编辑)
```
缺陷基因识别 → 治疗性蛋白设计
    → mRNA 序列优化
        → 蛋白折叠验证 (AF2/AF3/ESMFold)
```

### 5.3 Track C：放射治疗（质子/重离子/靶向放射性核素）

**适用条件**:
1. 没有可成药靶点 且 没有合适的新抗原
2. 肿瘤位置适合精准放疗
3. 与免疫治疗联合增效

**三种子方向**:

| 方式 | 原理 | 适用场景 | 计算需求 |
|------|------|---------|---------|
| **质子治疗 (Proton)** | Bragg 峰精准沉积 | 位置敏感肿瘤（脑、儿童） | 蒙特卡洛剂量规划 |
| **重离子治疗 (Carbon)** | 更高 RBE，复杂 DNA 损伤 | 放射抗性肿瘤 | 同上 |
| **靶向放射性核素 (TRT)** | 放射性同位素 + 靶向载体 | 有表面标记但不可小分子成药 | 药代动力学建模 |
| **FLASH 放疗** | 超高剂量率，保护正常组织 | 新兴技术，临床试验阶段 | 剂量-时间建模 |

**与免疫治疗联合** (参考 PMC11885950):
- 质子治疗可诱导免疫原性细胞死亡 (ICD)
- 与 checkpoint inhibitor 联合有协同效应
- 可与 Track B RNA 疫苗联合

**T2Lead 中的定位**:
> Track C 主要是临床决策层面的推荐，计算管线的角色是：
> 1. 识别适合 TRT 的表面靶标（如果有）
> 2. 输出"无法通过小分子/疫苗治疗"的判断 + 推荐放疗方案
> 3. 不需要自建放疗剂量规划系统

### 5.4 决策树与分流逻辑

```
                        患者多组学数据输入
                              │
                    ┌─────────┴──────────┐
                    │   多组学预处理       │
                    │   突变/表达/蛋白/代谢 │
                    └─────────┬──────────┘
                              │
                ┌─────────────┼──────────────┐
                ▼             ▼              ▼
        ┌──────────┐  ┌───────────┐  ┌──────────────┐
        │ 靶点识别  │  │ 新抗原预测 │  │ 表型签名提取  │
        │(变异+表达)│  │(WES+HLA)  │  │(RNA-seq/scRNA)│
        └────┬─────┘  └─────┬─────┘  └──────┬───────┘
             │              │               │
             ▼              ▼               ▼
      ┌──────────┐   ┌──────────┐    ┌────────────┐
      │可成药靶点?│   │强新抗原? │    │表型匹配化合物│
      │druggable │   │immunogenic│    │DrugReflector│
      └──┬───┬───┘   └──┬───┬───┘    └──────┬─────┘
         │   │          │   │               │
        Yes  No        Yes  No              │
         │   │          │   │               │
         ▼   │          ▼   │               ▼
    ┌────────┐│    ┌────────┐│        ┌──────────┐
    │Track A ││    │Track B ││        │ 合并 Hit │
    │小分子  ││    │RNA疫苗 ││        │→ H2L→LO  │
    │H2L→LO ││    │设计    ││        └──────────┘
    └────────┘│    └────────┘│
              │              │
              ▼              ▼
        ┌──────────────────────┐
        │ 都不可行?             │
        │ → Track C: 放疗推荐  │
        │ → 联合免疫治疗评估    │
        │ → 临床试验匹配       │
        └──────────────────────┘

注意: Track A / B / C 不互斥，可以并行/联合
```

**分流条件定义** (需讨论确认):

| 判断点 | 条件 | 分流方向 |
|--------|------|---------|
| 可成药靶点 | DGIdb hit + 结构可预测 + 结合口袋存在 | → Track A |
| 强新抗原 | NetMHCpan IC50 < 500nM + 表达量 > 阈值 + 克隆频率 > 0.1 | → Track B |
| DrugReflector 命中 | 排名 top-K 化合物 + 表型相关性 > 阈值 | → Track A (合并 hit) |
| 无靶点 + 无新抗原 | 以上均不满足 | → Track C |
| 联合方案 | Track A 候选 + Track B 候选同时存在 | → 联合推荐 |

---

## 6. 整体架构蓝图（新）

```
┌─────────────────────────────────────────────────────────────────────┐
│                    T2Lead v2 — 多模态精准治疗管线                     │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                    Layer 0: 数据输入层                        │   │
│  │  WES/WGS  RNA-seq  scRNA-seq  蛋白质组  代谢组  临床数据     │   │
│  │  (VCF)   (FASTQ)  (h5ad)    (MaxQuant) (mzML)  (表型/分期) │   │
│  └──────────────────────┬───────────────────────────────────────┘   │
│                         │                                           │
│  ┌──────────────────────▼───────────────────────────────────────┐   │
│  │                 Layer 1: 多组学整合层                         │   │
│  │  omics_integrator: 标准化 → 质控 → 整合 → 特征提取           │   │
│  │  输出: 突变列表, 基因签名, 表达矩阵, 通路活性, 新抗原候选     │   │
│  └──────────┬──────────────┬────────────────┬───────────────────┘   │
│             │              │                │                       │
│  ┌──────────▼──────┐ ┌────▼────────┐ ┌─────▼──────────────┐       │
│  │  Layer 2A:      │ │ Layer 2B:   │ │ Layer 2C:          │       │
│  │  靶点发现       │ │ 表型驱动    │ │ 新抗原预测         │       │
│  │  (OpenTargets   │ │ (Drug-      │ │ (pVACseq +         │       │
│  │   + 变异分析)   │ │  Reflector) │ │  NetMHCpan)        │       │
│  └────────┬────────┘ └──────┬──────┘ └────────┬────────────┘       │
│           │                 │                  │                     │
│  ┌────────▼─────────────────▼──────┐  ┌───────▼────────────┐       │
│  │  Layer 3: Hit Discovery + H2L   │  │ Layer 3B:          │       │
│  │  ChEMBL + DR hits → 合并        │  │ RNA 疫苗设计       │       │
│  │  → CReM + REINVENT4 类似物       │  │ mRNA 序列优化      │       │
│  │  → MPO 多目标优化                │  │ 免疫原性评分       │       │
│  └────────────────┬────────────────┘  └───────┬────────────┘       │
│                   │                           │                     │
│  ┌────────────────▼────────────────┐          │                     │
│  │  Layer 4: Lead Optimization     │          │                     │
│  │  结构预测 (PDB/AF-DB/ColabFold  │          │                     │
│  │    /AF2/AF3/ESMFold)            │          │                     │
│  │  → Vina 对接                    │          │                     │
│  │  → DeepADMET                    │          │                     │
│  │  → OpenMM MD (隐式→显式)        │          │                     │
│  │  → 合成可及性评估               │          │                     │
│  └────────────────┬────────────────┘          │                     │
│                   │                           │                     │
│  ┌────────────────▼───────────────────────────▼────────────────┐   │
│  │  Layer 5: 治疗推荐层                                        │   │
│  │  ┌─────────┐  ┌──────────┐  ┌──────────┐  ┌────────────┐  │   │
│  │  │Track A  │  │Track B   │  │Track C   │  │联合方案    │  │   │
│  │  │小分子药 │  │RNA疫苗   │  │放疗推荐  │  │A+B / B+C  │  │   │
│  │  │排名报告 │  │序列报告  │  │临床推荐  │  │综合方案   │  │   │
│  │  └─────────┘  └──────────┘  └──────────┘  └────────────┘  │   │
│  └─────────────────────────────────────────────────────────────┘   │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │  Layer 6: 患者报告 (patient_aggregation)                    │   │
│  │  汇聚所有 Track 结果 → 综合治疗推荐 → 可视化报告             │   │
│  └─────────────────────────────────────────────────────────────┘   │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│  基础设施: OpenMM | REINVENT4 | ColabFold | AF2/3 | ESMFold |     │
│            pVACseq | DeepADMET | Vina | CReM | Sarek              │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 7. 与现有代码的兼容性评估

### 7.1 可直接复用的部分

| 现有模块 | 在新架构中的角色 | 需要改动 |
|---------|----------------|---------|
| `variant_analysis/` | Layer 1 的突变分析组件 | 小改：增加输出接口 |
| `target_discovery/` | Layer 2A | 基本不变 |
| `target_to_hit/` | Layer 3 的 ChEMBL 通路 | 增加 DrugReflector 并行分支 |
| `hit_to_lead/` | Layer 3 | 增强 REINVENT4 默认启用 |
| `lead_optimization/` | Layer 4 | 增强结构预测分级策略 |
| `patient_aggregation` | Layer 6 | 扩展支持多 Track 合并 |
| `pipeline.py` | 顶层编排 | 重构为多轨调度 |
| `paths.py` | 路径管理 | 扩展 Track B/C 目录 |
| `config.py` / YAML | 配置 | 新增组学/DrugReflector/疫苗配置段 |

### 7.2 需要新增的部分

| 新模块 | 优先级 | 复杂度 | 估计工时 |
|-------|--------|--------|---------|
| `omics_integration/` | P0 | 高 | 2-3 周 |
| `drugreflector_bridge.py` | P0 | 中 | 1 周 |
| `rna_vaccine/` | P1 | 高 | 3-4 周 |
| ColabFold 本地集成 | P1 | 中 | 1 周 |
| AF3 集成 | P2 | 中 | 1-2 周 (需权重申请) |
| `treatment_router.py` (分流逻辑) | P1 | 中 | 1 周 |
| 放疗推荐模块 | P2 | 低 | 2-3 天 (主要是知识库) |

### 7.3 风险评估

| 风险 | 影响 | 缓解 |
|------|------|------|
| DrugReflector 训练数据偏向特定细胞系 | 非造血/肿瘤场景效果不确定 | 先在乳腺癌验证，逐步扩展 |
| AF3 权重获取不确定 | 复合物预测能力受限 | ColabFold+AF2 作为 fallback |
| 多组学数据质量参差 | 签名质量影响 DrugReflector 效果 | 严格 QC 流程 |
| RNA 疫苗模块复杂度高 | 开发周期长 | 分阶段：先新抗原预测，再 mRNA 设计 |
| 三轨合并报告复杂 | 临床可解释性挑战 | 参考 Knowledge Connector 决策框架 |

---

## 8. 待讨论问题清单

### 高优先级

- [ ] **Q1**: DrugReflector 在乳腺癌场景的适用性如何？其训练数据中造血细胞系较多，需要评估迁移性
- [ ] **Q2**: 多组学输入的最小可行集合是什么？是否先只做 WES + RNA-seq？
- [ ] **Q3**: RNA 疫苗路径是否需要完整的 mRNA 设计，还是只到新抗原预测为止？
- [ ] **Q4**: 三轨是否真正并行执行，还是有条件串行（先判断再走特定轨道）？
- [ ] **Q5**: 计算资源限制：ColabFold 本地部署需要多少 GPU 内存？AF3 需要 A100?

### 中优先级

- [ ] **Q6**: DrugReflector 的 CMap 化合物库 (~20K) 是否覆盖我们关注的治疗领域？
- [ ] **Q7**: 是否考虑 ADC (antibody-drug conjugate) 作为第四轨道？
- [ ] **Q8**: REINVENT4 是否应该完全替代 CReM，还是保持双路？
- [ ] **Q9**: 现有的 Xena 公共数据测试框架是否可以扩展来验证新管线？
- [ ] **Q10**: patient_aggregation 的权重体系是否需要重新设计？

### 低优先级

- [ ] **Q11**: 是否需要 GUI / Web 界面让临床人员直接使用？
- [ ] **Q12**: 数据隐私合规（患者多组学数据）
- [ ] **Q13**: 是否考虑联邦学习来保护患者数据？

---

## 9. 参考文献

1. **DrugReflector**: Cellarity. "Compound ranking predictions from gene expression signatures using ensemble neural network models." *Science* (2025). [GitHub](https://github.com/Cellarity/drugreflector)
2. **REINVENT4**: Loeffler et al. "Reinvent 4: Modern AI-driven generative molecule design." *J. Cheminformatics* (2024). [GitHub](https://github.com/MolecularAI/REINVENT4)
3. **OpenMM**: Eastman et al. "OpenMM: A toolkit for molecular simulation." [GitHub](https://github.com/openmm/openmm)
4. **ESM-2 / ESMFold**: Lin et al. "Evolutionary-scale prediction of atomic-level protein structure with a language model." *Science* 379(6637):1123-1130 (2023). [GitHub](https://github.com/facebookresearch/esm)
5. **AlphaFold 2**: Jumper et al. "Highly accurate protein structure prediction with AlphaFold." *Nature* 596:583-589 (2021). [GitHub](https://github.com/google-deepmind/alphafold)
6. **AlphaFold 3**: Abramson et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493-500 (2024). [GitHub](https://github.com/google-deepmind/alphafold3)
7. **ColabFold**: Mirdita et al. "ColabFold: Making protein folding accessible to all." *Nature Methods* (2022). [Protocol](https://github.com/steineggerlab/colabfold-protocol)
8. **Knowledge Connector**: "The Knowledge Connector decision support system for multiomics-based precision oncology." *Nature Communications* (2026).
9. **OmniNeo**: "OmniNeo: a multi-omics pipeline incorporating proteomics and AI selection for neoantigen optimization." *Frontiers in Immunology* (2025).
10. **质子治疗免疫效应**: "Immunological Effects of Proton Radiotherapy." PMC11885950 (2025).
11. **靶向放射性核素**: "The molecular blueprint of targeted radionuclide therapy." *Nature Reviews Clinical Oncology* (2025).
12. **AI 多组学精准肿瘤学**: "AI-driven multi-omics integration in precision oncology." *Clinical and Experimental Medicine* (2025).

---

## 10. 基础设施优化：立即可做的改进

> 以下是基于代码审计和日志分析的具体改进建议，不需要等架构重构就能做。

### 10.1 UniProt 本地化 — 消灭重复网络请求

**现状问题**: 同一 UniProt ID 在多个变异中被重复 GET 5 次 + 重试，纯浪费时间。

**UniProt 大小参考**:

| 数据库 | 压缩大小 | 解压大小 | 说明 |
|--------|---------|---------|------|
| **Swiss-Prot FASTA** | **89 MB** | ~250 MB | 已审核蛋白，~57 万条 |
| **Swiss-Prot DAT** | 660 MB | ~4 GB | 含注释、交叉引用等 |
| TrEMBL FASTA | 49 GB | ~120 GB | 全部未审核蛋白 |
| **人类子集 (idmapping)** | ~200 MB | ~500 MB | 只要人类的就够了 |

**推荐方案**: 下载 Swiss-Prot FASTA + 人类 idmapping 即可，总共不到 1 GB。

```bash
# 下载 Swiss-Prot FASTA (89MB 压缩)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# 下载人类蛋白子集 ID 映射 (Gene → UniProt → PDB)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
```

**代码改进** — 增加本地 UniProt 查找缓存:

```python
# src/drugpipe/utils/uniprot_cache.py (新增)
class UniProtLocalCache:
    """本地 Swiss-Prot FASTA + ID mapping 缓存，消灭重复网络请求"""

    def __init__(self, db_dir: Path):
        self._fasta_index = {}   # uniprot_id -> sequence
        self._gene_map = {}      # gene_symbol -> [uniprot_ids]
        self._load(db_dir)

    def get_sequence(self, uniprot_id: str) -> Optional[str]:
        return self._fasta_index.get(uniprot_id)

    def gene_to_uniprot(self, gene: str) -> List[str]:
        return self._gene_map.get(gene, [])
```

### 10.2 ESMFold 本地部署 — 告别远程 413

**现状问题**: 远程 ESMFold API 对长序列 (>400aa) 频繁返回 413，重试 5 次每次都失败，浪费 ~5 分钟/序列。

**你的硬件**: RTX 4090 24GB VRAM

| 序列长度 | 预计 VRAM | 4090 能跑? |
|---------|----------|-----------|
| ≤300 aa | ~16 GB | ✅ |
| ≤500 aa | ~20 GB | ✅ |
| ≤700 aa | ~24 GB | ⚠️ 边界 (开 chunk) |
| >700 aa | >24 GB | ❌ 需 ColabFold 接管 |

**部署方案** (推荐 LiteFold Docker):

```bash
# 方案 A: LiteFold (Docker + FastAPI，开箱即用)
git clone https://github.com/LiteFold/litefold.git
cd litefold && docker compose up -d
# 本地 API: http://localhost:8000/predict

# 方案 B: 直接 Python 调用 (更灵活)
pip install "fair-esm[esmfold]"
pip install "openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059"
```

**代码改进** — `structure_bridge.py` 分级策略:

```python
# 改动: _predict_esmfold 方法
_ESMFOLD_LOCAL_URL = "http://localhost:8000/predict"  # 本地部署
_ESMFOLD_MAX_LOCAL = 700   # 本地 4090 最大序列长度
_ESMFOLD_SKIP_THRESHOLD = 1200  # 超过此长度直接跳过 ESMFold

def _predict_esmfold(self, mp, out_dir):
    seq_len = len(mp.mutant_sequence)

    if seq_len > _ESMFOLD_SKIP_THRESHOLD:
        logger.info("%s: length %d > %d, skip ESMFold → next fallback",
                     mp.gene_symbol, seq_len, _ESMFOLD_SKIP_THRESHOLD)
        return None  # 直接进 ColabFold/AF-DB fallback

    if seq_len <= _ESMFOLD_MAX_LOCAL:
        result = self._call_local_esmfold(mp, out_dir)  # 本地推理
        if result:
            return result

    # 超过本地限制 → ColabFold (而不是远程 ESMFold)
    return self._predict_colabfold(mp, out_dir)
```

### 10.3 全局结构缓存 — 消灭重复计算

**现状**: `_pdb_cache` 只缓存 RCSB 搜索结果，不跨变异共享。同一基因的 5 个突变各自下载 + 预测一遍。

**改进**: 增加文件级全局缓存，按 UniProt+突变 去重:

```python
class StructureBridge:
    # 增加类级别的全局缓存
    _global_structure_cache: Dict[str, Path] = {}

    def _resolve_one(self, mp, out_dir):
        cache_key = f"{mp.uniprot_id}_{mp.mutation_label}"
        if cache_key in self._global_structure_cache:
            cached = self._global_structure_cache[cache_key]
            if cached.exists():
                logger.info("Cache hit: %s → %s", cache_key, cached)
                return {..., "method": "cached", "pdb_path": str(cached)}

        result = self._resolve_one_uncached(mp, out_dir)
        if result:
            self._global_structure_cache[cache_key] = Path(result["pdb_path"])
        return result
```

### 10.4 长序列快速失败策略

**现状**: 日志证明 >400aa 序列远程 ESMFold 几乎必定 413，但代码仍然尝试 5 次重试 × 120s timeout = 浪费 10 分钟。

**立即改动** (不需要本地部署也能做):

```python
# structure_bridge.py 的 _predict_esmfold:
def _predict_esmfold(self, mp, out_dir):
    seq = mp.mutant_sequence
    # 超过阈值直接跳过，不浪费重试时间
    if len(seq) > self._esmfold_max_len:  # 配置项，默认 600
        logger.info(
            "%s %s: seq length %d > %d, skipping remote ESMFold (known 413)",
            mp.gene_symbol, mp.mutation_label, len(seq), self._esmfold_max_len,
        )
        return None
    # ... 正常流程
```

---

## 11. 三人分工方案

> 前提：你（总负责人）负责任务下发、进度追踪、最终集成。三个研究所各自独立开发，**保证输出能成为别人的输入**（接口标准化）。

### 11.1 总体分工逻辑

基于管线数据流的自然切分点：

```
患者多组学数据
    │
    ├──→ [研究所 A] DNA/变异 → 靶点识别 → 可成药性判断
    │         输出: 靶点列表 + 突变蛋白序列 + 新抗原候选
    │              ↓                    ↓
    │    [研究所 C]                [研究所 B]
    │    靶点→小分子配体           新抗原→mRNA疫苗
    │    (含DrugReflector          设计+验证
    │     + H2L + LO)
    │         输出: 候选药物排名    输出: 候选疫苗序列
    │              ↓                    ↓
    └──→ [你] 集成层: 合并 Track A/B/C → 综合治疗推荐
```

---

### 研究所 A — "上游：多组学分析 + 靶点识别 + 新抗原预测"

**一句话**: 从患者原始数据到「靶点列表 + 突变序列 + 新抗原候选」

#### 负责范围

| 模块 | 说明 | 现有基础 |
|------|------|---------|
| 多组学数据预处理 | WES/WGS → VCF; RNA-seq → 差异表达; 蛋白质组/代谢组标准化 | T2Lead 已有 `variant_analysis/` VCF 部分 |
| 靶点识别 | 驱动突变识别 + OpenTargets 排名 + **可成药性评估** | T2Lead 已有 `target_discovery/` |
| 新抗原预测 | pVACseq / NetMHCpan → HLA 结合 → 免疫原性评分 | **新建** |
| 基因表达签名提取 | RNA-seq / scRNA-seq → gene signature (供研究所 C 的 DrugReflector 使用) | **新建** |

#### 输出接口定义 (供下游消费)

```yaml
# target_report.yaml — 研究所 A 的标准化输出
patient_id: "PT-001"
cancer_type: "breast_invasive_ductal"

targets:
  - gene: "PIK3CA"
    uniprot_id: "P42336"
    mutation: "H1047R"
    druggable: true
    druggability_score: 0.92
    wildtype_sequence: "MPPRP..."
    mutant_sequence: "MPPRP...R..."
    evidence: ["hotspot", "oncogene", "FDA-approved-target"]

neoantigens:
  - peptide: "KQFIDNTML"
    hla_allele: "HLA-A*02:01"
    binding_affinity_nM: 45.3
    immunogenicity_score: 0.87
    source_mutation: "PIK3CA_H1047R"
    expression_tpm: 23.5

gene_signatures:
  - name: "tumor_vs_normal"
    up_genes: ["ERBB2", "MYC", ...]     # top 100 上调
    down_genes: ["TP53", "RB1", ...]     # top 100 下调
    signature_type: "differential_expression"

treatment_routing:
  has_druggable_targets: true
  has_strong_neoantigens: true
  recommended_tracks: ["A", "B"]  # A=小分子, B=mRNA疫苗
```

#### 关键工具栈

- Sarek (FASTQ→VCF，T2Lead 已有)
- DESeq2 / edgeR (差异表达)
- Scanpy (scRNA-seq)
- pVACseq + NetMHCpan 4.1 (新抗原)
- OpenTargets API (靶点排名)
- DGIdb (可成药性)

#### 交付里程碑

| 阶段 | 内容 | 时间 |
|------|------|------|
| M1 | WES→VCF→靶点 + 可成药性判断（复用 T2Lead 现有代码） | 2 周 |
| M2 | RNA-seq 差异表达 → gene signature 输出 | 2 周 |
| M3 | 新抗原预测流程 (pVACseq + NetMHCpan) | 3 周 |
| M4 | 多组学整合 + 标准化输出接口 | 2 周 |

---

### 研究所 B — "RNA/mRNA 疫苗设计"

**一句话**: 从新抗原候选到「候选 mRNA 疫苗序列 + 免疫原性评分」

#### 负责范围

| 模块 | 说明 | 现有基础 |
|------|------|---------|
| 新抗原筛选与排序 | 接收研究所 A 的 neoantigen 列表 → 多维评分排序 | **新建** (消费 A 的输出) |
| mRNA 序列设计 | 密码子优化 + 5'UTR/3'UTR 选择 + polyA 尾 | **新建** |
| mRNA 结构预测 | RNA 二级结构 (ViennaRNA/LinearFold) + 稳定性评估 | **新建** |
| 免疫原性验证 | T细胞识别潜力建模 + 人群覆盖率 (HLA 频率) | **新建** |
| [可选] LNP 递送优化 | 脂质纳米颗粒配方参数 | 后续阶段 |

#### 输入接口 (消费研究所 A 的输出)

```yaml
# 研究所 B 直接读取研究所 A 的 target_report.yaml 中的 neoantigens 段
input_from_A:
  neoantigens:
    - peptide: "KQFIDNTML"
      hla_allele: "HLA-A*02:01"
      binding_affinity_nM: 45.3
      source_mutation: "PIK3CA_H1047R"
```

#### 输出接口定义

```yaml
# vaccine_report.yaml — 研究所 B 的标准化输出
patient_id: "PT-001"

vaccine_candidates:
  - rank: 1
    target_neoantigens: ["KQFIDNTML", "YLAEIFKAL"]
    mrna_sequence: "AUGGCUAAG..."        # 完整 mRNA 序列
    codon_optimization_index: 0.91
    mrna_length_nt: 1247
    predicted_stability_half_life_h: 8.2
    secondary_structure_mfe: -342.5       # kcal/mol
    population_coverage: 0.78             # HLA 人群覆盖率
    immunogenicity_composite_score: 0.85

  - rank: 2
    ...

design_parameters:
  utr_5: "Kozak_optimized_v2"
  utr_3: "beta_globin"
  poly_a_length: 120
  cap_structure: "CleanCap_AG"
```

#### 关键工具栈

- LinearFold / ViennaRNA (RNA 结构)
- CodonOptimizer (密码子优化)
- NetMHCpan 4.1 (HLA 结合亲和力，与 A 共用)
- **AlphaFold 3** (RNA 三维结构预测 — AF3 支持 RNA 输入)
- IEDB 工具集 (免疫原性)
- MHCflurry (机器学习 HLA 预测，补充)

#### 交付里程碑

| 阶段 | 内容 | 时间 |
|------|------|------|
| M1 | 新抗原排序流程 + 输入接口对接 | 2 周 |
| M2 | mRNA 序列设计 (密码子优化 + UTR) | 3 周 |
| M3 | RNA 结构预测 + 稳定性评估 | 2 周 |
| M4 | 人群覆盖率 + 综合评分 + 输出接口 | 2 周 |

---

### 研究所 C — "小分子配体：靶点到候选药物"

**一句话**: 从靶点/基因签名到「候选小分子药物排名」

#### 负责范围

| 模块 | 说明 | 现有基础 |
|------|------|---------|
| DrugReflector 表型驱动筛选 | 接收 gene signature → 化合物排名 | **新建** (整合 [Cellarity/drugreflector](https://github.com/Cellarity/drugreflector)) |
| 虚拟筛选 (ChEMBL) | 靶点 → ChEMBL 已知活性化合物 | T2Lead 已有 `target_to_hit/` |
| Hit-to-Lead | CReM 类似物 + REINVENT4 scaffold hopping + MPO | T2Lead 已有 `hit_to_lead/` |
| Lead Optimization | 结构预测 + Vina 对接 + ADMET + OpenMM MD | T2Lead 已有 `lead_optimization/` |
| 结构预测分级 | PDB → AF-DB → **本地 ESMFold** → ColabFold → AF3 | T2Lead 部分有，需**升级** |

#### 输入接口 (消费研究所 A 的输出)

```yaml
# 研究所 C 读取研究所 A 的 target_report.yaml 中的 targets + gene_signatures
input_from_A:
  targets:
    - gene: "PIK3CA"
      uniprot_id: "P42336"
      mutation: "H1047R"
      druggable: true
      mutant_sequence: "..."
  gene_signatures:
    - name: "tumor_vs_normal"
      up_genes: [...]
      down_genes: [...]
```

#### 输出接口定义

```yaml
# drug_report.yaml — 研究所 C 的标准化输出
patient_id: "PT-001"

drug_candidates:
  - rank: 1
    smiles: "CC(=O)Oc1ccccc1C(=O)O"
    name: "Alpelisib"
    source: "drugreflector+chembl"        # 来源标记
    target_gene: "PIK3CA"
    vina_score_kcal: -9.2
    admet:
      oral_bioavailability: "high"
      herg_risk: "low"
      hepatotoxicity: "low"
    md_binding_energy_kcal: -45.3
    synthetic_accessibility: 3.2
    composite_score: 0.91
    structure_method: "experimental_mutated_localopt(3HHM)"

  - rank: 2
    smiles: "..."
    source: "reinvent4_denovo"            # REINVENT4 生成
    ...

pipeline_metadata:
  structure_prediction_tool: "ColabFold_AF2"
  docking_tool: "AutoDock_Vina_1.2"
  md_tool: "OpenMM_8.1_implicit"
  scoring: "Vina + ADMET + MD composite"
```

#### 关键工具栈

- **DrugReflector** (表型驱动化合物排名)
- T2Lead `target_to_hit/` + `hit_to_lead/` + `lead_optimization/` (核心管线)
- **REINVENT4** (de novo 分子生成)
- **OpenMM** (MD 模拟)
- **ColabFold** (本地结构预测 — 重点部署)
- **ESMFold** (本地部署，快速预筛)
- AutoDock Vina (对接)
- DeepADMET (ADMET 预测)

#### 研究所 C 额外基础设施任务

由于 C 承接了最重的计算负载，还需要做：

1. **UniProt 本地缓存** → 消灭重复网络请求 (见 §10.1)
2. **ESMFold 本地部署** → 告别 413 (见 §10.2)
3. **ColabFold 本地部署** → 长序列 + 高精度结构
4. **全局结构缓存** → 跨变异去重 (见 §10.3)

#### 交付里程碑

| 阶段 | 内容 | 时间 |
|------|------|------|
| M1 | 基础设施: UniProt 本地 + ESMFold 本地 + 结构缓存 + 长序列快速失败 | 1.5 周 |
| M2 | DrugReflector 整合 (bridge + gene signature → 化合物排名) | 2 周 |
| M3 | REINVENT4 从可选提升为默认 + 多目标奖励函数 | 2 周 |
| M4 | 结构预测分级策略 (PDB → AF-DB → ESMFold → ColabFold) | 2 周 |
| M5 | 完整管线联调 + 输入/输出接口标准化 | 1.5 周 |

---

### 11.2 你（总负责人）的职责

| 职责 | 具体内容 |
|------|---------|
| **接口标准化** | 定义并维护 `target_report.yaml` / `vaccine_report.yaml` / `drug_report.yaml` 的 schema |
| **集成层** | 开发 `treatment_router.py`（分流逻辑）+ `patient_report_generator.py`（合并三轨） |
| **进度追踪** | 每周同步，按里程碑验收 |
| **最终管线编排** | 改造 `pipeline.py` 支持三轨调度 |
| **Track C 兜底** | 放疗推荐模块（知识库 + 规则，工作量小） |
| **测试数据** | 提供乳腺癌测试集 (Xena 已有) 供三方联调 |

### 11.3 协作接口图

```
┌─────────────────────────────────────────────────┐
│                   你（总负责）                    │
│  接口定义 / 集成层 / 进度 / Track C 兜底         │
└───────┬──────────────┬──────────────┬────────────┘
        │              │              │
        ▼              ▼              ▼
┌──────────────┐ ┌──────────────┐ ┌──────────────┐
│ 研究所 A      │ │ 研究所 B      │ │ 研究所 C      │
│ 多组学→靶点   │ │ 新抗原→mRNA  │ │ 靶点→小分子   │
│ +新抗原预测   │ │ 疫苗设计     │ │ +DrugReflector│
├──────────────┤ ├──────────────┤ ├──────────────┤
│输入:          │ │输入:          │ │输入:          │
│ 患者FASTQ/VCF│ │ A的neoantigen│ │ A的targets   │
│ RNA-seq      │ │              │ │ A的gene_sig  │
│ 蛋白/代谢组  │ │              │ │              │
├──────────────┤ ├──────────────┤ ├──────────────┤
│输出:          │ │输出:          │ │输出:          │
│target_report │ │vaccine_report│ │drug_report   │
│ .yaml        │ │ .yaml        │ │ .yaml        │
└──────┬───────┘ └──────────────┘ └──────────────┘
       │                ▲                ▲
       ├────────────────┘                │
       └─────────────────────────────────┘
       (A 的输出 feed B 和 C)
```

### 11.4 T2Lead 代码拆分建议

**不建议拆仓库**，而是按模块权限分配：

| 模块目录 | 主要负责 | 其他人 |
|---------|---------|--------|
| `src/drugpipe/omics_integration/` (新) | 研究所 A | — |
| `src/drugpipe/target_discovery/` | 研究所 A | C 只读 |
| `src/drugpipe/variant_analysis/` | 研究所 A | C 只读 |
| `src/drugpipe/rna_vaccine/` (新) | 研究所 B | — |
| `src/drugpipe/target_to_hit/` | 研究所 C | — |
| `src/drugpipe/hit_to_lead/` | 研究所 C | — |
| `src/drugpipe/lead_optimization/` | 研究所 C | — |
| `src/drugpipe/drugreflector_bridge.py` (新) | 研究所 C | — |
| `src/drugpipe/pipeline.py` | 你 | 三方 PR |
| `src/drugpipe/treatment_router.py` (新) | 你 | — |
| `src/drugpipe/patient_aggregation/` | 你 | 三方扩展 |
| `configs/` | 你 | 三方提交子配置 |

每人在自己负责的目录内自由开发，通过标准化 YAML 接口串联。Git 分支策略: 各研究所在 `feat/institute-X-xxx` 分支开发，PR 合并到 `dev`。

---

## 12. Q&A 第二轮（2026-04-15 讨论记录）

### Q: ColabFold 本地数据库真的 ~1TB 吗？

**是的，大致准确**。ColabFold 完全本地化需要 `colabfold_search` + MMseqs2 数据库：

| 数据库 | 压缩 | 解压+索引 |
|--------|------|----------|
| UniRef30 (2023-02) | 96 GB | ~182 GB |
| ColabFoldDB (envdb) | 110 GB | ~597 GB |
| **合计** | **~206 GB** | **~780 GB** |

实际占用约 **800 GB**（不到 1TB，但确实很大）。

**替代方案**: ColabFold 支持 MMseqs2 API 模式 — 不下载数据库，MSA 搜索发到 Söding Lab 免费服务器，只本地跑结构预测。这样只需 ~5 GB（模型权重）。但有网络依赖和排队延迟。

**建议**: 先用 API 模式跑起来（0 存储成本），如果吞吐量成瓶颈再考虑下载完整数据库。

### Q: AF3 要不要本地部署？谁来做？

**AF3 本地部署对你的硬件不现实**：

| 需求 | AF3 要求 | 你的机器 | 结论 |
|------|---------|---------|------|
| GPU | A100 80GB 或 H100 80GB | RTX 4090 24GB | ❌ 差太多 |
| 磁盘 | ~630 GB (数据库) + ~1 TB (含索引) | — | 可以满足 |
| RAM | 64 GB+ | 看配置 | 可能不够 |
| 权重 | 需向 Google 申请 (2-3 工作日) | — | 需单独申请 |

**建议**:
1. **先不本地部署 AF3** — 用 [AlphaFold Server](https://alphafoldserver.com) 的免费在线服务验证少量关键复合物
2. AF3 本地部署留给有 A100 的研究所（如果某个研究所有高端 GPU 集群）
3. AF3 的独特价值是蛋白-配体/RNA 复合物预测 + 甲基化/糖基化 — 只对最终 top candidates 用，不需要批量
4. **谁做**: 如果三个研究所中有一个有 A100 集群，归那个人；否则暂用在线服务

### Q: DrugReflector 可以拆开吗？多组学→签名 和 签名→化合物排名 能分给不同人吗？

**可以拆**。DrugReflector 的工作流分两步，天然可以在中间切开：

```
步骤 1 (归研究所 A):
  多组学数据 → 差异表达分析 → gene signature (上调/下调基因列表)
  这是标准的生信分析，不依赖 DrugReflector

步骤 2 (归研究所 C):
  gene signature → DrugReflector MLP 推理 → 化合物排名
  这是 DrugReflector 的核心部分
```

**接口切点**: `gene_signature` (一个基因列表 + fold-change 向量) — 研究所 A 产生，研究所 C 消费。
这正好对应 `target_report.yaml` 中的 `gene_signatures` 段。DrugReflector 本身不做差异表达，它只接受签名输入。

### Q: mRNA 疫苗方向，人下周才到，我先怎么探路？

**你本周可以做的探路工作**（不需要研究所 B 的人）：

#### 第 1 天：环境 + 工具调研
```bash
# 1. 安装 pVACtools (新抗原预测核心工具)
conda create -n pvac python=3.10
conda activate pvac
pip install pvactools

# 2. 快速验证
pvacseq run --help

# 3. 安装 ViennaRNA (mRNA 二级结构)
conda install -c bioconda viennarna
```

#### 第 2-3 天：用已有 Xena 乳腺癌数据跑一个 demo
```bash
# 从已有 variant_analysis 的 VCF 输出中提取突变
# → 生成新抗原候选（pVACseq 需要 VEP 注释的 VCF + HLA 分型）

# 简化 demo: 手动选 2-3 个已知乳腺癌新抗原做 mRNA 设计原型
# 参考: PIK3CA H1047R 是已知热点突变
```

#### 第 4-5 天：搭建最小 mRNA 设计原型
```python
# src/drugpipe/rna_vaccine/__init__.py
# src/drugpipe/rna_vaccine/neoantigen_predictor.py  ← 封装 pVACseq
# src/drugpipe/rna_vaccine/mrna_designer.py          ← 密码子优化 + UTR

# 最小原型: 给定一个肽段 → 反向翻译 → 密码子优化 → 加 UTR → 输出 mRNA 序列
```

这样研究所 B 的人到了可以直接在你的骨架上继续。

### Q: 同一个 GitHub 仓库不会 README/pyproject.toml 乱套吗？

**合理的担心**。推荐 **monorepo + 子包结构**：

```
T2Lead/                           # 顶层: 你维护
├── README.md                     # 项目总览（你写）
├── CONTRIBUTING.md               # 贡献指南（你写）
├── pyproject.toml                # 顶层包: drugpipe（你维护）
├── configs/                      # 统一配置（你维护 schema）
│
├── src/drugpipe/
│   ├── pipeline.py               # 你维护
│   ├── treatment_router.py       # 你维护
│   │
│   ├── omics_integration/        # 研究所 A 的领地
│   │   ├── README.md             # A 自己的模块文档
│   │   └── ...
│   │
│   ├── rna_vaccine/              # 研究所 B 的领地
│   │   ├── README.md             # B 自己的模块文档
│   │   └── ...
│   │
│   ├── target_to_hit/            # 研究所 C 的领地
│   ├── hit_to_lead/              # 研究所 C
│   ├── lead_optimization/        # 研究所 C
│   └── drugreflector_bridge.py   # 研究所 C
│
└── docs/
    ├── interfaces/               # 接口定义（你维护）
    │   ├── target_report_schema.yaml
    │   ├── vaccine_report_schema.yaml
    │   └── drug_report_schema.yaml
    └── ...
```

**关键原则**:
1. **pyproject.toml 只有一份**，由你维护。各研究所在自己的 `requirements-X.txt` 里声明额外依赖
2. **README.md 顶层**由你写，各模块有自己的子 README（不冲突）
3. **CODEOWNERS 文件**指定目录归属，PR 自动 assign reviewer
4. **分支保护**: `main` 只有你能合并，`dev` 三方 PR

```
# .github/CODEOWNERS
src/drugpipe/omics_integration/     @institute-a-lead
src/drugpipe/target_discovery/      @institute-a-lead
src/drugpipe/rna_vaccine/           @institute-b-lead
src/drugpipe/target_to_hit/         @institute-c-lead
src/drugpipe/hit_to_lead/           @institute-c-lead
src/drugpipe/lead_optimization/     @institute-c-lead
src/drugpipe/pipeline.py            @you
configs/                            @you
docs/interfaces/                    @you
```

### Q: 429 封禁风险 — 几十个靶点 × 成百上千变体 → 数千请求

**确实是问题**。当前代码的 `structure_polite_sleep: 0.05` (50ms) 太激进。

**已实现的缓解措施** (§10 代码已落地):
1. ✅ **全局结构缓存**: 同一 `UniProt+突变` 不会重复请求
2. ✅ **ESMFold 远程默认关闭**: `remote_enabled: false`
3. ✅ **长序列直接跳过**: 不再对 >800aa 序列发远程请求

**额外建议** (还没改但应该做):
- 将 `structure_polite_sleep` 改为 `0.5` (500ms) 对 RCSB
- RCSB 搜索结果（按 UniProt → PDB IDs）已有缓存 `_pdb_cache`
- 下载的 PDB 文件已有磁盘缓存 (`if pdb_path.exists(): return`)

---

## 13. 你（总负责人）的完整任务清单

### 13.1 角色定位

你是 **架构师 + 集成者 + 进度管理者**。不直接做三个研究所的核心模块，但负责：
1. 框架搭建和接口定义
2. 三轨合并和最终推荐输出
3. 进度追踪和技术仲裁
4. mRNA 疫苗方向的先行探路（直到研究所 B 的人到位）

### 13.2 你的具体任务表

#### Phase 0：立即 (本周)

| # | 任务 | 状态 | 说明 |
|---|------|------|------|
| P0-1 | §10 四项基础设施改进 | ✅ 已完成 | UniProt 缓存 / ESMFold 分级 / 全局去重 / 配置 |
| P0-2 | 下载 Swiss-Prot FASTA + 人类 idmapping | 待做 | `wget` 两个文件，~300 MB |
| P0-3 | ESMFold 本地部署 (LiteFold Docker) | 待做 | 4090 能跑 ≤700aa |
| P0-4 | 配置验证：用 Xena 乳腺癌数据跑一轮确认 §10 改进生效 | 待做 | 对比日志：缓存命中数 / 跳过数 |

#### Phase 1：接口定义 + 探路 (下周)

| # | 任务 | 状态 | 说明 |
|---|------|------|------|
| P1-1 | 定义 `target_report.yaml` JSON Schema | 待做 | 研究所 A 的输出格式 |
| P1-2 | 定义 `vaccine_report.yaml` JSON Schema | 待做 | 研究所 B 的输出格式 |
| P1-3 | 定义 `drug_report.yaml` JSON Schema | 待做 | 研究所 C 的输出格式 |
| P1-4 | 创建 `src/drugpipe/rna_vaccine/` 骨架 | 待做 | 最小原型：肽段 → mRNA 序列 |
| P1-5 | 安装 pVACtools + ViennaRNA 验证 | 待做 | mRNA 探路用 |
| P1-6 | 创建 `src/drugpipe/treatment_router.py` 骨架 | 待做 | 三轨分流逻辑 |

#### Phase 2：集成层开发 (第 2-3 周)

| # | 任务 | 状态 | 说明 |
|---|------|------|------|
| P2-1 | `treatment_router.py` 完整实现 | 待做 | 可成药性判断 + 新抗原评估 → 分流 |
| P2-2 | `pipeline.py` 重构：支持三轨调度 | 待做 | 现有 4-stage → 多轨编排 |
| P2-3 | `patient_aggregation` 扩展 | 待做 | 合并 drug_report + vaccine_report |
| P2-4 | Track C 放疗推荐知识库 | 待做 | 规则+文献，不需要计算 |
| P2-5 | CODEOWNERS + 分支保护配置 | 待做 | GitHub 权限管理 |

#### Phase 3：联调 (第 4-5 周)

| # | 任务 | 状态 | 说明 |
|---|------|------|------|
| P3-1 | 三方接口联调（模拟数据） | 待做 | 用 Xena 乳腺癌做端到端测试 |
| P3-2 | 综合治疗推荐报告生成 | 待做 | HTML/PDF 可视化报告 |
| P3-3 | 性能测试 + 瓶颈优化 | 待做 | 确认本地缓存/ESMFold 的加速效果 |

### 13.3 你本周的执行清单 (按天)

**周三 (今天剩余)**:
- [x] §10 代码已落地
- [ ] 下载 Swiss-Prot: `wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -P /root/autodl-fs/uniprot/`
- [ ] 下载 idmapping: `wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz -P /root/autodl-fs/uniprot/`
- [ ] 更新 `default_config.yaml` 中 `uniprot_local.db_dir` 指向下载目录

**周四**:
- [ ] ESMFold 本地部署: `git clone https://github.com/LiteFold/litefold && cd litefold && docker compose up -d`
- [ ] 更新 config: `esmfold.local_url: http://localhost:8000/predict`
- [ ] 用 Xena 乳腺癌跑一轮验证

**周五**:
- [ ] 安装 pVACtools: `conda create -n pvac python=3.10 && conda activate pvac && pip install pvactools`
- [ ] 安装 ViennaRNA: `conda install -c bioconda viennarna`
- [ ] 创建 `src/drugpipe/rna_vaccine/__init__.py` 骨架
- [ ] 把这份文档发给三个研究所的负责人

**下周一-二**:
- [ ] 定义三份 YAML Schema (P1-1/2/3)
- [ ] 创建 `treatment_router.py` 骨架 (P1-6)
- [ ] 研究所 B 的人到位 → 交接 mRNA 探路成果

---

> **下一步**: 
> 1. 是否现在就执行 UniProt 下载和 ESMFold 本地部署？
> 2. 三个研究所的负责人 GitHub 账号是什么？我帮你配 CODEOWNERS
> 3. 是否需要我先创建 `rna_vaccine/` 的最小骨架代码？
