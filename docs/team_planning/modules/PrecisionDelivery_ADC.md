# PrecisionDelivery 模块任务书（精密递送 / 靶向偶联家族）

**负责人：** 新人（蛋白质工程 / AI 蛋白设计方向）  
**模块定位：** 面向肿瘤细胞表面异常蛋白的"靶向投递"家族，第一阶段聚焦 ADC（抗体-药物偶联物），第二/三阶段扩展到 SMDC 和 BiTE。  
**一句话目标（第一阶段）：** 拿 BioDriver 的胞外靶点清单，产出 ADC 候选（抗体序列 + linker + payload）并送 Simulation Hub 做复合物稳定性验证。

---

## 0) 一分钟导读

> **这个模块在干什么？** 给肿瘤细胞"装定位炸弹"。BioDriver 帮你选好了"哪个表面蛋白只肿瘤有，正常组织没有"—— 你负责给它设计一个抗体（导航系统）+ linker（引信）+ payload（炸弹），让药物只在肿瘤里释放毒性。第一阶段只做 ADC（成药最成熟），之后再扩到 SMDC（小分子做导航）和 BiTE（招 T 细胞来杀）。
>
> - **上游**：BioDriver（`surface_targets.csv` + 胞外域 PDB）
> - **下游**：Simulation Hub（抗体-抗原复合物 PDB，做蛋白-蛋白 MD）
> - **第一阶段 smoke 目标**：复现 1 个已知 ADC（如 T-DM1 / T-DXd）做 positive control，再为 BioDriver 提供的 1 个新靶点设计 3-5 条候选抗体。
> - **关键风险**：AI 抗体设计"可成药性"仍开放问题；payload 脱靶毒性无法纯算法预测；第一阶段**不做**递送/药代药动学。

---

## 1) 输入数据（从 BioDriver 来）

```
deliveries/<run_id>/to_adc/
├── surface_targets.csv        # gene_symbol, uniprot_id, surface_confidence, tumor_expr, normal_expr, lfc
├── ectodomains/<uniprot>.pdb  # 每个靶点的胞外域（extracellular domain）3D 结构
└── meta.json
```

---

## 2) 中间过程（工具链，允许自主选型）

### 2.1 靶点质量把关（Target Qualification）
- 胞外域必须存在（跨膜或 GPI 锚定的才适合 ADC）
- 肿瘤高表达、正常组织低表达（避免正常组织毒性）
- 内吞能力（ADC 依赖抗体-抗原复合物被内吞到溶酶体释放 payload）
  → 查文献 / UniProt subcellular location

### 2.2 抗体设计
- 从头生成：RFantibody（RoseTTAFold 系列）/ DiffAb / IgFold
- 结构预测：AlphaFold-Multimer / IgFold
- CDR 优化：SAbPred / ABlooper

### 2.3 抗原-抗体对接
- HADDOCK / ClusPro / AlphaFold-Multimer 直接生成复合物
- 表位（epitope）定位：对肿瘤特异表位（mutant/glycosylated）打标签

### 2.4 Linker + Payload 选择
- **Payload 候选**：MMAE / MMAF / DM1 / 拓扑异构酶 I 抑制剂（DXd）
- **Linker 类型**：可裂解（cleavable, Val-Cit）vs 不可裂解
- **DAR（Drug-Antibody Ratio）**：标准 2 / 4 / 8，建议 4 作为默认

### 2.5 开发性评估（Developability）
- 聚集倾向、等电点、糖基化位点、免疫原性（免疫原性工具：DeepImmuno / NetMHCpan II）

---

## 3) 输出数据

### 3.1 结果文件（给 Lead 的报告产物）

```
results/<run_id>/
├── adc_candidates.csv          # candidate_id, target, ab_vh_seq, ab_vl_seq, payload, linker, dar, developability_score
├── complexes/<candidate>_Ag_Ab.pdb
├── payloads/<payload_id>.sdf
├── linkers/<linker_id>.sdf
├── figures/
│   ├── binding_interface.png
│   └── target_expression_matrix.png
└── REPORT.md
```

### 3.2 给 Simulation Hub 的交付包

你这条路**要给 SimHub 提交两种包**：

**3.2.1 抗体-抗原复合物包（主 MD）—— 走 `antibody_antigen` 分支**

**完整字段级契约 → [`../SIMHUB_CONTRACT.md` §4 分支 C](../SIMHUB_CONTRACT.md#4-分支-cantibody_antigen)。**

速查：

```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb           # Chain H (VH) + Chain L (VL) + Chain A (抗原胞外域)
│                         # BiTE 额外：Chain H2/L2 + Chain T (CD3ε)
│                         # CDR 6 条 loop 不可缺失；SSBOND 必须完整
├── cdr_annotation.json   # 可选
└── meta.json             # molecule_type: "antibody_antigen"
```

> ⚠️ 抗体是大蛋白，**不要传 SDF**。整条抗体塞进 SDF + AM1-BCC 直接让 `sqm` 跑爆。

**3.2.2 payload + linker 小分子包（独立评估）—— 走 `small_molecule` 分支**

> 原因：抗体 ~150 kDa + payload ~1 kDa 的系统规模差 150 倍，塞进一次 MD 时间/空间尺度不匹配。分开做。

**完整字段级契约 → [`../SIMHUB_CONTRACT.md` §2 分支 A](../SIMHUB_CONTRACT.md#2-分支-asmall_molecule)。**

速查：

```
deliveries/<run_id>/to_simhub/<payload_case_id>/
├── receptor.pdb   # 可选：若评估 payload 与靶点（如 tubulin）结合
├── ligand.sdf     # payload 或 payload-linker-lys 3D
└── meta.json      # molecule_type: "small_molecule"; parent_hit_id: <抗体 case_id>
```

> 第一阶段只做 3.2.1 的抗体-抗原 MD，payload 路径等抗体验证通过后启动。

---

## 4) 建议执行步骤

1. **基线验证**：挑 1 个文献 ADC（如 Trastuzumab-DM1 / HER2-DXd）作为 positive control，跑通全流程。
2. **靶点筛选**：从 BioDriver 清单里，按"胞外 + 高特异性 + 可内吞"筛出 Top 3-5 个。
3. **抗体从头设计**：对每个靶点设计 3-5 条候选抗体 VH/VL。
4. **复合物预测 + 界面分析**：挑界面最合理的送 Simulation Hub。
5. **Linker/Payload 配对**：对 Top 候选给出 2-3 种组合。

---

## 5) 交付物

- 独立 GitHub 仓库（含 README / 环境 / 运行命令）
- `REPORT.md`：5-10 页等价
- 上述结果文件 + Simulation Hub 交付包

---

## 6) 验收标准

- **可复现性**：新环境按 README 可完整跑通。
- **完整性**：候选 ADC 设计 + 复合物结构 + linker/payload 三件套齐全。
- **可解释性**：每个 Top 候选要有靶点合理性 + 界面合理性 + 开发性打分。
- **Positive Control**：至少有 1 个已知 ADC 通过流程且被识别为"高分"。

---

## 7) 风险与边界

- ADC 毒性主要来自 payload 脱靶，ADC 设计工具无法完全预测，要**明确写入风险章节**。
- AI 抗体设计的"可成药性"仍是开放问题，湿实验验证必要。
- **不做**：小分子药物筛选、siRNA 设计、mRNA 疫苗、分子动力学本身（走 Simulation Hub）。

---

## 8) 可探索方向（Roadmap：ADC → SMDC → BiTE 家族扩展）

### 8.1 第一阶段：ADC（核心，Must Do）
- 已在上文展开，走抗体 + linker + payload 路线。

### 8.2 第二阶段：SMDC（Small Molecule Drug Conjugate，Should Do）

> SMDC = 把 ADC 里的"抗体导航"换成"小分子导航"。合成成本低、组织穿透更好；但特异性不如抗体，需要导航配体对靶点有**天然高亲和**。

- **经典案例**：叶酸-药物偶联（folate receptor α 肿瘤过表达）、PSMA-617（前列腺癌，已上市 Pluvicto）
- **与 ADC 的差别**：
  - 导航：小分子配体（folate / PSMA 配体 / 整合素 RGD / ...）
  - Linker：常用可裂解酯键 / 腙键（pH 敏感）
  - Payload：同 ADC（MMAE / DM1 / DXd 等）
  - 不依赖内吞（部分 SMDC 在细胞外切割后扩散入胞）
- **工具链**：
  - 配体选择参考 ChEMBL 已知 binder + 复现 [`src/drugpipe`](../../../src/drugpipe) 的小分子 pipeline
  - 配体-payload 连接的 3D 构象用 AutoDock Vina / RDKit

### 8.3 第三阶段：BiTE（Bispecific T-cell Engager，Could Do）

> BiTE = 双特异抗体。一条臂抓肿瘤抗原，另一条臂抓 T 细胞上的 CD3 → **不靠 payload 杀伤，靠把 T 细胞拉到肿瘤旁边杀伤**。机制完全不同于 ADC/SMDC，更像免疫疗法。

- **经典案例**：Blinatumomab（CD19×CD3，ALL 白血病）、Tarlatamab（DLL3×CD3，SCLC）
- **与 ImmunoGen 的关系**：都属免疫疗法，但 BiTE 是"现成抗体架构"，ImmunoGen 是"让病人自己的免疫系统识别新抗原"。在本项目里 BiTE 归 PrecisionDelivery，**不归 ImmunoGen**（因为核心仍是抗体工程 + 结构优化）。
- **工具链**：
  - scFv 或 knob-into-hole 设计
  - 双臂亲和力平衡（两端不能差太多，否则偏向某侧导致免疫突触不形成）
  - Ab-Ab 相对取向 MD 评估（Sim Hub 蛋白-蛋白 MD 可复用）
- **风险**：细胞因子风暴（CRS）、T 细胞耗竭 —— 纯算法无法预测，湿实验必要。

### 8.4 通用方向（任意阶段）
- Bispecific ADC（抗体部分也做双特异）
- Conditional-active（肿瘤微环境 pH / 蛋白酶响应）
- 自动化多靶点批处理 + 开发性综合评分系统
- 和 BioDriver 的 "surface_targets.csv" 联动：按"最佳疗法线"自动分流（适合 ADC / SMDC / BiTE 中哪一个）
