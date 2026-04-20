# RNAi-Therapy 模块任务书（RNA 干扰疗法）

**负责人：** 新人（核酸药物 / 生信方向）  
**模块定位：** 基于致病基因设计 siRNA 候选序列。  
**一句话目标：** 拿 BioDriver 的致病靶点清单 + mRNA 序列，产出高质量 siRNA 候选 + 脱靶风险评估。

---

## 0) 一分钟导读

> **这个模块在干什么？** 设计"让致病基因闭嘴"的 siRNA 药物。给定 BioDriver 挑出的致病基因 + 对应 transcript 序列，在 mRNA 上找 21nt 最佳切割位点、做 BLAST 脱靶、做二级结构可及性检查、再给化学修饰建议，输出一批候选 siRNA 序列。
>
> - **上游**：BioDriver（`pathogenic_targets.csv` + transcripts FASTA）
> - **下游**：Simulation Hub（第一阶段**可选**，标准 MD 对 siRNA-Ago2 复合物成本高，第一阶段只交脱靶 + 可及性报告给 Lead 即可）
> - **第一阶段 smoke 目标**：跑通 1 个靶点 → 输出 Top-10 siRNA 候选 CSV（序列 + 修饰 + 脱靶评分）。
> - **关键风险**：siRNA 脱靶不是靠一个工具就能杜绝的，必须多工具交叉；递送（LNP/GalNAc）是现实瓶颈，但这阶段**不要求**做递送设计。

---

## 1) 输入数据（从 BioDriver 来）

```
deliveries/<run_id>/to_rnai/
├── pathogenic_targets.csv     # gene_symbol, ensembl_gene_id, canonical_transcript_id, priority_score, reason
├── transcripts.fasta          # 对应 transcript 的 mRNA 序列
└── meta.json
```

**不需要也不要**直接对接 BGI。原始数据问题统一提给 BioDriver 负责人。

---

## 2) 中间过程（工具链，可自主扩展）

### 2.1 siRNA 候选生成
- 规则型：Reynolds / Ui-Tei / Amarzguioui 规则
- 在线/本地工具：siDESIGN Center、DSIR、siRNA Selection Server
- ML 评分：DeepSHAP / siVirus / SMEsiR 等

### 2.2 脱靶风险评估（Off-target）
- BLAST / Bowtie 对人类全转录组做 seed region 匹配
- miRNA-like 机制评估（3' UTR seed match）
- 必要时：人类 RefSeq / GENCODE 全 transcript DB 索引

### 2.3 二级结构分析
- RNAfold（ViennaRNA）对目标 mRNA 做二级结构预测
- 避开高度折叠区域（降低可及性）

### 2.4 化学修饰建议
- 2'-OMe / 2'-F / LNA 位点建议
- PS 骨架修饰降解抗性

### 2.5 （可选）递送与剂量
- LNP 递送策略综述（本模块第一阶段可不深入）

---

## 3) 输出数据

### 3.1 结果文件（给 Lead 的报告产物）

```
results/<run_id>/
├── sirna_candidates.csv        # id, target_gene, sense_seq, antisense_seq, position_on_mrna, score, gc_content
├── offtarget_report.csv        # sirna_id, potential_offtarget_genes, seed_match_count, risk_level
├── secondary_structure/<gene>_mrna_fold.ps
├── figures/
│   └── target_accessibility.png
└── REPORT.md
```

### 3.2 给 Simulation Hub 的交付包（第一阶段**可选**）

> 第一阶段**建议不走标准 MD**，因为核酸力场复杂且优先级较低。  
> 改为把"**脱靶 + 二级结构评估 + siRNA 评分**"打包交给 Lead 直接复核。

如果后期做 siRNA-Ago2 装载 MD，再按下面契约提交：

**你走的是 `molecule_type = sirna_ago2` 分支。第一阶段默认不进 SimHub**（RNAi 交脱靶 + 二级结构给 Lead 即可）。

**第二阶段启用时的完整字段级契约 → [`../SIMHUB_CONTRACT.md` §5 分支 D](../SIMHUB_CONTRACT.md#5-分支-dsirna_ago2第二阶段启用)。**

速查（第二阶段才需要）：

```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb         # Chain S (sense RNA) + Chain A (antisense RNA) + Chain P (Ago2)
│                       # RNA 残基名用 A/U/G/C；化学修饰用 HETATM + CONECT
├── modifications.json  # 每个修饰位置的化学结构 + SMILES
└── meta.json           # molecule_type: "sirna_ago2"
```

> ⚠️ siRNA 是核酸不是小分子，**不要传 `ligand.sdf`**。必须用 AMBER OL3（核酸）+ AMBER14SB（蛋白）组合，走 `complex.pdb` 方案。

---

## 4) 建议执行步骤

1. **基线验证**：挑 1 个已知有效 siRNA（文献报道过的，如 Patisiran 的 TTR siRNA）跑通全流程，作为 positive control。
2. **真实靶点**：对 BioDriver 给的 Top 3-5 个靶点分别设计 5-10 条候选 siRNA。
3. **排序**：按"活性预测分数 × (1 - 脱靶风险) × 可及性"综合排名。
4. **报告**：REPORT.md 写明为什么选这些靶点、为什么选这些 siRNA。

---

## 5) 交付物

- 独立 GitHub 仓库（含 README / 环境 / 运行命令）
- `REPORT.md`：5-10 页等价内容
- 上述 CSV / 二级结构图 / 合理性分析

---

## 6) 验收标准

- **可复现性**：新环境按 README 可完整跑通。
- **完整性**：候选 siRNA + 脱靶评估 + 二级结构三件套齐全。
- **可解释性**：Top 候选要有打分依据，不是只列一堆序列。
- **Positive Control**：至少有 1 个已知 siRNA 作为 sanity check 且被识别为"高分"。

---

## 7) 风险与边界

- 脱靶数据库要与 BioDriver 用的人类基因注释版本一致（GENCODE 版本号写进 meta）。
- siRNA 活性预测模型泛化能力有限，Top-N 结果最终要湿实验验证。
- **不做**：小分子药物筛选、蛋白结构预测、ADC 设计、免疫表位预测。

---

## 8) 可探索方向（Roadmap）

- ASO（反义寡核苷酸）模块作为补充
- crisprRNA / sgRNA 设计作为衍生路线
- 自动化多靶点批处理
