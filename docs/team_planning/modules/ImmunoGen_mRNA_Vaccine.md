# ImmunoGen 模块任务书（mRNA 疫苗 / 免疫基因）

**负责人：** 新人（免疫学 / AI 序列工程方向）  
**模块定位：** 基于病人特异性新抗原（Neoantigen）设计多价 mRNA 疫苗。  
**一句话目标：** 拿 BioDriver 的新抗原候选 + HLA 分型，输出候选 mRNA 全序列（含 UTR、密码子优化、多价串联设计）。

---

## 0) 一分钟导读

> **这个模块在干什么？** 做**病人个性化 mRNA 肿瘤疫苗**。从 BioDriver 拿到"这个病人哪些突变肽能被他自己 HLA 递呈"的候选清单，挑免疫原性最高的 10-20 条突变多肽，用 linker 串起来做成"多价表位"，前后加 UTR、做密码子优化、做二级结构检查，最后输出一条完整的、可直接合成的 mRNA 序列。
>
> - **上游**：BioDriver（`neoantigen_candidates.csv` + `hla_typing.json`）
> - **下游**：Simulation Hub（peptide-MHC 复合物 PDB，做结合稳定性 MD）
> - **第一阶段 smoke 目标**：跑通 1 个病人 → 挑 Top-10 neoantigen → 串成 1 条多价 mRNA → 至少 1 个 peptide-MHC 交付包送 Sim Hub。
> - **关键风险**：免疫原性预测（NetMHCpan / DeepImmuno）有假阳；mRNA 二级结构过紧会影响翻译效率；**疫苗设计到 LNP 递送仍有鸿沟**（第一阶段不做 LNP 设计）。

---

## 1) 输入数据（从 BioDriver 来）

```
deliveries/<run_id>/to_immunogen/
├── neoantigen_candidates.csv   # mutation, mut_peptide, wt_peptide, transcript_id, variant_vaf
├── hla_typing.json             # { "HLA-A": [...], "HLA-B": [...], "HLA-C": [...] }
└── meta.json
```

> **HLA 分型不需要你算**，BioDriver 从 WES 用 OptiType/HLA-HD 算好给你。

---

## 2) 中间过程（工具链）

### 2.1 表位筛选（Epitope Prediction）
- MHC-I 结合预测：**NetMHCpan-4.1** / **MHCflurry-2.0** / **BigMHC**
- MHC-II 结合预测（CD4 辅助）：NetMHCIIpan
- 免疫原性（Immunogenicity）评分：DeepImmuno / PRIME / Repitope
- T 细胞表位工具集：IEDB Analysis Resources

### 2.2 候选多肽排序
- 综合打分：`rank_score = w1·HLA_affinity + w2·immunogenicity + w3·VAF + w4·wt_peptide_dissimilarity`
- 必须过滤：与 wild-type 过于相似的（避免中枢免疫耐受）

### 2.3 多价串联设计（Multivalent Design）
- 选 Top 10-20 peptide，通过 linker 串联（常用 AAY、GS-linker）
- 加入信号肽（MITD 或经典信号肽）+ 跨膜/胞内定位序列（可选）

### 2.4 mRNA 全序列设计
- 5' UTR：HBA1 / HBB / 优化后的合成 UTR
- 3' UTR：α-globin 系列 + 多 A 尾
- 密码子优化：**LinearDesign**（清华）/ COOL / CodonOpt
- 修饰碱基策略（设计说明）：N1-methylpseudouridine

### 2.5 二级结构与稳定性
- RNAfold / LinearFold 看 MFE、GC 含量
- （可选）mRNA 稳定性预测工具：Saluki / RNAsnp

### 2.6 （可选）LNP 递送方案
- 综述 SM-102 / ALC-0315 等脂质组合（第一阶段可不深入）

---

## 3) 输出数据

### 3.1 结果文件（给 Lead 的报告产物）

```
results/<run_id>/
├── peptide_mhc_ranking.csv     # peptide, hla_allele, affinity_nM, immunogenicity, rank_score
├── selected_peptides.csv        # 最终进入多价疫苗的 top-N peptide
├── mrna_vaccine.fasta           # 完整 mRNA 序列（多价串联后）
├── mrna_design.json             # 5'UTR / ORF / 3'UTR / poly-A / 修饰方案说明
├── figures/
│   ├── binding_affinity_heatmap.png
│   └── mrna_secondary_structure.png
└── REPORT.md
```

### 3.2 给 Simulation Hub 的交付包

**你走的是 `molecule_type = peptide_mhc` 分支，纯蛋白多链方案，禁止使用 SDF**。

**完整字段级契约 → 请严格按 [`../SIMHUB_CONTRACT.md` §3 分支 B](../SIMHUB_CONTRACT.md#3-分支-bpeptide_mhc) 交付。**

速查：

```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb         # Chain M (MHC α) + Chain B (β2m) + Chain P (peptide)
│                       # 初始构象用 AlphaFold-Multimer 或 PANDORA 生成
├── hla_allele.txt      # 可选，如 "HLA-A*02:01"
└── meta.json           # molecule_type: "peptide_mhc"
```

> ⚠️ 强行用 SDF + AM1-BCC 会让 SimHub 的 `sqm` 算几小时后崩溃，物理结果全错。多肽必须保留氨基酸残基拓扑，走 AMBER14SB 纯蛋白力场。

> Simulation Hub 对 peptide-MHC 走标准 AMBER14SB 蛋白力场 MD，评估：
> - 界面 RMSD（peptide 是否稳定结合槽内）
> - 氢键/盐桥保持
> - 相对结合自由能（对比 wt_peptide 作 Baseline）

---

## 4) 建议执行步骤

1. **基线验证**：挑 1 个已知免疫原性强的 neoantigen（如 KRAS G12D 来源的经典 peptide）跑通全流程。
2. **真实数据**：对 BioDriver 给的 Top 30-50 候选 peptide 做 MHC 结合预测 + 免疫原性评分。
3. **筛选**：按 `rank_score` 取 Top 10-20 进入多价疫苗。
4. **多价设计**：串联成一条完整 ORF，做密码子优化 + UTR 拼接 + 二级结构检查。
5. **送 Simulation Hub**：选 Top 3-5 peptide-MHC 做 MD 稳定性验证。

---

## 5) 交付物

- 独立 GitHub 仓库（含 README / 环境 / 运行命令）
- `REPORT.md`：5-10 页等价
- 上述所有结果文件 + Simulation Hub 交付包

---

## 6) 验收标准

- **可复现性**：新环境按 README 可完整跑通。
- **完整性**：peptide-MHC 排名 + 多价串联设计 + 完整 mRNA 三件套齐全。
- **可解释性**：Top 候选要有 HLA / 免疫原性 / VAF 三方面证据。
- **Positive Control**：至少有 1 个已知 neoantigen 通过流程且被识别为"高分"。

---

## 7) 风险与边界

- MHC-II 预测准确率显著低于 MHC-I，要在报告中注明。
- AI 免疫原性预测仍有较大不确定性，最终要湿实验（ELISpot / 四聚体染色）。
- **不做**：小分子药物筛选、siRNA 设计、ADC 设计。

---

## 8) 与任务书其他模块的关系

- 这条线**才是**"把 Top 突变串进一条 mRNA 做多价疫苗"的落脚点；DrugReflector 是小分子老药新用，两者不要混。

---

## 9) 可探索方向（Roadmap）

- 个体化 HLA-II 深度优化（CD4 辅助信号）
- 自体树突细胞疫苗联合路线
- 自动化多病人批处理 + 临床级质控
