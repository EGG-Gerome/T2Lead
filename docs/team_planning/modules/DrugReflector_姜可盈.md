# DrugReflector 模块任务书（姜可盈）

**负责人：** 姜可盈  
**模块定位：** 独立表型药物发现模块（老药新用 / 小分子重定位），并作为"智能衍生"的 hit 来源入口  
**一句话目标：** 拿 BioDriver 的疾病表达签名，用 [Cellarity DrugReflector](https://github.com/Cellarity/drugreflector) 做小分子候选排序，产出 Top-N 候选清单；对其中进入 Simulation Hub 后表现不足的 hit，后续触发 Level 2 类似物生成。

---

## 0) 一分钟导读

> **这个模块在干什么？** 小分子"表型药物发现"。不靠某个特定靶点，而是让单细胞疾病签名直接"照镜子"—— 在 DrugReflector 这个现成模型里找最能把疾病状态扭回正常的老药，排序输出 Top-N 候选，挑其中 3-5 个送 Simulation Hub 做 MD 验证。
>
> - **上游**：BioDriver（单细胞预处理 .h5ad + disease_signature.csv）
> - **下游**：Simulation Hub（小分子 + 靶点 PDB 的验证包）
> - **第一阶段 smoke 目标（L1 只做老药新用）**：跑通 1 个病人的签名 → DrugReflector Top-N → 给 Sim Hub 至少 1 个合格交付包。
> - **第二阶段（L2 智能衍生）**：对 Sim Hub 判定"弱"的 hit，基于 scaffold 触发类似物生成（引擎沿用 T2Lead 资产，Lead 维护）。
> - **关键风险**：DrugReflector 只认 978 landmark genes，HGNC 映射不全会掉覆盖度；商业授权受限（学术 OK）。

---

## 0.1) 分层能力（Level 1 → Level 2）

模块分两层能力逐步上线，**不要混在一起做**：

| Level                              | 能力                                                         | 交付范围                               | 阶段             |
| ---------------------------------- | ------------------------------------------------------------ | -------------------------------------- | ---------------- |
| **L1 · 重定位（Repositioning）**   | 用单细胞签名匹配现成小分子（DrugReflector）                  | Top-N 老药候选 + 给 Sim Hub 的交付包   | **第一阶段必达** |
| **L2 · 智能衍生（AI Derivation）** | 基于 L1 的 hit，用**骨架感知生成 + pIC50 ML + ADMET + MPO**（复用 T2Lead 现有资产）产出类似物 | 每个弱 hit 对应一批类似物的 SDF + 打分 | 第二阶段拓展     |

**触发逻辑（由 Simulation Hub 驱动，自动闭环）：**

```
DrugReflector L1 → Top-N hits → Simulation Hub MD
                                     │
                ┌────────────────────┴──────────────┐
                ▼                                    ▼
      ΔG 达标 & RMSD 稳定                    ΔG 不足 / RMSD 漂移
            │                                        │
     进"决策级候选档案"                      触发 L2：Analog Generation
     → 广科医                                 （以该 hit 为 scaffold）
                                                     │
                                                     ▼
                                          生成类似物 → 回送 Simulation Hub
```

> **第一阶段你只要做 L1**。L2 的引擎沿用 T2Lead 原有代码（`src/drugpipe/target_to_hit`、`hit_to_lead`、`lead_optimization`），由 Lead 协助维护；你只需要保证 L1 输出里带上"骨架可追溯"字段（SMILES、母骨架标识），让 L2 引擎能拿来启动生成。

---

## 1) 输入数据（从 BioDriver 来，不直接对接 BGI）

你只需要 BioDriver 交付给你的这一包：

```
deliveries/<run_id>/to_drugreflector/
├── preprocessed.h5ad
├── disease_signature.csv     # gene_symbol, v_score, lfc, adj_p
└── meta.json
```

> 你自己也可以对 `preprocessed.h5ad` 做"最后一公里"的签名细化（因为 DrugReflector 更偏向用**全基因 v-score 向量**，而不是仅一个短 DEGs 列表）。

**BGI 数据需求清单不由你对接**，参见 [`../BGI_DATA_REQUEST.md`](../BGI_DATA_REQUEST.md)，只作为背景阅读。

---

## 2) 中间过程（工具链）

- 读入：Scanpy 读 `.h5ad`；pandas 读 DEA CSV
- DEA / 签名：可以用 BioDriver 给的，也可以用 `dr.compute_vscores_adata()` 自己重算
    （两者对比做一致性检查，作为健壮性证据）
- 通路富集（合理性检查）：`gseapy` / `enrichr`
- DrugReflector 推理：
    - 下载 checkpoints（Zenodo DOI 10.5281/zenodo.16912444）
    - `dr.DrugReflector(checkpoint_paths=...).predict(vscores, n_top=50)`
    - **必须做** gene name HGNC 规范化 + 基因覆盖度检查（`check_gene_coverage`）
- 结果注释：对 Top-N 化合物查 DrugBank / ChEMBL / PubChem，补机制与适应症

---

## 3) 输出数据（给 Simulation Hub 与 Lead）

### 3.1 给 Lead 的报告产物

```
results/<run_id>/
├── DEA_results.csv              # gene, log2FC, p-value, adj_p-value
├── disease_signature.csv        # 用于 DrugReflector 的签名向量
├── drugreflector_top_hits.csv   # compound_id, rank, score, prob, mechanism_note,
│                                # smiles, scaffold_smiles  ← L2 所需
├── figures/
│   ├── umap.png
│   ├── volcano.png
│   ├── pathway_enrichment.png
│   └── top_drugs_ranking.png
└── REPORT.md
```

> **L2 前置要求**：`drugreflector_top_hits.csv` 必须含 `smiles` 与 `scaffold_smiles`（用 RDKit 的 Murcko scaffold 提取即可，一行代码），这是 Level 2 衍生引擎唯一需要的额外字段。

### 3.2 给 Simulation Hub 的交付包（每个要验证的小分子一个目录）

**你走的是 `molecule_type = small_molecule` 分支**，是 4 条下游中**唯一合法使用 `ligand.sdf`** 的路径。

**完整字段级契约 → 请严格按 [`../SIMHUB_CONTRACT.md` §2 分支 A](../SIMHUB_CONTRACT.md#2-分支-asmall_molecule) 交付。**

速查：

```
deliveries/<run_id>/to_simhub/<case_id>/
├── receptor.pdb    # 靶蛋白（清理后）
├── ligand.sdf      # 3D 小分子（已对接、pH 7.4 加氢、白名单元素）
└── meta.json       # molecule_type: "small_molecule"
```

> ⚠️ **如果你看到其他模块（ImmunoGen / ADC / RNAi）也在用 SDF，那是错的**。纯蛋白或核酸系统必须走 `complex.pdb` 方案。

> 注意：你不一定自己做对接。可以先选 Top-5/10 候选，和 Lead 协商由 Simulation Hub 侧做对接 + MD，或你先用 AutoDock Vina 出个初始 3D pose 再送 MD。

---

## 4) 建议执行步骤

1. **数据接收与审计**：拿到 BioDriver 的交付包，检查字段完整性与质量。
2. **基线复跑**：先跑 DrugReflector 官方 PBMC 例子，确认环境、checkpoint、推理正常。
3. **真实数据接入**：用 BioDriver 的 `preprocessed.h5ad` 跑一遍端到端。
4. **结果解释**：选 Top-5/10 做机制注释，挑 1-3 个标记为"最值得送 MD 验证"。
5. **交付 Simulation Hub**：按统一契约打包。

---

## 5) 交付物

- 独立 GitHub 仓库（含 README / 环境 / 运行命令）
- `REPORT.md`（5-10 页等价内容）：数据来源、方法、结果、风险、下一步
- 4 份结果表（上面列过）
- 至少 1 份给 Simulation Hub 的合规交付包

---

## 6) 验收标准

- **可复现性**：新环境按 README 可完整跑通。
- **完整性**：DEA + 签名 + DrugReflector 结果三件套齐全。
- **可解释性**：Top 候选药物有机制 / 文献级解释，不是纯黑盒分数。
- **质量控制**：报告中明确展示 QC 与批次影响处理思路。
- **边界清晰**：声明该模块**不等同**于新抗原疫苗设计（那条线归 ImmunoGen）。

---

## 7) 风险与边界

- **基因覆盖度问题**：DrugReflector 只认 978 landmark genes 的子集；HGNC 映射不全会掉覆盖度，必须做日志与报告。
- **签名质量依赖上游**：如果 BioDriver 分群有问题，你的签名会带偏差 → 在 REPORT 里要和 BioDriver 的 QC 报告交叉引用。
- **商业使用受限**：DrugReflector 许可证禁止商业使用，研究/学术 OK。
- **L1 与 L2 的边界**：L2 是第二阶段能力，**不要第一阶段就去折腾生成模型**，会拖慢端到端打通。第一阶段只要保证 hits 里的 `smiles / scaffold_smiles` 写对。

---

## 8) Level 2（第二阶段拓展，不在第一阶段验收内）

这部分复用 T2Lead 现有代码资产，**你了解即可，不强制你第一阶段实现**。

- **输入**：L1 的 `drugreflector_top_hits.csv` 中任一 hit 的 `scaffold_smiles` + Sim Hub 反馈（ΔG 不足 / RMSD 漂移）
- **引擎**（T2Lead 现成）：
    - 骨架感知生成 / 类似物生成（`src/drugpipe/target_to_hit`）
    - pIC50 回归（RF+MLP，ChEMBL IC50 训练）
    - ADMET / QED 过滤
    - MPO 多目标打分（`src/drugpipe/hit_to_lead/mpo.py`）
- **输出**：每个弱 hit 对应一个 `analogs_<hit_id>/` 目录，内含 SDF + 打分表
- **回送**：新分子按标准契约包入 Simulation Hub 再评估
- **责任**：引擎维护由 Lead 负责；你负责触发与结果解读

---

## 9) 精简版

> 负责搭建一个独立的 DrugReflector 表型药物发现模块：从 BioDriver 给你的单细胞预处理包（`.h5ad` + 疾病签名 CSV）出发，做必要的 QC / 签名细化，运行 DrugReflector 输出 Top-N 小分子候选，并为 Top 3-5 准备好 Simulation Hub 的输入契约包。**第一阶段只做 Level 1（老药新用）**；如果某个 hit 在 Sim Hub 里表现不足，第二阶段会触发 Level 2 类似物生成（引擎由 Lead 维护，沿用 T2Lead 资产），你只需在输出里带上 `smiles / scaffold_smiles` 两列以便启动。技术细节自由探索，但要保证可复现、可解释。