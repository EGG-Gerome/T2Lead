# Simulation Hub 模块任务书

**负责人：** 卢振涛  
**模块定位：** 物理仿真 / 验证中台。对 4 下游的所有候选分子做统一 MD + 结合自由能 + 稳定性评估，输出最终可递交临床的候选清单。  
**复用基础：** T2Lead 当前的 Stage4（Docking + MD + RMSD），已完成 AM1-BCC 恢复、MD reliability gating、多 GPU 调度雏形。

---

## 0) 一分钟导读

> **这个模块在干什么？** 当**全团队的"物理质检员"**。4 条下游管线设计出来的分子（小分子 / siRNA / 抗体-抗原 / peptide-MHC），不管用什么工具做的，都到 Sim Hub 里过一遍 MD + ΔG + RMSD，给出"稳定 / 不稳定"、"值得湿实验 / 不值得"的判定，并把最终排名 + 证据打包送给广科医。额外承担 **L2 触发器**：小分子路径里判定"弱"的 hit，自动触发类似物生成（沿用 T2Lead 资产）。
>
> - **上游**：4 下游模块，按 `molecule_type` 分派契约（小分子走 `receptor.pdb + ligand.sdf`；多肽/抗体/siRNA 走 `complex.pdb` 多链方案，**不能用 SDF**）
> - **下游**：广科医（决策级候选档案）+ Analog Generation（阶段二触发）
> - **第一阶段 smoke 目标**：T2Lead Stage4 抽离成独立 `simulation-hub` 仓库 → 4 下游各跑通 1 个 case → 产出 `FINAL_RANKING.csv` + HTML dashboard。
> - **关键风险**：MD / MM-GBSA 是**近似相对排序**而非真值；严格契约是 Sim Hub 活下来的底线，**上游不满足契约直接拒收**，不要擦屁股。

---

## 1) 输入数据（从 4 下游来，**按 molecule_type 走不同契约**）

### 1.0 为什么不能"一个契约走天下"

**核心教训**：`.sdf` 是为**小分子化学药**（Small Molecules）设计的，只记录原子坐标 + 化学键，**不保留氨基酸残基名、链 ID、二级结构拓扑**。

如果把 9-mer 多肽或抗体 Fab 塞进 `ligand.sdf`：
1. 底层 AMBER / GAFF 会把它当成"极其巨大、畸形的未知小分子"
2. 对每个重原子尝试 AM1-BCC 半经验量子化学计算 → **`sqm` 要么跑几小时要么 segfault**
3. 即使算完，静电参数也全错（多肽的静电来自主链 peptide bond + 侧链，GAFF 根本不对）
4. 输出的自由能/轨迹**物理上毫无意义**

所以必须**按 `molecule_type` 分派不同契约**。下面四种分支，SimHub CLI (`simhub run`) 读 `meta.json.molecule_type` 后分派到对应 pipeline；**格式不匹配直接拒收**。

---

### 1.1 分支 A：`small_molecule`

**来源**：DrugReflector（重定位） / AnalogGen（L2 类似物生成）

**契约**：
```
deliveries/<run_id>/to_simhub/<case_id>/
├── receptor.pdb      # 靶蛋白单体/结构域；清理后（无结晶水、原配体、无关金属离子、硫酸根）
├── ligand.sdf        # 3D 小分子；已对接到口袋；pH 7.4 加氢；仅 C/H/O/N/S/P/F/Cl/Br/I
└── meta.json
```

**校验规则（拒收项）**：
- `receptor.pdb` 中不得含 HETATM 配体（除必要金属辅因子如 Zn/Mg 且在 meta 声明）
- `ligand.sdf` 必须是 3D（V2000/V3000），**不接受 2D 坐标**
- `ligand` 质心到 `receptor` 最近原子距离 < 10Å（否则 `E_LIGAND_NOT_DOCKED`）
- 元素白名单：`{C,H,O,N,S,P,F,Cl,Br,I}`；出现其他原子 → `E_UNSUPPORTED_ELEMENT`

**力场**：
- 受体：**AMBER14SB**
- 配体：**OpenFF-2.x**（首选） / **GAFF2**（fallback）+ **AM1-BCC** 电荷
- 水：TIP3P；离子：NaCl 平衡

**SimHub pipeline**（现成 T2Lead Stage4 逻辑）：
```
受体 PDBFixer → 加氢 → 加水 → 加离子
  + 配体 OpenFF/GAFF 参数化 + AM1-BCC 电荷
  → 合并系统 → 能量最小化 → NVT/NPT 平衡 → 生产 MD
  → MM-GBSA + RMSD profile + md_reliable 判定
```

---

### 1.2 分支 B：`peptide_mhc`（纯蛋白-蛋白）

**来源**：ImmunoGen（mRNA 疫苗设计里的免疫原性验证）

**契约**：
```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb       # 纯蛋白，多链：
│                     #   Chain M  → MHC I α chain  (~275 AA)
│                     #   Chain B  → β2-microglobulin (~99 AA)
│                     #   Chain P  → peptide (8-11 AA，通常 9-mer)
│                     # 全部标准氨基酸残基名（ALA/GLY/... 20 种）
│                     # 禁止 UNK / 非标残基 / 糖基化（第一阶段不做糖基化）
├── hla_allele.txt    # 可选，例如 "HLA-A*02:01"（便于报告，不影响 MD）
└── meta.json
```

**坚决摒弃 `.sdf`**。多肽必须保留完整氨基酸残基拓扑。

**校验规则（拒收项）**：
- 链 ID 必须齐全：至少 `M / P` 两条；缺 `P` → `E_CHAIN_ID_MISSING`
- `P` 链残基数在 `[8, 11]` 之间
- 所有残基名必须在 AMBER14SB 标准 20 残基表里（含 HIS 的三种质子化形态 HID/HIE/HIP）
- 不得出现 HETATM（除结构水 HOH，但默认我们去掉）

**力场**：
- **AMBER14SB 纯蛋白力场**
- 水：TIP3P；离子：NaCl
- **绝对不走 GAFF / AM1-BCC**

**SimHub pipeline**（新子流程 `simhub/pipelines/protein_protein.py`）：
```
complex.pdb
  → PDBFixer 补缺失原子（不补大段缺失残基；缺失 > 5 个连续残基直接拒）
  → Modeller 加氢（pH 7.0）
  → 加水盒子 + 离子
  → AMBER14SB 参数化（所有残基直接查参数）
  → 能量最小化 → NVT/NPT 平衡 → 生产 MD
  → 分析：
      - 界面 RMSD（P 链相对于 M 链）
      - P 链 Cα RMSD 时间序列
      - 氢键保持率（P 链与 MHC 结合沟之间）
      - MM-GBSA 结合自由能（把 P 链当"ligand"去做能量分解）
```

**关键评估指标**：peptide 在结合沟里稳定 >100 ns 不解离 → 通过。

---

### 1.3 分支 C：`antibody_antigen`（纯蛋白-蛋白）

**来源**：PrecisionDelivery（ADC / SMDC / BiTE 的抗体部分）

**契约**：
```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb       # 纯蛋白，多链：
│                     #   ADC/SMDC:  Chain H (VH, ~120 AA) + Chain L (VL, ~110 AA) + Chain A (抗原胞外域)
│                     #   BiTE:      额外 Chain H2 + Chain L2 + Chain T (CD3ε 胞外域)
│                     # CDR loop 不可缺失（H1/H2/H3 + L1/L2/L3 共 6 条）
│                     # SSBOND 记录必须完整（抗体至少 4 对二硫键）
├── cdr_annotation.json   # 可选：CDR 残基编号（Kabat 或 IMGT），用于报告可视化
└── meta.json
```

**同样坚决不传 `.sdf`**：ADC 的 payload + linker 是**另外一回事**，走分支 A（small_molecule）独立评估，**不**试图把整个 "抗体+linker+payload" 塞进一次 MD，那不现实（系统太大、时间尺度不匹配）。

**校验规则（拒收项）**：
- 至少 3 条 chain ID（ADC）或 6 条（BiTE）
- `SSBOND` 记录数 ≥ 4；缺失 → `E_SSBOND_MISSING`（抗体没有二硫键会 5 ns 内结构散架）
- CDR H3 残基不能缺失（CDR H3 是主要识别区）
- 抗原链（Chain A）必须是该抗原的 **功能域**，不能是全长膜蛋白（跨膜段走不了水盒子 MD）

**力场**：
- **AMBER14SB 纯蛋白力场**
- **绝对不走 GAFF / AM1-BCC**

**SimHub pipeline**（复用分支 B 的 `protein_protein.py`，只是链 ID 约定不同）：
```
complex.pdb
  → PDBFixer（保留 SSBOND）
  → Modeller 加氢
  → 加水盒子 + 离子
  → AMBER14SB 参数化
  → 能量最小化 → 平衡 → 生产 MD
  → 分析：
      - 界面 RMSD（H+L 相对 A）
      - 抗原-抗体界面原子接触保持率
      - 氢键/盐桥保持率
      - 界面 SASA 变化
      - （可选）MM-PBSA 界面结合自由能分解到每个 CDR
```

**关键评估指标**：CDR-抗原接触在 100 ns 内不解离，界面 RMSD < 3 Å。

---

### 1.4 分支 D：`sirna_ago2`（蛋白-核酸，**第一阶段不跑**）

**来源**：RNAi-Therapy（第一阶段只交脱靶 + 二级结构给 Lead，**不进 SimHub**）

**契约（第二阶段启用时）**：
```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb       # 多链：
│                     #   Chain S → sense RNA 链（反义链的互补，用于建立 duplex）
│                     #   Chain A → antisense RNA 链（真正执行切割的引导链）
│                     #   Chain P → Ago2 PIWI 域 或全长 Ago2
│                     # 残基名：RNA 用 A/U/G/C（不是 DA/DT/DG/DC）
│                     # 修饰位点（2'-OMe / PS / LNA）用 HETATM + CONECT 表示
├── modifications.json     # 每个修饰位置的化学结构描述 + SMILES
└── meta.json
```

**为什么第一阶段不跑**：
- AMBER OL3（RNA）+ AMBER14SB（蛋白）组合的参数兼容性维护成本高
- 化学修饰核苷（2'-OMe / PS / LNA）不在标准残基库里，每个都要手工参数化
- 生物学上：Ago2 切割发生在皮秒-纳秒量级 vs 我们的百纳秒 MD 几乎抓不到构象变化
- ROI 低 → **第一阶段只做脱靶评分（序列比对）+ 二级结构（RNAfold）**，不做 MD

**第二阶段启用条件**：有具体 siRNA 候选表现出"脱靶分析矛盾"需要用结构证据拍板时再开。

---

### 1.5 分派表（SimHub 内部 dispatcher）

| molecule_type | 契约关键文件 | 来源模块 | 力场 | 第一阶段 | 成熟度 |
|--------------|------------|---------|------|---------|--------|
| `small_molecule` | `receptor.pdb` + `ligand.sdf` | DrugReflector / AnalogGen | AMBER14SB + OpenFF + AM1-BCC | ✅ 必做 | P0（T2Lead Stage4 已有） |
| `peptide_mhc` | `complex.pdb` (M/B/P) | ImmunoGen | AMBER14SB 纯蛋白 | ⚠️ 若 ImmunoGen 入选则做 | P1 扩展 |
| `antibody_antigen` | `complex.pdb` (H/L/A) | PrecisionDelivery | AMBER14SB 纯蛋白 | ⚠️ 若 ADC 入选则做 | P1 扩展 |
| `sirna_ago2` | `complex.pdb` + `modifications.json` | RNAi-Therapy | AMBER OL3 + AMBER14SB | ❌ 不跑 MD | P2 后期 |

---

### 1.6 拒收原因码与反馈路径

SimHub 拒收时写 `rejections/<case_id>/reason.md`，包含：

```markdown
case_id: <id>
upstream_module: <module>
rejected_at: <ISO timestamp>
error_code: E_XXX
error_detail: |
  具体哪个文件、哪一行、哪个字段违规
fix_suggestion: |
  上游怎么修（具体到"把 complex.pdb 第 234 行的 UNK 改成 TYR"）
```

拒收列表：详见 [`../ARCHITECTURE_AND_DATAFLOW.md`](../ARCHITECTURE_AND_DATAFLOW.md#36-拒收原因码simhub-统一返回) 第 3.6 节。

---

## 2) 中间过程（工具链）

### 2.1 蛋白准备（Protein Prep）

**第一阶段（锁定现状，稳定跑通）：**

- PDBFixer：补全缺失残基、缺失原子
- OpenMM Modeller：加氢（默认 pH 7.0）、加水盒子、加离子
- 优点：轻量、纯 Python、自动化友好
- 缺点（已知）：质子化状态处理粗糙；影响结合口袋氢键网络

**第二阶段（灰度升级）：**

- 保留 PDBFixer 做缺失残基补全
- 把加氢步骤交给 **PDB2PQR + PROPKA**（基于 pKa 的残基质子化状态判断）
- 只在跑完第一阶段端到端 + 收集到 PDBFixer 失败案例后再切换，避免为理论正确性牺牲稳定性

### 2.2 配体参数化（Ligand Parameterization）——**仅 `small_molecule` 分支走这条**

> ⚠️ **红线**：分支 B / C（peptide_mhc / antibody_antigen）**永远不进这条路径**。纯蛋白-蛋白系统的所有残基都在 AMBER14SB 参数库里，直接查表就行，不需要量子化学电荷。把多肽扔给 AM1-BCC 是前一版设计的重大错误（已于 v0.2 修复）。

**第一阶段（锁定现状）：**

- OpenFF + GAFFTemplateGenerator 作为接口层
- 底层调用 AmberTools（`antechamber` + `sqm`）计算 AM1-BCC 电荷
- 已知痛点：`sqm` 慢 + 脆弱，输入分子稍有扭曲就不收敛 / Segfault

**第二阶段（灰度升级）：Espaloma Fallback**

- 当 `sqm` 报错或超时 → 自动切到 **Espaloma**（GNN-based ML 电荷）
- **关键规则**：任何走 Espaloma fallback 的分子，最终报告必须打上红色标签 `[WARNING: ML_CHARGE]`
- 触发时机：等第一阶段跑通端到端 + 日志里收集到真实坏分子样本，再正式引入

### 2.3 溶剂 & 离子

- **第一阶段**：OpenMM `Modeller.addSolvent`（TIP3P 水 + NaCl）已足够，不替换
- **第二阶段（仅特殊场景）**：
    - 跨膜蛋白模拟时，考虑 PACKMOL 搭膜（POPC/POPE）
    - 这不是 P0，写进 Roadmap

### 2.4 MD 运行与分析

- OpenMM + CUDA 平台，单 GPU baseline + 多 GPU 调度（已在 `feat/stage4-4gpu-xena-autoadapt` 做过雏形）
- 标准流程：能量最小化 → NVT 平衡 → NPT 平衡 → 生产 MD
- 分析：
    - RMSD 稳定性（heavy atom RMSD 轨迹）
    - 结合自由能：MM-GBSA（近似，快）；后期可上 FEP/TI
    - 氢键/盐桥保持率（蛋白-蛋白复合物尤其重要）

### 2.5 质量把关（Gating）

- 已在 T2Lead Stage4 实现：`md_reliable` 逻辑（ΔG 合理范围、RMSD 范围）
- 对不同 `molecule_type` 设不同阈值（ΔG 上限、RMSD 上限）

### 2.6 L2 触发器（Analog Generation，第二阶段启用）

Simulation Hub 不是单纯"打分器"，还是**小分子路径的 L2 触发器**：

```python
if mol["molecule_type"] == "small_molecule" \
   and mol["generation_source"] == "repositioning" \
   and (dG > DG_THRESHOLD  or  rmsd_mean > RMSD_THRESHOLD):
     enqueue_analog_generation(
         parent_hit_id = mol["case_id"],
         scaffold_smiles = mol["scaffold_smiles"],
     )
```

- **引擎**：沿用 T2Lead 现有资产（`target_to_hit` 骨架生成 + `hit_to_lead/mpo.py` 打分 + `lead_optimization` 的 ADMET/QED 过滤）。
- **循环保护**：回送的类似物 `generation_source=analog_generation`，不再进入 L2 判定（最多一次衍生，避免无限循环）。
- **第一阶段不启用**，只要在 ranking CSV 里先预留 `would_trigger_l2=True/False` 标记，方便评估引入 L2 后的预期收益。

---

## 3) 输出数据（给广科医 / Lead 报告）

```
deliveries/<run_id>/simhub_results/<case_id>/
├── trajectory.dcd            # MD 轨迹
├── complex_final.pdb         # 最终构象
├── energy_report.json        # ΔG（MM-GBSA）、能量分解
├── rmsd_profile.csv          # 时间序列 RMSD
├── qc_flags.json             # md_reliable / warnings / ml_charge 标签
└── summary.md                # 可读摘要
```

每个 run 末尾生成一份**总汇总表**：

```
deliveries/<run_id>/FINAL_RANKING.csv
# case_id, upstream_module, target, dG, rmsd_mean, qc_flag, final_rank
```

---

## 4) 跨模块交付契约：硬规定

1. **上游不满足契约 = 直接拒收**。不要 Simulation Hub 来擦屁股。
2. 每个上游模块在自己的 README 里必须引用 [`../ARCHITECTURE_AND_DATAFLOW.md#3-simulation-hub-统一输入契约`](../ARCHITECTURE_AND_DATAFLOW.md#3-simulation-hub-统一输入契约)。
3. 拒收要走正式渠道：`rejections/<case_id>/reason.md` 反馈给上游，而不是私下讨论。

---

## 5) 错误处理（已知痛点 → 已有方案）

| 场景                         | 现状                  | 计划                                             |
| ---------------------------- | --------------------- | ------------------------------------------------ |
| 上游给 SMILES 不给 3D        | 拒收                  | 统一契约规定必须 SDF                             |
| 上游 PDB 不干净（杂水/离子） | 人工清理很痛          | 卡在契约层拒收                                   |
| AM1-BCC 失败                 | 当前直接 raise        | 切 Espaloma + `[WARNING: ML_CHARGE]`（第二阶段） |
| MD 发散（RMSD 失控）         | md_reliable flag 标记 | 自动重跑 / 降级到更短 MD 窗口                    |
| OpenMM CUDA 版本/依赖地狱    | 手动锁版本            | 写死 Dockerfile + conda env.yml                  |

---

## 6) 第一阶段任务清单（Lead 自己负责）

1. **把 T2Lead Stage4 抽离**成独立仓库 `simulation-hub/`
    - 清理掉 T2Lead 特有的 disease pipeline 耦合
    - 暴露标准 CLI：`simhub run --input deliveries/<run_id>/to_simhub/<case_id>/`
2. **统一 `molecule_type` 分派逻辑**
    - 小分子：现成 Stage4
    - 蛋白-蛋白（peptide-MHC / antibody-antigen）：扩展一个蛋白-蛋白 MD 子流程
3. **锁定版本 + Docker**
    - OpenMM 版本（当前用的）、AmberTools、OpenFF Toolkit 全部写进 Dockerfile
4. **端到端 smoke test**
    - 4 下游各给 1 个 case_id，Simulation Hub 全部跑通
5. **总汇总报告**（FINAL_RANKING.csv + PDF/HTML dashboard）

---

## 7) 两阶段升级路线

### 第一阶段（锁定现状 + 严格契约，Must Do）

- 路线：`OpenFF + PDBFixer + Modeller`
- 重点：**严格执行交付契约** + 端到端跑通
- 不做重构，不做复杂替换

### 第二阶段（灰度升级，Should Do）

- **基于第一阶段真实失败案例**再升级：
    - `sqm` 失败日志积累 → Espaloma Fallback + `[WARNING: ML_CHARGE]`
    - PDBFixer 失败案例 → PDB2PQR + PROPKA 替换加氢
    - 长 MD 性能 → 多 GPU 调度正式化
- **不做的事**（第一阶段不碰）：
    - 不换 OpenMM `addSolvent`（够用）
    - 不上 PACKMOL（没到跨膜场景）
    - 不上 FEP/TI（MM-GBSA 先稳住）

---

## 8) 未解决的讨论题（Lead 自己要决定）

1. **数据流向**：上下游用文件交接（当前方案）还是 Protocol Buffers？  
    → 一阶段坚持文件交接，简单 > 完美。二阶段若模块间调用频繁再考虑。
2. **交付给广科医的形式**：HTML dashboard？PDF 报告？或者只给 CSV？  
    → 建议 HTML dashboard（最直观），底层数据 CSV 附带。
3. **开发标准**：OpenMM 版本锁定、Python 版本、最小单元测试集？  
    → 写进 `simulation-hub/CONTRIBUTING.md`。
4. **向 BGI / 广科医的沟通节奏**：每月一次 / 里程碑驱动？  
    → Lead 个人决定，不在本文档讨论。

---

## 9) 风险与边界

- MD / MM-GBSA 是近似方法，绝对数值不是真值。**Simulation Hub 提供的是相对排序 + 稳定性筛子**，不是"证明有效"的金标准。
- 加了 guardrail（md_reliable）仍是经验阈值，不是物理真值。
- 不支持共价抑制剂（共价结合需特殊流程）、变构位点 MD（需更长时间尺度）——**明确写入风险章节**。