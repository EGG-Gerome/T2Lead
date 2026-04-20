# Simulation Hub 输入契约（SINGLE SOURCE OF TRUTH）

> **这是 SimHub 输入契约的唯一真相源**。所有下游模块（DrugReflector / RNAi-Therapy / PrecisionDelivery / ImmunoGen）在自己的 README / 任务书里**只引用此文档**，不要复制契约内容到其他地方 —— 防止契约漂移。
>
> 如需更新契约：**只改本文档**，然后通知所有下游模块负责人。

---

## 0) 铁律

1. **按 `molecule_type` 分派格式**。一个契约走天下会让力场 pipeline 崩（多肽 / 抗体塞 `.sdf` → AM1-BCC 量子化学计算直接 segfault 或算几小时垃圾数值）。
2. **SDF 只能装小分子**。它不保留氨基酸残基名 / 链 ID / SSBOND。把多肽或抗体塞进 `.sdf` 会让 AMBER 把它当"超大畸形小分子"用 GAFF 参数化 —— 物理上错、计算上炸。
3. **上游不满足契约 = SimHub 直接拒收**。SimHub 不擦屁股，只按 `rejections/<case_id>/reason.md` 返回错误码 + 修复建议。

---

## 1) 总览表

| molecule_type | 来源模块 | 文件契约 | 主力场 | 配体参数化 | 第一阶段 |
|--------------|---------|---------|-------|-----------|---------|
| `small_molecule` | DrugReflector / AnalogGen | `receptor.pdb` + `ligand.sdf` + `meta.json` | AMBER14SB（蛋白） + OpenFF/GAFF（配体） | **AM1-BCC** 量子化学 | ✅ 必做 |
| `peptide_mhc` | ImmunoGen | `complex.pdb` + `hla_allele.txt`（可选）+ `meta.json` | AMBER14SB（纯蛋白） | **直接查残基库** | ⚠️ ImmunoGen 入选则做 |
| `antibody_antigen` | PrecisionDelivery（ADC / SMDC / BiTE 抗体部分） | `complex.pdb` + `cdr_annotation.json`（可选）+ `meta.json` | AMBER14SB（纯蛋白） | **直接查残基库** | ⚠️ ADC 入选则做 |
| `sirna_ago2` | RNAi-Therapy | `complex.pdb` + `modifications.json` + `meta.json` | AMBER OL3（RNA）+ AMBER14SB（蛋白） | **直接查核苷酸/残基库**，化学修饰例外 | ❌ 默认不跑 MD |

> ADC 的 **payload + linker 另走 `small_molecule` 分支独立评估**，不要把整个抗体+linker+payload 塞进一次 MD（系统规模 150 倍差，尺度不匹配）。

---

## 2) 分支 A：`small_molecule`

**唯一合法使用 SDF 的分支**。

```
deliveries/<run_id>/to_simhub/<case_id>/
├── receptor.pdb      # 靶蛋白单体/结构域
│                     #   - 清理后：无结晶水、原配体、无关金属离子、硫酸根
│                     #   - 必要金属辅因子 (Zn/Mg) 可保留，但要在 meta.json 声明
├── ligand.sdf        # 3D 小分子
│                     #   - V2000 或 V3000，禁止 2D
│                     #   - 已对接到口袋（质心到受体最近原子 < 10 Å）
│                     #   - pH 7.4 加氢
│                     #   - 元素白名单：{C, H, O, N, S, P, F, Cl, Br, I}
└── meta.json         # 见 §6
```

**力场**：AMBER14SB（受体）+ OpenFF-2.x（首选）/ GAFF2（fallback）+ AM1-BCC 电荷。

---

## 3) 分支 B：`peptide_mhc`

```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb       # 纯蛋白多链
│                     #   Chain M  → MHC I α chain (~275 AA)
│                     #   Chain B  → β2-microglobulin (~99 AA)
│                     #   Chain P  → peptide (8-11 AA，通常 9-mer)
│                     #   全部标准氨基酸残基（ALA/GLY/...）
│                     #   禁止 UNK / HETATM 非标残基 / 糖基化（第一阶段）
├── hla_allele.txt    # 可选："HLA-A*02:01" 等
└── meta.json         # 见 §6
```

**力场**：AMBER14SB 纯蛋白；**绝对不走 GAFF / AM1-BCC**。

**校验规则**：
- 链 ID：至少 `M` + `P`；缺 `P` → `E_CHAIN_ID_MISSING`
- `P` 链残基数 ∈ [8, 11]
- 残基名必须在 AMBER14SB 标准 20 残基表内（含 HIS 三种质子化态 HID/HIE/HIP）

---

## 4) 分支 C：`antibody_antigen`

```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb       # 纯蛋白多链
│                     #   ADC / SMDC:  Chain H (VH) + Chain L (VL) + Chain A (抗原胞外域)
│                     #   BiTE:        额外 Chain H2 + Chain L2 + Chain T (CD3ε 胞外域)
│                     #   CDR 6 条 loop 不可缺失
│                     #   SSBOND 记录必须完整（抗体至少 4 对二硫键）
├── cdr_annotation.json    # 可选：CDR 残基编号 (Kabat / IMGT)
└── meta.json         # 见 §6
```

**力场**：AMBER14SB 纯蛋白；**绝对不走 GAFF / AM1-BCC**。

**校验规则**：
- 至少 3 条 chain ID（ADC）或 6 条（BiTE）
- `SSBOND` 记录数 ≥ 4；缺失 → `E_SSBOND_MISSING`
- CDR H3 残基不可缺失
- 抗原链（Chain A）必须是功能域，不能是全长膜蛋白（跨膜段走不了水盒子 MD）

---

## 5) 分支 D：`sirna_ago2`（第二阶段启用）

```
deliveries/<run_id>/to_simhub/<case_id>/
├── complex.pdb       # 多链
│                     #   Chain S → sense RNA 链
│                     #   Chain A → antisense RNA 链
│                     #   Chain P → Ago2 PIWI 域或全长 Ago2
│                     #   RNA 残基名：A / U / G / C（不是 DA/DT/DG/DC）
│                     #   化学修饰用 HETATM + CONECT 表示
├── modifications.json     # 每个修饰位置的化学结构 + SMILES（2'-OMe / PS / LNA 等）
└── meta.json         # 见 §6
```

**力场**：AMBER OL3（RNA）+ AMBER14SB（Ago2）。

**第一阶段默认不跑 MD**：RNAi 第一阶段交脱靶评分 + 二级结构即可，SimHub 这侧不接。

---

## 6) 所有分支共享的 `meta.json`

```json
{
  "case_id": "string (e.g., DR-20260420-0001)",
  "upstream_module": "DrugReflector|RNAi-Therapy|PrecisionDelivery|ImmunoGen|AnalogGen",
  "molecule_type": "small_molecule|peptide_mhc|antibody_antigen|sirna_ago2",
  "target_name": "string (e.g., KRAS_G12C, HLA-A*02:01)",
  "target_uniprot": "string | null",
  "upstream_score": "float (模块内部打分)",
  "priority": "high|medium|low",
  "notes": "free text",

  "parent_hit_id": "string | null",
  "scaffold_smiles": "string | null (仅 small_molecule 用)",
  "generation_source": "repositioning | analog_generation | de_novo"
}
```

`parent_hit_id / scaffold_smiles / generation_source` 三个字段仅小分子路径用于支撑 L1 ↔ L2 闭环；蛋白 / 核酸路径忽略即可。

---

## 7) 拒收原因码

| 代码 | 含义 |
|------|------|
| `E_MISSING_FILE` | 契约要求的文件缺失 |
| `E_WRONG_FORMAT` | 小分子传了 PDB，或多肽/抗体传了 SDF |
| `E_NONSTANDARD_RESIDUE` | `complex.pdb` 里出现 UNK / 非标残基但未在 `modifications.json` 声明 |
| `E_CHAIN_ID_MISSING` | 必需的 chain ID 不全 |
| `E_UNSUPPORTED_ELEMENT` | 小分子含 GAFF 不支持的元素 |
| `E_LIGAND_NOT_DOCKED` | `ligand.sdf` 远离 `receptor.pdb` 口袋（初始距离 > 10 Å） |
| `E_SSBOND_MISSING` | 抗体 PDB 没写 SSBOND |
| `E_META_INCOMPLETE` | `meta.json` 必填字段缺失 |
| `E_WRONG_RESIDUE_NAMING` | RNA 链用了 DNA 残基名（DA/DT/DG/DC） |

SimHub 拒收 → 写 `rejections/<case_id>/reason.md` → 上游修材料重新提交。

---

## 8) 版本历史

| 版本 | 日期 | 变更 | 影响 |
|------|------|------|------|
| v0.1 | 2026-04 | 初版：`protein.pdb + ligand.sdf + meta.json` 统一契约 | 【已废弃】多肽/抗体强行塞 SDF 会崩 |
| **v0.2** | **2026-04-17** | 按 `molecule_type` 分派 4 个分支契约；纯生物大分子路径改为 `complex.pdb` 方案 | 当前版本 |

**契约变更流程**：
1. 改本文档 + 版本号
2. 通知所有下游模块负责人
3. 给 1 周过渡期，下游更新 CI 自检脚本
4. SimHub 同步升级 dispatcher 实现
