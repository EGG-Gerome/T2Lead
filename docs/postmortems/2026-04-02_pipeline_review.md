# T2Lead Pipeline Review & Improvement Plan

> Generated: 2026-04-02
> Status: Action items for Codex 5.3 Extra High
> Priority: P0 = must fix before tonight's rerun; P1 = important; P2 = enhancement

---

## 1. Opt(struct) 分数为什么上市药物这么低？（已诊断，需修复）

### 现象
| Drug | Dock | MD ΔG | RMSD | Opt(struct) |
|---|---|---|---|---|
| GEFITINIB | -7.75 | 3565.0 | 11.1 | 0.984 |
| AFATINIB DIMALEATE | -7.84 | 3579.0 | 22.0 | 0.306 |
| NERATINIB | -7.31 | — | — | 0.213 |
| OSIMERTINIB | -7.14 | — | — | 0.183 |
| MOBOCERTINIB | -6.52 | — | — | 0.074 |
| INAVOLISIB (breast) | -3.82 | — | — | 0.242 |

### 根因分析

**Opt(struct) 只看三个维度**（权重在 `default_config.yaml` line 218-223）：
- `w_docking: 0.25` — 对接打分
- `w_md_energy: 0.40` — MM-GBSA ΔG
- `w_stability: 0.35` — MD RMSD

它**不看** pIC50、QED、MPO、SA、ADMET。这些在 `Fast score` 里。

问题 1: **只有 top-N（md_top_n=2）做了 MD**，剩余 3 个药 MD 列全空。`_composite_score` 用 min-max 归一化，只有 2 个值参与归一化 → 排名几乎全看这 2 个值之间的相对差异。

问题 2: **MD ΔG 数值本身不正确**。GEFITINIB 的 `md_binding_energy=3565` kcal/mol（正值！），这说明 MM-GBSA complex 能量爆了（通常应该是 -20 ~ -80 kcal/mol 范围内的负值）。这个异常大正值是 `complex_energy - receptor_energy - ligand_energy` 计算溢出导致的，所以 Opt(struct) 的绝对值没有物理意义。

问题 3: **min-max 归一化是批内相对排名**，不是绝对打分。只有 5 个药参与归一化，最好的=1，最差的=0。所以 `0.074` 不代表"差"，只代表"这批里最差"。

### 解决方案（P0）
- [ ] `benchmark.md_top_n` 改为 `5`（或 `all`），让**全部**上市药都跑 MD
- [ ] 在 `md_simulation.py` 的 `_gb_energy_complex` 中加"能量异常检测"：如果 `ΔG > 0` 或 `|ΔG| > 5000`，标记为 `NaN` 而非写入假值
- [ ] Opt(struct) 归一化改为"绝对区间映射"而非 min-max，或在 dashboard 明确标注"此分数为批内相对排名"
- [ ] Dashboard 补充文字说明：Opt(struct) 是 Dock + MD 结构维度的排名分，不等于药物优劣

---

## 2. 靶点匹配：上市药是否针对正确靶点？（已修复，需验证）

### 现状
- **Lung**: Stage 1 选中 EGFR (CHEMBL203)，benchmark 搜的是 EGFR 上市药（GEFITINIB、AFATINIB、OSIMERTINIB、MOBOCERTINIB、NERATINIB）✅
- **Breast**: Stage 1 选中 PIK3CA (CHEMBL4005)，benchmark 搜的是 PIK3CA 上市药（ALPELISIB、COPANLISIB、INAVOLISIB）✅

### 已实现的机制
`pipeline.py` → `_apply_target_specific_stage4_pdb()` 根据 `chembl_to_pdb` 映射表自动把 Stage 4 受体 PDB 切换到正确结构：
```yaml
auto_target_pdb:
  chembl_to_pdb:
    CHEMBL203: "1M17"   # EGFR
    CHEMBL4005: "4JPS"  # PIK3CA
```

### 遗留风险（P1）
- [ ] 如果 Stage 1 选了一个不在 `chembl_to_pdb` 映射里的靶点，回退到 `pdb_id: "4JPS"`（默认值），可能用错受体
- [ ] 映射表目前只有 2 个靶点，需要扩充或改为 API 自动查 PDB（UniProt → PDBe → 挑共结晶配体最好的）

---

## 3. 目录结构混乱：EGFR_L858R 出现在 breast_cancer 下（已确认）

### 现象
```
breast_cancer/stage4_optimization/
├── optimized_leads.csv          ← 常规 Stage1→4 路径 (PIK3CA)
├── benchmark_drugs.csv
├── EGFR_L858R/                  ← WES variant 路径
│   ├── optimized_leads.csv
│   ├── benchmark_drugs.csv
│   ├── dashboard.html
│   └── ...
└── PIK3CA_H1047R/               ← 另一个 WES variant
    └── stage2_hits/
```

### 根因
当使用 `--vcf-path + --auto-stage23` 运行 variant-driven 管线时，`resolve_variant_stage4_out()` 把输出放在 `stage4_optimization/<GENE>_<MUTATION>/` 下。这和常规 Stage1→4 流程的 `stage4_optimization/` 共用同一个父目录。

### 问题
1. 用户无法从目录名区分"常规靶点发现路径"和"WES 变异驱动路径"
2. Dashboard 生成时可能读错 `optimized_leads.csv`（根目录还是子目录的？）
3. 未来如果多个患者数据，无法区分来源

### 解决方案（P1）
- [ ] **方案 A：按运行路径分层**
  ```
  breast_cancer/
  ├── standard/           ← Stage1 → PIK3CA → Stage4 (标准路径)
  │   └── stage4_optimization/
  └── variants/           ← WES-driven (变异路径)
      ├── EGFR_L858R/
      └── PIK3CA_H1047R/
  ```
- [ ] **方案 B：按患者+日期隔离**
  ```
  breast_cancer/
  ├── runs/
  │   ├── 20260401_standard_PIK3CA/
  │   └── 20260402_WES_patient001/
  │       ├── EGFR_L858R/
  │       └── PIK3CA_H1047R/
  ```
- [ ] **方案 C（推荐）：run manifest + hash ID**
  每次运行生成一个 `run_manifest.json`，包含 patient_id（可选）、pipeline_mode（standard/variant）、target_list、timestamp。目录名用 `<date>_<mode>_<short_hash>`。

### 关于患者隔离
用户提到未来需要每位患者独立（因为要加代谢数据权重）。建议：
- 用 `patient_id` 或 `sample_id`（而非姓名）作为目录标识
- 在 config 中新增 `pipeline.sample_id` 字段
- 输出路径: `<disease>/<sample_id>/runs/<timestamp>/`

---

## 4. Stage 1 选靶与 Stage 4 受体 PDB 不一致（已修复）

### 原始问题
Stage 1 选了 EGFR，但 Stage 4 一直用 `pdb_id: "4JPS"`（PIK3CA 的结构）。

### 修复措施（已实施）
`pipeline.py` 新增 `_apply_target_specific_stage4_pdb()`，在进入 Stage 4 之前根据 `target_chembl_id` 自动查表覆盖 `pdb_id`。

`default_config.yaml` 新增：
```yaml
lead_optimization:
  auto_target_pdb:
    enabled: true
    override_config_pdb: true
    chembl_to_pdb:
      CHEMBL203: "1M17"
      CHEMBL4005: "4JPS"
```

### 验证
Lung 的 dashboard 现在显示 `PDB: 1M17`，而不是 `4JPS`。✅

### 遗留（P1）
- [ ] `chembl_to_pdb` 映射需要扩展（目前只有 2 个）
- [ ] 长远应自动通过 ChEMBL → UniProt → PDBe 获取最佳共结晶结构

---

## 5. MD 环境依赖问题（已修复）

### 原始问题
当前 `t2lead` conda 环境中 OpenMM/openmmforcefields/openff-toolkit/mdtraj 虽然安装了，但：
1. 直接调用 `t2lead/bin/python` 时 `PATH` 缺少 `t2lead/bin`，导致 `antechamber/sqm` 不在搜索路径 → AM1-BCC 充电失败
2. `conda run -n t2lead bash -lc ...` 有时解析成 base env 的 python

### 修复
调用时显式设置 `PATH="/root/miniconda3/envs/t2lead/bin:$PATH"`。

### 遗留（P0）
- [ ] 在 pipeline 启动脚本（或 Makefile / shell wrapper）中固定 `PATH`
- [ ] 或改为 `conda activate t2lead` 后再跑（推荐写一个 `run.sh`）
- [ ] `sqm` 对某些大分子（如 AFATINIB DIMALEATE 的盐型）会崩溃 → 已在 `scorer.py` 中加了 `_largest_fragment_smiles()` 预处理

---

## 6. 上市药分数"低"的根本解释（如实回答）

### 这不全是 bug，但部分是 bug

**合理的部分：**
- 上市药是通过临床试验验证的，但它们的**计算预测分数**不一定高。计算对接/MD 是近似方法，受 force field 精度、采样不足、溶剂模型限制影响
- 一些上市药（如 OSIMERTINIB）是不可逆共价抑制剂，Vina 对接不能模拟共价键形成 → 评分偏低
- INAVOLISIB 的 Dock=-3.82 确实偏低，但它是 2024 年批准的高选择性 PI3Kα 抑制剂，可能需要特殊的结合模式（allosteric）

**有 bug 的部分：**
- MD ΔG 正值（+3565 kcal/mol）是 MM-GBSA 计算溢出，应该被标记为 NaN
- 只跑了 2/5 个药的 MD → 归一化失真
- 盐型 SMILES 没预处理 → AM1-BCC 崩溃 → MD 失败（已修复）
- benchmark docking `exhaustiveness=8` 偏低（Stage 4 正式跑用的是 32）

### 解决方案（P0）
- [ ] MD ΔG 正值自动标 NaN + 报警
- [ ] 全部 benchmark 药都跑 MD（md_top_n = all）
- [ ] benchmark docking exhaustiveness 与 Stage 4 保持一致（32）
- [ ] 共价抑制剂标记（在 drug_finder 中标注 `mechanism_type: covalent`），dashboard 显示提示

---

## 7. Dashboard 显示问题

### 7a. 颜色规则文字太小（P1）
当前 `.rules` CSS `font-size: .78rem`，在大屏上可以但对部分用户偏小。
- [ ] 改为 `.88rem`，加粗关键词
- [ ] 或改为带色块图例的行内组件，而非纯文本

### 7b. Chart 纵坐标挤在一起（已修复，需验证）
已实现 `bounds()` 函数做动态 Y 轴范围。
- [ ] 验证今晚重跑后 chart 轴是否合理

### 7c. Composite score 糊在一起（已修复）
已拆分为 Fast score（pIC50/QED/MPO/Dock/SA/ADMET）和 Opt(struct)（Dock+MD+RMSD），分别展示。

---

## 8. 结果归档与复用（P2）

### 用户需求
跑完的结果能保存下来，后续直接查看，按化学式或靶点索引。

### 解决方案
- [ ] 每次 pipeline 完成后在输出目录写一个 `run_manifest.json`：
  ```json
  {
    "run_id": "20260402_lung_EGFR_abc123",
    "disease": "lung cancer",
    "target": "EGFR (CHEMBL203)",
    "pdb_id": "1M17",
    "timestamp": "2026-04-02T15:52:00Z",
    "n_leads": 10,
    "n_benchmark": 5,
    "dashboard": "dashboard.html",
    "pipeline_mode": "standard"
  }
  ```
- [ ] 在 `<out_dir>/` 下维护一个 `run_index.json`，列出所有历史运行
- [ ] Dashboard 直接是静态 HTML，拷到任何地方都能打开

---

## 9. 今晚重跑前需完成的代码改动清单（P0）

| # | 改动 | 文件 | 说明 |
|---|---|---|---|
| 1 | MD ΔG 正值 → NaN | `md_simulation.py` | `if dG > 0 or abs(dG) > 5000: dG = NaN` |
| 2 | benchmark md_top_n = all | `default_config.yaml` | `md_top_n: 5` (或 = max_drugs) |
| 3 | benchmark docking exhaustiveness = 32 | `default_config.yaml` | `docking_exhaustiveness: 32` |
| 4 | Dashboard 颜色规则字号 | `html_builder.py` CSS | `.rules` 字号 → `.88rem` |
| 5 | 启动脚本固定 PATH | 新建 `scripts/run.sh` | `export PATH=.../t2lead/bin:$PATH` |
| 6 | benchmark 药物全跑 MD | `scorer.py` | 确认 md_top_n 生效 |
| 7 | Opt(struct) 说明注释加大 | `html_builder.py` CSS | score_note 字号 → `.85rem` |

---

## 10. 更长期的改进路线（P2）

| # | 改进 | 说明 |
|---|---|---|
| 1 | PDB 自动查询 | ChEMBL target → UniProt → PDBe API → 选最佳共结晶结构 |
| 2 | 共价抑制剂识别 | ChEMBL mechanism_of_action 标注 → 用 CovDock 或 warhead-aware scoring |
| 3 | 患者级隔离 | `pipeline.sample_id` + 按患者分目录 |
| 4 | 代谢数据权重 | 接入个体化 PBPK 模型（未来） |
| 5 | run_manifest + 索引 | 每次运行自动归档，支持浏览历史 |
| 6 | Dashboard 3D 查看 | 嵌入 NGL/Mol* 查看对接 pose |
| 7 | Streamlit 交互版 | 在静态 HTML 之外提供可选的交互式面板 |

---

## 附录：对话中所有已解决问题检查表

| 问题 | 状态 | 备注 |
|---|---|---|
| ChEMBL ID 可点击跳转 | ✅ 已解决 | dashboard 中 chembl_url 生成 `<a>` 标签 |
| Benchmark 模块 | ✅ 已实现 | `src/drugpipe/benchmark/` |
| HTML 报告自动生成 | ✅ 已实现 | `src/drugpipe/report/` |
| 日志问题按严重程度排序 | ✅ 已实现 | Error > Warning > Info, 良性折叠 |
| 中英文切换 | ✅ 已实现 | `i18n.py` + dashboard 右上角按钮 |
| Dashboard 颜色风格统一 | ✅ 已实现 | teal/cyan 暗色主题 |
| 不要让用户输入 | ✅ 已实现 | 全自动搜索上市药 |
| Breast funnel 参数丢失 | ✅ 已修复 | `find_latest_full_log` 优先完整日志 |
| 图表 tooltips 显示药名 | ✅ 已修复 | hover 显示 name + 多维度分数 |
| Stage 1-4 布局 | ✅ 已恢复 | "Stage N — xxx" 格式 |
| 靶点-受体不匹配 (lung) | ✅ 已修复 | auto_target_pdb 机制 |
| Breast benchmark 选错靶点 | ✅ 已修复 | 解析日志选实际靶点 |
| Opt=1 全相同 | ✅ 已修复 | fallback 到 MPO/pIC50 spread |
| Fast vs Opt(struct) 混淆 | ✅ 已修复 | 拆成两列 + 说明文字 |
| 图表 Y 轴从 0 开始 | ✅ 已修复 | 动态 bounds() |
| 盐型 SMILES 导致 MD 崩溃 | ✅ 已修复 | _largest_fragment_smiles |
| OpenMM 环境缺失 | ✅ 已确认存在 | PATH 问题，非包缺失 |
| MD ΔG 正值异常 | 🔴 **未修复** | 需加异常检测 |
| benchmark 只跑 2/5 MD | 🔴 **未修复** | 需改 md_top_n |
| 颜色规则文字太小 | 🔴 **未修复** | 需改 CSS |
| 目录结构混乱 | 🟡 **设计中** | 方案见第 3 节 |
| 结果归档 | 🟡 **设计中** | run_manifest 方案见第 8 节 |
