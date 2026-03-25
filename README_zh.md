# T2Lead

[English (README.md)](README.md)

**Target → Hit → Lead → Optimized Lead**：端到端模块化药物发现流水线，输入疾病名称，输出经计算验证的优化候选先导化合物。

```
疾病名称 ──► 靶点发现 ──► 靶点到苗头 ──► 苗头到先导 ──► 先导优化 ──► 优化先导化合物 CSV
 (字符串)     (OpenTargets /   (ChEMBL 爬取 /     (骨架分析 /      (对接 / ADMET /
               OriGene)         ML + 虚筛 + ADMET)   CReM / REINVENT4) MD / MM-GBSA)
```

## 输入输出一览

| 阶段 | 输入 | 输出 |
|---|---|---|
| **1 — 靶点发现** | 疾病名称 `str`（如 `"肺癌"`）或 EFO ID（如 `EFO_0001378`） | 排序靶点列表：`chembl_id, symbol, name, score, source` |
| **2 — 靶点到苗头** | ChEMBL 靶点 ID `str`（如 `"CHEMBL4005"`） | `data/<疾病>/final_hit_candidates.csv` |
| **3 — 苗头到先导** | 阶段二输出的 Hit CSV（需含 `canonical_smiles`、`pred_pIC50_ens`） | `data/<疾病>/final_lead_candidates.csv` |
| **4 — 先导优化** | 阶段三输出的 Lead CSV + 靶蛋白 PDB ID | `data/<疾病>/optimized_leads.csv` |

## 重要说明：项目运行机制

本项目**不是预训练好直接推理的模型**，各阶段实际发生的事：

- **阶段一**：在线查询 OpenTargets 公共 API（秒级）。
- **阶段二**：首先爬取 ChEMBL 分子与 IC50 数据（视规模可能需要数小时），然后**现场训练**随机森林 + MLP 模型——没有预置模型，每次运行针对你的靶点从真实生物活性数据学习。
- **阶段三**：骨架/聚类分析（即时）；CReM 类似物生成（需片段库）；MPO 综合打分；可选 REINVENT4 强化学习分子生成。
- **阶段四**：准备靶蛋白结构（PDB），用 AutoDock Vina 对所有先导做分子对接，增强 ADMET 评估（SA 评分、hERG 心毒、CYP 抑制），可选 OpenMM MD 模拟 + MM-GBSA 结合自由能（GPU 加速）。

## 功能概览

- **阶段一 — 靶点发现**：按疾病名或 EFO ID 查询 OpenTargets GraphQL API；可选与 OriGene AI 代理合并结果。
- **阶段二 — 靶点到苗头**：爬取 ChEMBL 分子与 IC50 活性数据（支持断点续跑），构建 pIC50 数据集，训练随机森林 + Torch MLP 回归器（GPU 加速），全库虚拟筛选，ADMET / QED / PAINS 过滤。
- **阶段三 — 苗头到先导**：Murcko 骨架 SAR 分析、Butina 多样性聚类、CReM 类似物枚举、多参数优化（MPO）打分，可选 REINVENT4 强化学习分子生成（以效力预测器为奖励函数）。
- **阶段四 — 先导优化**：蛋白准备（RCSB PDB 下载 + PDBFixer 修复）、分子对接（AutoDock Vina）、增强 ADMET 评估（合成可及性评分、hERG 毒性基团、CYP 抑制风险、Veber 口服生物利用度规则）、可选 MD 模拟 + MM-GBSA（OpenMM + CUDA GPU 加速）、综合评分排序。
- **GPU 支持**：自动检测 CUDA（NVIDIA）、MPS（Apple Silicon）或回退到 CPU，通过 `pipeline.device` 配置。
- **可配置**：单一 YAML 文件控制全部参数；环境变量可覆盖任意字段。
- **模块化**：各阶段可独立运行，也可串联执行。

## 项目结构

```
T2Lead/
├── configs/default_config.yaml    # 全流水线参数
├── src/drugpipe/
│   ├── config.py                  # 配置加载
│   ├── pipeline.py                # 编排器与 CLI 入口
│   ├── target_discovery/          # 阶段一
│   │   ├── opentargets.py
│   │   ├── origene_client.py
│   │   └── target_ranker.py
│   ├── target_to_hit/             # 阶段二
│   │   ├── chembl_api.py
│   │   ├── dataset.py
│   │   ├── featurizer.py
│   │   ├── models.py
│   │   ├── screener.py
│   │   └── filters.py
│   ├── hit_to_lead/               # 阶段三
│   │   ├── scaffold.py
│   │   ├── clustering.py
│   │   ├── analog_gen.py
│   │   ├── mpo.py
│   │   ├── reinvent_bridge.py
│   │   └── lead_ranker.py
│   ├── lead_optimization/         # 阶段四
│   │   ├── protein_prep.py
│   │   ├── docking.py
│   │   ├── admet_deep.py
│   │   ├── md_simulation.py
│   │   └── lead_optimizer.py
│   └── utils/                     # 共享工具
│       ├── chem.py
│       ├── io.py
│       └── http.py
├── scripts/
│   ├── run_pipeline.py            # 运行完整流水线
│   └── run_stage.py               # 运行单阶段
└── data/                          # 输出目录（已加入 .gitignore）
    ├── logs/                      # 完整日志 + 精简日志
    ├── fp_cache/                  # 缓存的 Morgan 指纹
    └── <疾病>/                    # 各疾病输出子目录
```

## 安装

### 快速开始（Makefile）

```bash
# 核心安装
make install

# 可选：阶段三 REINVENT4
make install-reinvent4
make download-reinvent4-prior

# 运行全流程
make run DISEASE="breast cancer"
```

### 分步安装

### 环境要求

- Python >= 3.9
- RDKit（建议通过 conda 安装）
- 推荐 NVIDIA GPU（CUDA）用于 MLP 训练 + MD 模拟加速

> **GPU 兼容性**：RTX 50 系列（Blackwell，如 RTX 5090）需要 CUDA 12.8+ 的 PyTorch（`cu128` 版本）。RTX 40 系列（Ada）用 CUDA 12.4+（`cu124`）。RTX 30 系列（Ampere）用 CUDA 11.8+（`cu118`）。详见下方步骤 4。

### 安装步骤

```bash
# 1. 进入项目目录
cd T2Lead

# 2. 创建 conda 环境（推荐）
conda create -n t2lead python=3.11 -y
conda init && source ~/.bashrc	# 首次安装 conda 后要执行
conda activate t2lead
conda install -c conda-forge rdkit -y

# 3. 安装软件包与核心依赖（版本见 pyproject.toml）
pip install -e .

# 4. （可选）深度学习支持 — 根据你的 GPU 选择一行执行：
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128  # RTX 5090/5080 (Blackwell)
# pip install torch torchvision --index-url https://download.pytorch.org/whl/cu124  # RTX 4090/4080 (Ada)
# pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118  # RTX 3090/3080 (Ampere)
# pip install torch torchvision  # 仅 CPU（无 GPU）

# 5. （可选）CReM 用于阶段三类似物生成
pip install crem
bash scripts/download_crem_db.sh

# 6. （可选）REINVENT4 用于阶段三强化学习分子生成
git clone https://github.com/MolecularAI/REINVENT4.git
cd REINVENT4 && pip install -e . && cd ..
mkdir -p /root/REINVENT4/priors
curl -L "https://zenodo.org/api/records/15641297/files/reinvent.prior/content" \
  -o /root/REINVENT4/priors/reinvent.prior
# 可选：如果 REINVENT4 提示 prior 元数据 hash 无效
python scripts/fix_reinvent_prior_metadata.py /root/REINVENT4/priors/reinvent.prior

# 7. （可选）阶段四分子对接
pip install -e ".[docking]"   # vina、meeko、gemmi（见 pyproject.toml）

# 8. （可选）阶段四 MD 模拟（GPU 加速）
conda install -c conda-forge openmm pdbfixer mdtraj -y

# 9. （可选）MD 配体力场参数化（GAFF2）
conda install -c conda-forge openmmforcefields openff-toolkit -y
```

### 环境变量

```bash
cp .env.example .env
# 编辑 .env 填入所需密钥
# 推荐（AutoDL/容器环境根盘较小时）：把所有输出（数据集、缓存、日志、对接结果等）放到大盘挂载
export DP_PIPELINE__OUT_DIR=/autodl-fs/data/T2Lead
#
# （可选）复用已下载的 CReM 片段数据库（无需重新下载）
export DP_HIT_TO_LEAD__ANALOG_GEN__CREM_DB_PATH=/root/autodl-fs/crem_db/chembl33_sa25_f5.db
#
# Stage 3 使用 REINVENT4 时，必须使用生效的 DP_ 覆盖项：
export DP_HIT_TO_LEAD__REINVENT4__REINVENT_PATH=/root/miniconda3/envs/t2lead/bin/reinvent
export DP_HIT_TO_LEAD__REINVENT4__PRIOR_PATH=/root/REINVENT4/priors/reinvent.prior
```

## 使用说明

### 场景 A — 已知靶点（ChEMBL 中有数据，默认模式）

```bash
# 针对某疾病运行全部四阶段（自动发现靶点）
python scripts/run_pipeline.py --disease "breast cancer" -v

# 指定 ChEMBL 靶点 ID（跳过阶段一）
python scripts/run_pipeline.py --target CHEMBL4005 --stages target_to_hit hit_to_lead lead_optimization

# 自定义配置 + 详细日志
python scripts/run_pipeline.py -c my_config.yaml -v
```

在 `configs/default_config.yaml` 中设置蛋白结构用于对接：

```yaml
lead_optimization:
  pdb_id: "4JPS"    # RCSB PDB 实验晶体结构
```

### 单阶段运行

```bash
# 阶段一：靶点发现
python scripts/run_stage.py target_discovery --disease "breast cancer"

# 阶段二：靶点到苗头
python scripts/run_stage.py target_to_hit --target CHEMBL4005

# 阶段三：苗头到先导
python scripts/run_stage.py hit_to_lead

# 阶段四：先导优化（从 data/<疾病>/final_lead_candidates.csv 读取）
python scripts/run_stage.py lead_optimization
```

### 新靶点支持（不在 ChEMBL 中）

对于 ChEMBL 中数据不足或完全没有的新靶点，阶段二支持三种替代模式：

**场景 B — 用户自有 IC50 数据**（靶点有文献/BindingDB/自有实验数据，但不在 ChEMBL 中）：

```bash
# 提供自有 IC50 CSV（需包含列：molecule_chembl_id, target_chembl_id,
# standard_value, standard_type, standard_units）
python scripts/run_pipeline.py \
  --target MY_TARGET_ID \
  --activities-csv /path/to/my_ic50_data.csv \
  -v

# 还可以同时提供自有化合物筛选库
python scripts/run_pipeline.py \
  --target MY_TARGET_ID \
  --activities-csv /path/to/my_ic50_data.csv \
  --screening-library /path/to/my_compounds.csv \
  -v
```

**场景 C — 纯对接模式**（全新靶点，完全没有 IC50 活性数据）：

`--docking-only` 跳过 ML 流程（RandomForest + MLP 训练及虚拟筛选），因为没有 IC50 数据无法训练模型。如果你已经知道靶点，不需要 Stage 1（直接用 `--target`）。Stage 2 仅按类药性过滤候选分子。Stage 3 仍可通过 CReM 生成类似物。Stage 4 对接成为主要打分方式。

```bash
# 有已知 PDB 结构：
python scripts/run_pipeline.py \
  --target MY_NOVEL_TARGET \
  --docking-only \
  -v

# 只有氨基酸序列（无需 PDB —— ESMFold 自动预测 3D 结构）：
python scripts/run_pipeline.py \
  --target MY_NOVEL_TARGET \
  --docking-only \
  --protein-sequence "MTEYKLVVVGAVGVGKSALT..." \
  -v
```

Pipeline 按以下优先级自动获取蛋白结构：

| 优先级 | 来源 | 你需要提供的 |
|--------|------|-------------|
| 1 | RCSB PDB | 配置中设置 `pdb_id: "4JPS"` |
| 2 | 本地 PDB 文件 | 将 `.pdb` 文件放入输出目录 |
| 3 | **ESMFold API**（全自动） | `--protein-sequence "MTEYKLVV..."` 或在配置中设置 `protein_sequence` |

ESMFold（Meta AI）通过免费 REST API 从氨基酸序列预测接近 AlphaFold 质量的 3D 结构——无需注册、无需 GPU、已完全集成到流水线中。

> **序列长度限制**：ESMFold API 对 400 残基以下的蛋白效果最好。更长的蛋白建议从 [AlphaFold DB](https://alphafold.ebi.ac.uk/) 下载，或使用 [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)。

### 作为 Python 库调用

```python
from drugpipe.config import load_config
from drugpipe.pipeline import (
    run_target_discovery, run_target_to_hit,
    run_hit_to_lead, run_lead_optimization,
)

cfg = load_config(overrides={
    "target_discovery": {"disease": "breast cancer"},
    "pipeline": {"stages": ["target_discovery", "target_to_hit", "hit_to_lead", "lead_optimization"]},
    "lead_optimization": {"pdb_id": "4JPS"},
})

# 阶段一
targets = run_target_discovery(cfg)
target_id = targets[0]["chembl_id"]

# 阶段二
df_hits = run_target_to_hit(cfg, target_chembl_id=target_id)

# 阶段三
trainer = df_hits.attrs.get("_trainer")
featurizer = df_hits.attrs.get("_featurizer")
df_leads = run_hit_to_lead(
    cfg, df_hits,
    model_predict_fn=trainer.predict if trainer else None,
    featurizer_fn=featurizer.transform if featurizer else None,
)

# 阶段四
df_optimized = run_lead_optimization(cfg, df_leads)
```

## 配置说明

所有参数位于 `configs/default_config.yaml`：

| 配置节 | 说明 |
|---|---|
| `pipeline.stages` | 要运行的阶段（`target_discovery`、`target_to_hit`、`hit_to_lead`、`lead_optimization`） |
| `pipeline.seed` | 全局随机种子 |
| `pipeline.device` | 计算设备：`auto`、`cuda`、`mps`、`cpu` |
| `target_discovery.disease` | 输入疾病名或 EFO ID |
| `target_discovery.origene.enabled` | 是否调用 OriGene（需本地部署） |
| `target_to_hit.chembl.*` | ChEMBL 爬取参数 |
| `target_to_hit.model.*` | ML 模型超参 |
| `target_to_hit.filter.*` | ADMET 规则阈值 |
| `hit_to_lead.analog_gen.*` | CReM 类似物生成设置 |
| `hit_to_lead.mpo.*` | MPO 权重 |
| `hit_to_lead.reinvent4.*` | REINVENT4 设置（可执行路径、先验模型路径、步数） |
| `lead_optimization.pdb_id` | 靶蛋白 RCSB PDB ID |
| `lead_optimization.docking.*` | AutoDock Vina 对接参数 |
| `lead_optimization.admet_deep.*` | 增强 ADMET 评估设置 |
| `lead_optimization.md_simulation.*` | OpenMM MD 模拟设置 |
| `lead_optimization.scoring.*` | 综合评分权重 |

环境变量覆盖（`DP_` 前缀）：

```bash
export DP_TARGET_TO_HIT__CHEMBL__MAX_ACTIVITIES=50000
export DP_LEAD_OPTIMIZATION__PDB_ID=4JPS
```

## 各阶段详解

### 阶段一：靶点发现

**输入**：疾病名称字符串（如 `"breast cancer"`）或 EFO ID。

**过程**：
- 查询 OpenTargets Platform GraphQL API 获取疾病-靶点关联
- 将 Ensembl 基因 ID 映射为 ChEMBL 靶点 ID
- 可选合并 OriGene AI 推荐（需本地部署 [OriGene](https://github.com/GENTEL-lab/OriGene) 服务）

**输出**：带关联得分的 ChEMBL 靶点 ID 排序列表。

### 阶段二：靶点到苗头

**输入**：ChEMBL 靶点 ID 字符串。

**过程**：
1. **ChEMBL 爬取** — 下载分子与 IC50 活性数据（支持断点续跑）
2. **靶点选择** — 自动选择得分最高且 IC50 数据充足（>= 200 条）的靶点
3. **数据集构建** — 合并分子与活性，IC50 → pIC50，每分子取中位数聚合
4. **特征化** — 计算 2048 位 Morgan 指纹（ECFP4）
5. **训练** — 随机森林 + Torch MLP（GPU 加速）；集成预测
6. **虚拟筛选** — 对全库分子打分
7. **ADMET 过滤** — Lipinski、TPSA、QED >= 0.6、PAINS/Brenk、pIC50 >= 6.0

**输出**：`data/<疾病>/final_hit_candidates.csv`

### 阶段三：苗头到先导

**输入**：阶段二的 Hit CSV。

**过程**：
1. **骨架分析** — 提取 Murcko 骨架，按通用骨架分组，SAR 汇总
2. **多样性聚类** — 基于 Tanimoto 距离的 Butina 聚类
3. **类似物生成** — CReM 片段突变扩展化学空间
4. **MPO 打分** — 加权：40% 效力 + 25% QED + 20% ADMET + 15% 新颖性
5. **REINVENT4 RL** — 可选强化学习分子生成，以效力预测器作为奖励
6. **先导排序** — 最终排序输出 top-N

**输出**：`data/<疾病>/final_lead_candidates.csv`

### 阶段四：先导优化

**输入**：阶段三的 Lead CSV + 靶蛋白 PDB ID。

**过程**：
1. **蛋白准备** — 从 RCSB 下载 PDB，PDBFixer 修复（添加缺失原子、氢原子），自动检测结合位点，转 PDBQT
2. **分子对接** — AutoDock Vina 对接全部先导（CPU 多线程），计算结合亲和力（kcal/mol）
3. **增强 ADMET** — 合成可及性（SA）评分、hERG 心毒 SMARTS 检测、CYP3A4/2D6 抑制风险、Veber 口服生物利用度规则、综合风险评分
4. **MD 模拟** — （仅 top-N）OpenMM 隐式溶剂能量最小化 + MM-GBSA 结合自由能（CUDA GPU 加速）
5. **综合评分** — 对接分数 + ADMET 风险 + MD 结合能加权排序 → 最终候选

**输出**：`data/<疾病>/optimized_leads.csv`

## 输出文件

所有输出写入 `data/<疾病>/`（可通过 `pipeline.out_dir` 配置）：

| 文件 | 说明 | 产生于 |
|---|---|---|
| `fp_cache/morgan_*.npy` | 缓存的 Morgan 指纹（自动生成） | 阶段二 |
| `molecules_chemblid_smiles.csv` | 爬取的分子 | 阶段二 |
| `activities_ic50.csv` | 爬取的 IC50 活性数据 | 阶段二 |
| `dataset_target_ic50.csv` | 所选靶点的训练数据集 | 阶段二 |
| `scored_candidates.csv` | 带预测 pIC50 的全部分子 | 阶段二 |
| `final_hit_candidates.csv` | 过滤后的 hit 化合物 | 阶段二 |
| `h2l_scaffold_summary.csv` | 骨架 SAR 汇总 | 阶段三 |
| `h2l_cluster_summary.csv` | 多样性聚类汇总 | 阶段三 |
| `final_lead_candidates.csv` | 最终排序的候选先导化合物 | 阶段三 |
| `<PDB_ID>.pdb` | 下载的蛋白结构 | 阶段四 |
| `<PDB_ID>_receptor.pdbqt` | 准备好的受体（对接用） | 阶段四 |
| `docking_poses/` | 每个分子的最佳对接构象 | 阶段四 |
| `optimized_leads.csv` | 含全部评分的最终优化候选 | 阶段四 |

## 外部依赖

| 组件 | 是否必需 | 安装方式 | 用于 |
|---|---|---|---|
| RDKit | 是 | `conda install -c conda-forge rdkit` | 全部阶段 |
| PyTorch | 可选 | `pip install torch` | 阶段二（MLP）、阶段三（REINVENT4） |
| CReM | 可选 | `pip install crem` + [片段数据库](https://github.com/DrrDom/crem#databases) | 阶段三（类似物生成） |
| REINVENT4 | 可选 | `git clone` + `pip install -e .`，见 [REINVENT4](https://github.com/MolecularAI/REINVENT4) | 阶段三（RL 优化） |
| AutoDock Vina | 可选 | `pip install vina meeko gemmi` | 阶段四（分子对接） |
| OpenMM | 可选 | `conda install -c conda-forge openmm pdbfixer mdtraj` | 阶段四（MD 模拟） |
| openmmforcefields | 可选 | `conda install -c conda-forge openmmforcefields openff-toolkit` | 阶段四（MD 配体力场参数化） |
| OriGene | 可选 | 见 [OriGene 仓库](https://github.com/GENTEL-lab/OriGene) | 阶段一（AI 靶点推荐） |

所有可选依赖均优雅降级——未安装时跳过对应步骤并输出警告。

## 致谢与参考

| 阶段 | 参考项目 | 我们借鉴/学习了什么 |
|---|---|---|
| **阶段一** | [OriGene](https://github.com/GENTEL-lab/OriGene)（GENTEL-lab，CC-BY-NC-SA 4.0） | 多智能体 AI 治疗靶点发现方法；`origene_client.py` 提供对其服务的集成封装。 |
| **阶段一** | [OpenTargets Platform](https://platform.opentargets.org/) | 疾病-靶点关联、遗传证据、可成药性数据的公共 GraphQL API。 |
| **阶段二** | [ChEMBL REST API](https://www.ebi.ac.uk/chembl/) | EMBL-EBI 维护的公共化合物与生物活性数据库。 |
| **阶段三** | [CReM](https://github.com/DrrDom/crem)（Polishchuk, 2020） | 上下文感知片段替换，用于类似物枚举。 |
| **阶段三** | [REINVENT4](https://github.com/MolecularAI/REINVENT4)（MolecularAI，Apache 2.0） | 基于强化学习的生成式分子设计。 |
| **阶段四** | [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)（Scripps Research） | 分子对接引擎及 Python 绑定。 |
| **阶段四** | [OpenMM](https://openmm.org/)（Stanford，MIT License） | GPU 加速分子动力学工具包，用于 MM-GBSA 结合自由能计算。 |
| **全部** | [RDKit](https://www.rdkit.org/) | Murcko 骨架、Butina 聚类、PAINS/Brenk 过滤、QED、SA 评分、分子描述符。 |

## 贡献者

- **卢振涛**：主要开发者，负责项目整体架构与核心代码实现（实习期间于 hupper）
- **姜可盈**：开发者，负责项目整体架构与核心代码实现（实习期间于 hupper）

## 许可证

Apache-2.0 © 2026 hupper
