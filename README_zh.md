# T2Lead

[English (README.md)](README.md)

**Target → Hit → Lead**：端到端模块化药物发现流水线，输入疾病名称，输出排序后的候选先导化合物。

```
疾病名称 ──► 靶点发现 ──► 靶点到先导化合物 ──► 先导化合物优化 ──► 候选先导化合物 CSV
 (字符串)     (OpenTargets /   (ChEMBL 爬取 /     (骨架分析 /
               OriGene)         ML + 虚筛 + ADMET)   CReM / REINVENT4)
```

## 输入输出一览

| 阶段 | 输入 | 输出 |
|---|---|---|
| **靶点发现** | 疾病名称 `str`（如 `"肺癌"`）或 EFO ID（如 `EFO_0001378`） | 排序靶点列表：`chembl_id, symbol, name, score, source` |
| **靶点到先导化合物** | ChEMBL 靶点 ID `str`（如 `"CHEMBL204"`） | `data/final_hit_candidates.csv` — 列：`molecule_chembl_id, canonical_smiles, pred_pIC50_ens, pred_IC50_nM_ens, QED, MW, cLogP, TPSA, HBD, HBA, RotB, Rings, HeavyAtoms` |
| **先导化合物优化** | 阶段二输出的 Hit CSV（需含 `canonical_smiles`、`pred_pIC50_ens`） | `data/final_lead_candidates.csv` — 列：`canonical_smiles, source_smiles, origin, pred_pIC50_ens, QED, mpo_score, admet_pass, HasAlert, scaffold_smi, cluster_id` |

## 重要说明：项目运行机制

本项目**不是预训练好直接推理的模型**，各阶段实际发生的事：

- **阶段一**：在线查询 OpenTargets 公共 API（秒级）。
- **阶段二**：首先爬取 ChEMBL 分子与 IC50 数据（视 `MAX_MOLECULES` / `MAX_ACTIVITIES` 设置，可能需要数小时），然后在爬取的数据上**现场训练**随机森林 + MLP 模型——没有预置模型，每次运行针对你的靶点从真实生物活性数据学习。
- **阶段三**：骨架/聚类分析（即时完成）；若提供 CReM 片段数据库则生成类似物；用 MPO 综合打分。

## 功能概览

- **阶段一 — 靶点发现**：按疾病名或 EFO ID 查询 OpenTargets GraphQL API；可选与 OriGene AI 代理合并结果，实现 LLM 辅助靶点推荐。
- **阶段二 — 靶点到先导化合物**：爬取 ChEMBL 分子与 IC50 活性数据（支持断点续跑），构建 pIC50 数据集，训练随机森林 + Torch MLP 回归器，全库虚拟筛选，ADMET / QED / PAINS 过滤。
- **阶段三 — 先导化合物优化**：Murcko 骨架 SAR 分析、Butina 多样性聚类、CReM 类似物枚举、多参数优化（MPO）打分，可选 REINVENT4 强化学习分子生成。
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
│   └── utils/                     # 共享工具
│       ├── chem.py
│       ├── io.py
│       └── http.py
├── scripts/
│   ├── run_pipeline.py            # 运行完整流水线
│   └── run_stage.py               # 运行单阶段
├── data/                          # 输出目录（已加入 .gitignore）
└── notebooks/
    └── analysis.ipynb
```

## 安装

### 环境要求

- Python >= 3.9
- RDKit（建议通过 conda 安装）

### 安装步骤

```bash
# 1. 进入项目目录
cd T2Lead

# 2. 创建 conda 环境（推荐）
conda create -n t2lead python=3.11 -y
conda activate t2lead
conda install -c conda-forge rdkit -y

# 3. 安装核心依赖
pip install -r requirements.txt

# 4. 以可编辑模式安装
pip install -e .

# 5. （可选）深度学习支持
pip install torch

# 6. （可选）CReM 用于阶段三类似物生成
pip install crem
# 还需下载 CReM 片段数据库（约 2 GB）：
# https://github.com/DrrDom/crem#databases
```

### 环境变量

```bash
cp .env.example .env
# 编辑 .env 填入所需密钥
```

## 使用说明

### 完整流水线

```bash
# 针对某疾病运行全部三阶段
python scripts/run_pipeline.py --disease "lung cancer"

# 指定靶点（跳过阶段一）
python scripts/run_pipeline.py --target CHEMBL204 --stages target_to_hit hit_to_lead

# 自定义配置 + 详细日志
python scripts/run_pipeline.py -c my_config.yaml -v
```

### 单阶段运行

```bash
# 阶段一：靶点发现
python scripts/run_stage.py target_discovery --disease "breast cancer"

# 阶段二：靶点到先导化合物
python scripts/run_stage.py target_to_hit --target CHEMBL1862

# 阶段三：先导化合物优化
python scripts/run_stage.py hit_to_lead
python scripts/run_stage.py hit_to_lead --hits-csv data/final_hit_candidates.csv
```

### 作为 Python 库调用

```python
from drugpipe.config import load_config
from drugpipe.pipeline import run_target_discovery, run_target_to_hit, run_hit_to_lead

cfg = load_config(overrides={
    "target_discovery": {"disease": "liver cancer"},
    "pipeline": {"stages": ["target_discovery", "target_to_hit", "hit_to_lead"]},
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
```

## 配置说明

所有参数位于 `configs/default_config.yaml`：

| 配置节 | 说明 |
|---|---|
| `pipeline.stages` | 要运行的阶段 |
| `pipeline.seed` | 全局随机种子 |
| `target_discovery.disease` | 输入疾病名或 EFO ID |
| `target_discovery.origene.enabled` | 是否调用 OriGene |
| `target_to_hit.chembl.*` | ChEMBL 爬取参数 |
| `target_to_hit.model.*` | ML 模型超参 |
| `target_to_hit.filter.*` | ADMET 规则阈值 |
| `hit_to_lead.analog_gen.*` | CReM 类似物生成设置 |
| `hit_to_lead.mpo.*` | MPO 权重 |
| `hit_to_lead.reinvent4.*` | 可选 REINVENT4 |

环境变量覆盖（`DP_` 前缀）：

```bash
export DP_TARGET_TO_HIT__CHEMBL__MAX_ACTIVITIES=50000
export DP_HIT_TO_LEAD__MPO__W_POTENCY=0.5
```

## 各阶段详解

### 阶段一：靶点发现

**输入**：疾病名称字符串（如 `"lung cancer"`）或 EFO ID（如 `EFO_0001378`）。

**过程**：
- 查询 OpenTargets Platform GraphQL API 获取疾病-靶点关联
- 将 Ensembl 基因 ID 映射为 ChEMBL 靶点 ID
- 可选合并 OriGene AI 推荐

**输出**：带关联得分的 ChEMBL 靶点 ID 排序列表（控制台打印 + 传递给阶段二）。

**前置条件**：仅需联网，无需本地数据或模型。

### 阶段二：靶点到先导化合物

**输入**：ChEMBL 靶点 ID 字符串（来自阶段一或手动指定，如 `"CHEMBL204"`）。

**过程**：
1. **ChEMBL 爬取** — 通过 REST API 下载分子（SMILES）与 IC50 活性数据，支持断点续跑（视规模耗时数小时）
2. **数据集构建** — 合并分子与活性，IC50 转 pIC50，每分子中位数聚合
3. **特征化** — 计算 2048 位 Morgan 指纹（ECFP4）
4. **训练** — 在爬取数据上训练随机森林 + 可选 Torch MLP；集成预测
5. **虚拟筛选** — 用训练好的模型对全库分子打分
6. **ADMET 过滤** — Lipinski 规则、TPSA、QED >= 0.6、PAINS/Brenk 结构警示、pIC50 >= 6.0

**输出**：`data/final_hit_candidates.csv`

**前置条件**：联网访问 ChEMBL API。模型现场训练——无需预训练权重。

### 阶段三：先导化合物优化

**输入**：阶段二输出的 Hit CSV（或任何含 `canonical_smiles` 和 `pred_pIC50_ens` 列的 CSV）。

**过程**：
1. **骨架分析** — 提取 Murcko 骨架，按通用骨架分组，SAR 汇总
2. **多样性聚类** — 基于 Tanimoto 距离的 Butina 聚类
3. **类似物生成** — CReM 片段突变扩展化学空间（需 CReM + 片段数据库）
4. **MPO 打分** — 加权：40% 效力 + 25% QED + 20% ADMET + 15% 新颖性
5. **先导排序** — 最终排序输出 top-N

**输出**：`data/final_lead_candidates.csv`

**前置条件**：骨架/聚类/MPO 仅需 RDKit。CReM + 片段数据库为可选（类似物生成）。REINVENT4 为可选（RL 优化）。

## 输出文件

所有输出写入 `data/`（可通过 `pipeline.out_dir` 配置）：

| 文件 | 说明 | 产生于 |
|---|---|---|
| `molecules_chemblid_smiles.csv` | 爬取的分子 | 阶段二 |
| `activities_ic50.csv` | 爬取的 IC50 活性数据 | 阶段二 |
| `dataset_target_ic50.csv` | 所选靶点的训练数据集 | 阶段二 |
| `scored_candidates.csv` | 带预测 pIC50 的全部分子 | 阶段二 |
| `final_hit_candidates.csv` | 过滤后的 hit 化合物 | 阶段二 |
| `h2l_scaffold_summary.csv` | 骨架 SAR 汇总 | 阶段三 |
| `h2l_cluster_summary.csv` | 多样性聚类汇总 | 阶段三 |
| `final_lead_candidates.csv` | 最终排序的候选先导化合物 | 阶段三 |

## 外部依赖

| 组件 | 是否必需 | 安装方式 | 用于 |
|---|---|---|---|
| RDKit | 是 | `conda install -c conda-forge rdkit` | 全部阶段 |
| PyTorch | 可选 | `pip install torch` | 阶段二（MLP 模型） |
| CReM | 可选 | `pip install crem` + [片段数据库](https://github.com/DrrDom/crem#databases) | 阶段三（类似物生成） |
| REINVENT4 | 可选 | 见 [REINVENT4 仓库](https://github.com/MolecularAI/REINVENT4) | 阶段三（RL 优化） |
| OriGene | 可选 | 见 [OriGene 仓库](https://github.com/GENTEL-lab/OriGene) | 阶段一（AI 靶点推荐） |

## 致谢与参考

本项目在学习与研究以下项目的基础上构建而成，在此致以诚挚感谢：

| 阶段 | 参考项目 | 我们借鉴/学习了什么 |
|---|---|---|
| **阶段一** | [OriGene](https://github.com/GENTEL-lab/OriGene)（GENTEL-lab，CC-BY-NC-SA 4.0） | 多智能体 AI 治疗靶点发现方法；`origene_client.py` 提供对其服务的集成封装。将数据驱动（OpenTargets）与 AI 驱动（LLM）靶点推荐结合的思路来源于其架构。 |
| **阶段一** | [OpenTargets Platform](https://platform.opentargets.org/) | 疾病-靶点关联、遗传证据、可成药性数据的公共 GraphQL API。 |
| **阶段二** | 既有的 ChEMBL 工作流与教程 | 原始单文件 ChEMBL 工作流（API 爬取 → IC50 数据集 → Morgan 指纹 → RF/MLP → 虚拟筛选 → ADMET/QED 过滤）被重构为 6 个模块化类。核心算法逻辑（pIC50 转换、指纹特征化、RF+MLP 集成、Lipinski/PAINS 过滤）源于该脚本。 |
| **阶段二** | [ChEMBL REST API](https://www.ebi.ac.uk/chembl/) | EMBL-EBI 维护的公共化合物与生物活性数据库。 |
| **阶段三** | [CReM](https://github.com/DrrDom/crem)（Polishchuk, 2020） | 上下文感知片段替换，用于类似物枚举。 |
| **阶段三** | [REINVENT4](https://github.com/MolecularAI/REINVENT4)（MolecularAI，Apache 2.0） | 基于强化学习的生成式分子设计；我们的桥接模块生成 REINVENT4 配置并解析其输出。 |
| **阶段三** | [RDKit](https://www.rdkit.org/) | Murcko 骨架分解、Butina 聚类、PAINS/Brenk 过滤目录、QED、分子描述符。 |

## 贡献者

- **卢振涛**：主要开发者，负责项目整体架构与核心代码实现（实习期间于 hupper）
- **姜可盈**：开发者，负责项目整体架构与核心代码实现（实习期间于 hupper）

## 许可证

Apache-2.0 © 2026 hupper
