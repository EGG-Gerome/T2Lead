# T2Lead

[English (README.md)](README.md)

**主路径 — 体细胞变异 → 突变蛋白 → 对接 / MD：** 将肿瘤/正常外显子数据（或 VEP 注释 VCF）与小分子先导库配对，针对**患者特异性**蛋白结构筛选候选化合物。

**次路径 — 疾病 → 靶点 → 苗头 → 先导：** 从疾病名称出发，经 ChEMBL、ML 虚拟筛选、CReM / REINVENT4 到阶段四优化的经典模块化流水线。

```
肿瘤/正常 FASTQ 或 VCF ──► variant_analysis（VCF → 突变 → 结构）
                                    │
            Lead CSV（阶段三）─────► 阶段四：Vina + ADMET + MM-GBSA（OpenMM）

疾病名称 ──► 阶段一至三（OpenTargets、ChEMBL、MPO / RL）──► 阶段四
```

## 快速开始

```bash
cd T2Lead
make install          # conda 环境 + RDKit/OpenMM + REINVENT4 + .env（下载量较大）
make test             # pytest（无需 GPU）

# 次路径（疾病驱动）
make run DISEASE="breast cancer"

# 主路径（变异驱动阶段四，需先准备 final_lead_candidates.csv，见文档）
# 在 configs/default_config.yaml 中设置 variant_analysis.enabled: true 与 vcf_path 或 FASTQ 路径
python scripts/run_pipeline.py -c my_variant_config.yaml -v
```

详细安装：[docs/guide/quickstart_zh.md](docs/guide/quickstart_zh.md)。

## 输出（按疾病 / 运行）

默认 `pipeline.output_layout.use_stage_subdirs` 下，产物位于 `<pipeline.out_dir>/<disease_slug>/`：

| 目录 | 内容 |
|------|------|
| `shared/chembl_stage2/`（默认，见 `target_to_hit.shared_library_dir`） | 跨疾病共用的 ChEMBL 分子/IC50 爬取、`crawl_state.json`、`fp_cache/` |
| `stage1_targets/` | `ranked_targets.csv` |
| `stage2_hits/` | 该疾病的训练集、`model_cache/`、虚拟筛选与 `final_hit_candidates.csv` 等（爬取与指纹缓存默认在共享目录） |
| `stage3_leads/` | `final_lead_candidates.csv`、REINVENT4 文件 |
| `stage4_optimization/` | PDB/PDBQT、`docking_poses/`、`optimized_leads.csv`、MD 轨迹 |
| `logs/` | 完整日志与精简日志 |

用户输入目录约定见 [data/user_inputs/README.md](data/user_inputs/README.md)（如有）。

## 阶段四评分

- **ADMET 硬过滤**：hERG 心毒、Veber 规则不通过、高合成难度（SA）化合物在 MD 前**直接剔除**。
- **CYP 风险**：仅作软权重（`lead_optimization.scoring.w_cyp_soft`）。
- **综合排名**：AutoDock Vina（25%）+ 隐式溶剂 MM-GBSA 项（40%）+ 构象/轨迹 RMSD 稳定性（35%）+ 可选 CYP 软项 — 权重见 [configs/default_config.yaml](configs/default_config.yaml)。
- **集成 MM-GBSA**：短独立 MD 采样复合物能量，报告近似 ΔG 的 **均值 ± 标准差**（详见 [docs/guide/configuration_zh.md](docs/guide/configuration_zh.md)）。

*评分对计算假设进行排名，不代表实验效价。*

## nf-core / sarek（可选上游）

WES/WGS 变异检出见 [docs/guide/variant_pipeline_zh.md](docs/guide/variant_pipeline_zh.md)。Makefile 目标：`make install-sarek`、`make run-sarek INPUT=samplesheet.csv`。**不影响** `make test`（单 Makefile）。

## 容器运行环境

阶段一至四的药物发现流水线**无需容器**，conda 环境直接运行即可。容器仅在运行 nf-core/sarek（变异检测上游）时需要。

### Docker — 本地开发 / 有 Docker 权限的服务器

适用：个人笔记本、Docker Desktop、有 sudo 的云主机。

```bash
# 运行 T2Lead 自身镜像
docker compose build t2lead
docker compose run --rm t2lead python -m pytest /app/tests -q

# 运行 sarek
nextflow run nf-core/sarek -profile docker ...
```

GPU：取消 [docker-compose.yml](docker-compose.yml) 中 `deploy` 段的注释。默认镜像使用 CPU 版 PyTorch。

### Apptainer (Singularity) — AutoDL / HPC / 受限容器环境

适用：AutoDL 实例、高校 HPC 集群、无法启动 `dockerd` 的任何环境。这些环境通常没有 Docker 所需的内核网络权限（iptables / 虚拟网桥），Apptainer 无需守护进程和网络特权，是 HPC 场景的标准方案。

```bash
# 安装 Apptainer（Ubuntu 22.04）
wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer_1.4.5_amd64.deb
sudo apt install -y ./apptainer_1.4.5_amd64.deb

# 运行 sarek
nextflow run nf-core/sarek -profile singularity ...
```

> **提示**：nf-core 官方同时维护 `-profile docker` 和 `-profile singularity`，两者功能完全等价，仅运行时后端不同。不是每个人都有本地 GPU，使用 AutoDL 等云 GPU 服务器时请选 Singularity 路径。

## 文档索引

| 文档 | 说明 | 英文 |
|------|------|------|
| [docs/guide/quickstart_zh.md](docs/guide/quickstart_zh.md) | 安装、环境变量、首次运行 | [quickstart.md](docs/guide/quickstart.md) |
| [docs/guide/variant_pipeline_zh.md](docs/guide/variant_pipeline_zh.md) | VCF/FASTQ → 结构 → 阶段四 | [variant_pipeline.md](docs/guide/variant_pipeline.md) |
| [docs/guide/disease_pipeline_zh.md](docs/guide/disease_pipeline_zh.md) | 疾病名驱动阶段一至四 | [disease_pipeline.md](docs/guide/disease_pipeline.md) |
| [docs/guide/configuration_zh.md](docs/guide/configuration_zh.md) | YAML / `DP_*` 参考 | [configuration.md](docs/guide/configuration.md) |
| [docs/guide/output_reference_zh.md](docs/guide/output_reference_zh.md) | 输出文件逐项说明 | [output_reference.md](docs/guide/output_reference.md) |
| [docs/reproducibility/reproduction_steps_zh.md](docs/reproducibility/reproduction_steps_zh.md) | 复现检查清单 | [reproduction_steps.md](docs/reproducibility/reproduction_steps.md) |
| [docs/reproducibility/test_data_preparation_zh.md](docs/reproducibility/test_data_preparation_zh.md) | WES/WGS 测试数据准备 | [test_data_preparation.md](docs/reproducibility/test_data_preparation.md) |
| [docs/research/somatic_variant_pipeline_feasibility_zh.md](docs/research/somatic_variant_pipeline_feasibility_zh.md) | 体细胞变异路径可行性 | [somatic_variant_pipeline_feasibility.md](docs/research/somatic_variant_pipeline_feasibility.md) |
| [docs/postmortems/20260330_md_openmm_quantity_type_mismatch.md](docs/postmortems/20260330_md_openmm_quantity_type_mismatch.md) | 复盘：OpenMM Quantity 类型兼容性 | — |

## 项目结构（简）

```
T2Lead/
├── configs/default_config.yaml
├── src/drugpipe/          # pipeline、各阶段、variant_analysis、paths.py
├── scripts/run_pipeline.py
├── tests/
└── docs/guide/ … docs/reproducibility/
```

## 致谢与参考

| 阶段 | 参考项目 | 借鉴内容 |
|------|---------|---------|
| **阶段一** | [OpenTargets Platform](https://platform.opentargets.org/) | 疾病-靶点关联公共 GraphQL API |
| **阶段二** | [ChEMBL REST API](https://www.ebi.ac.uk/chembl/) | 公共化合物与生物活性数据库 |
| **阶段三** | [CReM](https://github.com/DrrDom/crem)（Polishchuk, 2020） | 上下文感知片段替换，类似物枚举 |
| **阶段三** | [REINVENT4](https://github.com/MolecularAI/REINVENT4)（MolecularAI，Apache 2.0） | 强化学习生成式分子设计 |
| **阶段四** | [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)（Scripps Research） | 分子对接引擎及 Python 绑定 |
| **阶段四** | [OpenMM](https://openmm.org/)（Stanford，MIT License） | GPU 加速 MD，MM-GBSA 结合自由能 |
| **全部** | [RDKit](https://www.rdkit.org/) | 骨架、聚类、过滤、QED、SA、描述符 |

## 贡献者

- **卢振涛**：主要开发者，负责项目整体架构与核心代码实现（实习期间于 hupper）
- **姜可盈**：开发者，负责项目整体架构与核心代码实现（实习期间于 hupper）

## 许可证

本仓库为 **hupper 内部商用项目**，默认采用 **All rights reserved**（保留所有权利）。

未经 hupper 书面授权，禁止复制、分发、修改或用于任何公开发布。
