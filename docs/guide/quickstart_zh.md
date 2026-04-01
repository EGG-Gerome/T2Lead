# 快速开始

[English (quickstart.md)](quickstart.md)

## 环境要求

- Linux 或 macOS（Windows 建议使用 WSL2）。
- **Conda**（Miniconda / Mambaforge），用于安装 RDKit 与 OpenMM。
- **Java 21+**，若通过 Nextflow 运行 nf-core/sarek。
- **容器运行时**（仅 sarek 变异检测需要，阶段一至四无需容器）：
  - 本地开发 / 有 Docker 权限的服务器 → **Docker**（`-profile docker`）
  - AutoDL / HPC / 无法启动 dockerd 的环境 → **Apptainer**（`-profile singularity`），见下方安装说明
- NVIDIA GPU 推荐用于 MLP 训练与 OpenMM 分子动力学。

### Apptainer 安装（AutoDL / HPC 用户）

AutoDL 等云 GPU 容器和高校 HPC 集群通常无法运行 Docker（缺少内核网络权限）。安装 Apptainer 后即可通过 Singularity 后端运行 sarek：

```bash
wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer_1.4.5_amd64.deb
sudo apt install -y ./apptainer_1.4.5_amd64.deb
apptainer --version   # 验证
```

## 安装

在仓库根目录执行：

```bash
make install
```

将创建 `t2lead` conda 环境，安装科学计算栈、以可编辑方式安装 `drugpipe`、PyTorch（默认：CUDA 12.4 轮子）、下载 CReM 数据库、REINVENT4 与 prior，并写入包含 `DP_*` 路径的 `.env`。

覆盖 PyTorch 索引，例如仅 CPU：

```bash
make install TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu
```

若不使用 Make，步骤与上述一致；可参考 `Makefile` 各目标或 git 历史中的 README。

## 环境变量

- `DP_PIPELINE__OUT_DIR` — 输出根目录（建议放在大容量磁盘上，例如 AutoDL：`/autodl-fs/data/T2Lead`）。
- `DP_*` — 嵌套覆盖配置；详见 [configuration_zh.md](configuration_zh.md)。

## 用户数据目录约定

（以下为建议路径，可在 YAML/CLI 中指向任意位置。）

- `data/user_inputs/activities/` — 用于自定义训练的 IC50 CSV。
- `data/user_inputs/screening_library/` — 待筛 SMILES 化合物库。
- `data/user_inputs/fastq/` — 肿瘤/正常样本 FASTQ，用于变异检出。
- `data/user_inputs/vcf/` — 带注释的 VCF，对应配置项 `variant_analysis.vcf_path`。

若仓库中有用户数据说明，参见项目内 `data/user_inputs/README.md`（如有）。

## 冒烟测试

```bash
make test
# 或
python -m pytest tests/ -v
```

## 下一步

- 以变异为中心的工作流：[variant_pipeline_zh.md](variant_pipeline_zh.md)
- 以疾病为中心的工作流：[disease_pipeline_zh.md](disease_pipeline_zh.md)
