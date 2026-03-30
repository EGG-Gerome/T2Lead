# 体细胞变异检测 → 药物发现流水线：可行性与工具全景分析

**日期：** 2026-03-30  
**背景：** T2Lead 项目——将肿瘤/正常配对变异检测整合到现有四阶段药物发现流水线（详见周进展报告 2026-03-27 第 5.3 节）

> **英文版：** [somatic_variant_pipeline_feasibility.md](somatic_variant_pipeline_feasibility.md)

---

## 1. 流水线可行性评估

### 1.1 拟议流水线

```
肿瘤 FASTQ ──┐                                          ┌─→ 突变 PDB（ESMFold/Boltz-2）
             ├─→ BWA-MEM2 → SAMtools → GATK Mutect2 → VEP/SnpEff ──┤
正常 FASTQ ──┘       （比对）  （排序/去重）  （体细胞检测）  （注释）  └─→ 变异报告
                                                                               │
                                                                     T2Lead 阶段四
                                                                   （对接 + MD，基于
                                                                    突变蛋白结构）
```

### 1.2 结论：完全可行

每个步骤均采用生产级、维护活跃且互操作性已验证的工具。GATK Best Practices 体细胞 SNV/Indel 流程是各大癌症中心的事实标准。T2Lead 的创新贡献在于**下游连接**：将注释的错义突变输入结构预测，再进入对接/MD——这是大多数临床流水线尚未覆盖的环节。

**主要风险：** 变异到结构的桥接（VEP 注释 → 突变蛋白序列 → ESMFold → 对接）需要定制胶水代码，目前没有现成工具覆盖全链路。

---

## 2. 开源工具：版本、链接与推荐

### 2.1 序列比对 — BWA-MEM2

| 条目 | 详情 |
|------|------|
| **最新版本** | v2.3（2025-06-30） |
| **GitHub** | https://github.com/bwa-mem2/bwa-mem2 |
| **许可证** | MIT |
| **安装** | BioConda：`conda install -c bioconda bwa-mem2` |
| **参考基因组** | GRCh38（hg38）——使用带 ALT-aware contigs 的 GATK 资源包版本 |

**为何选 BWA-MEM2 而非原版 BWA-MEM：** 通过 SIMD 向量化（AVX-512/AVX2/SSE4.1）在 x86 平台实现 3–5× 加速，SAM 输出与原版完全兼容，下游工具无需变更。

**值得关注的替代方案：** `mm2-fast`（https://github.com/bwa-mem2/mm2-fast）——minimap2 + BWA-MEM2 的 SIMD 加速，短读取场景仍处于实验阶段。

### 2.2 BAM 处理 — SAMtools / Sambamba

| 条目 | 详情 |
|------|------|
| **SAMtools 最新版** | v1.21（2024-09） |
| **GitHub** | https://github.com/samtools/samtools |
| **许可证** | MIT/Expat |
| **关键操作** | `samtools sort`、`samtools markdup`（或 Picard/GATK MarkDuplicatesSpark） |

**实践建议：** 对于 T2Lead 的典型规模（5 GB + 5 GB FASTQ，≈30× WES 或 ~5× WGS），`samtools markdup` 比 Picard 更快。对于 30× 以上全基因组，可考虑 `sambamba markdup`（更好并行化）或 GATK Spark 版去重。

### 2.3 体细胞变异检测 — GATK Mutect2

| 条目 | 详情 |
|------|------|
| **最新版本** | GATK 4.6.2.0（2025-04-14） |
| **GitHub** | https://github.com/broadinstitute/gatk |
| **Docker 镜像** | `broadinstitute/gatk:4.6.2.0` |
| **许可证** | BSD-3-Clause |
| **必需资源** | `af-only-gnomad.hg38.vcf.gz`、`small_exac_common_3.hg38.vcf.gz`、正常样本 Panel（可选但推荐） |

**完整体细胞检测流程：**

```bash
# 1. 体细胞变异检测（肿瘤/正常配对模式）
gatk Mutect2 \
  -R hg38.fa \
  -I tumor.bam \
  -I normal.bam \
  -normal normal_sample_name \
  --germline-resource af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals pon.vcf.gz \
  --f1r2-tar-gz f1r2.tar.gz \
  -O unfiltered.vcf.gz

# 2. 学习读取方向模型（交联伪影过滤）
gatk LearnReadOrientationModel -I f1r2.tar.gz -O artifact-priors.tar.gz

# 3. 估算污染
gatk GetPileupSummaries -I tumor.bam -V small_exac_common.vcf.gz \
  -L small_exac_common.vcf.gz -O tumor_pileups.table
gatk CalculateContamination -I tumor_pileups.table -O contamination.table

# 4. 过滤
gatk FilterMutectCalls \
  -R hg38.fa \
  -V unfiltered.vcf.gz \
  --contamination-table contamination.table \
  --ob-priors artifact-priors.tar.gz \
  -O filtered.vcf.gz
```

### 2.4 变异注释 — VEP 与 SnpEff

#### Ensembl VEP

| 条目 | 详情 |
|------|------|
| **最新版本** | release/115.2（2025-09-26） |
| **GitHub** | https://github.com/Ensembl/ensembl-vep |
| **Docker 镜像** | `ensemblorg/ensembl-vep:release_115.2` |
| **许可证** | Apache 2.0 |
| **关键插件** | CADD、REVEL、SpliceAI、AlphaMissense、ClinVar（115.0 起含体细胞分类） |

**T2Lead 集成关键：** VEP 的 `--hgvsp` 输出给出蛋白质层面的变异后果（如 `ENSP00000275493.2:p.Leu858Arg`），可解析为突变氨基酸序列用于 ESMFold 输入。

#### SnpEff / SnpSift

| 条目 | 详情 |
|------|------|
| **最新版本** | v5.4c（2026-02-23） |
| **官网** | https://pcingola.github.io/SnpEff |
| **许可证** | MIT |
| **依赖** | Java 21+ |

**T2Lead 的 VEP vs SnpEff 选择：**

| 评判标准 | VEP | SnpEff |
|---------|-----|--------|
| 蛋白质层面 HGVS | 详细，含转录本选择 | 较好，含 HGVS.p |
| 插件生态 | 更丰富（CADD、AlphaMissense、REVEL） | SnpSift 用于数据库过滤 |
| 速度（WES VCF） | ~5–15 分钟 | ~1–3 分钟 |
| 安装难度 | 较重（Perl + 缓存下载） | 较简单（单个 JAR） |
| **推荐** | **首选**——转录本详情更佳，更利于突变序列提取 | 次选 / 交叉验证 |

### 2.5 结构预测 — ESMFold 及替代方案

| 工具 | 版本 | GitHub | 核心特性 |
|------|------|--------|---------|
| **ESMFold** | esm v2.0.0 | https://github.com/facebookresearch/esm | 比 AF2 快 10–30×；无需 MSA；单序列输入，非常适合快速突变建模 |
| **ColabFold** | v1.6.1（2026-03） | https://github.com/sokrypton/ColabFold | AF2 + MMseqs2 MSA；远同源物精度高于 ESMFold |
| **OpenFold3** | v0.4.0（2026-03） | https://github.com/aqlaboratory/openfold-3 | 开源 AF3 再实现；支持蛋白+配体+核酸复合体；Apache 2.0 |
| **Boltz-2** | 见项目仓库 | MIT 许可声明 | 结构 + 亲和力；与 T2Lead hit-to-lead 阶段最相关 |

**针对突变结构预测的实践建议：**
- ESMFold 是务实之选——将突变全长序列（野生型 + VEP 点突变）输入，在单张 A100 上几秒内得到 PDB。
- AlphaFold2/ColabFold 精度更高，但需要 MSA 计算（每条序列约 10–30 分钟）。
- 对 T2Lead 的典型场景（已知癌症基因的单点错义突变，如 PIK3CA、EGFR、BRAF），RCSB 的野生型 PDB + FoldX/Rosetta 点突变建模的实际精度可能**高于** ESMFold 全从头预测，因为大多数常见肿瘤驱动基因已有晶体结构。

**T2Lead 实践推荐：**

```
若野生型 PDB 存在于 RCSB/AlphaFold DB：
    → FoldX BuildModel（快速点突变，秒级）
    → 或 Rosetta ddg_monomer
否则：
    → ESMFold 输入突变序列（当前 T2Lead 默认）
    → 若算力充足，ColabFold 精度更高
```

### 2.6 对接与 MD — 现有 T2Lead 阶段四工具 + 新增

| 工具 | 用途 | GitHub/URL |
|------|------|-----------|
| **AutoDock Vina** | 经典对接 | https://github.com/ccsb-scripps/AutoDock-Vina |
| **DiffDock** | 基于 ML 的盲对接（top-1 38% vs Vina 23%） | https://github.com/gcorso/DiffDock |
| **GROMACS** | MD 模拟 | https://github.com/gromacs/gromacs |
| **OpenMM** | MD 模拟（T2Lead 已集成） | https://github.com/openmm/openmm |
| **P2Rank** | 结合口袋预测 | https://github.com/rdk/p2rank |

### 2.7 集成流水线

#### nf-core/sarek（变异检测）

| 条目 | 详情 |
|------|------|
| **最新版本** | 3.8.1（2026-02-12） |
| **GitHub** | https://github.com/nf-core/sarek |
| **框架** | Nextflow ≥ 25.10.2 |
| **集成的体细胞检测器** | Mutect2、Strelka2、Manta、FreeBayes、TIDDIT、CNVkit |
| **注释** | VEP（115.0）、SnpEff、SnpSift |

Sarek 以单条 `nextflow run nf-core/sarek` 命令覆盖拟议流水线的 1–4 步（比对到注释），处理肿瘤/正常配对并产出带 VEP 注释的过滤 VCF。

**T2Lead 集成策略：** 使用 Sarek 处理上游变异检测，以其注释 VCF 输出作为与 T2Lead 结构预测和药物发现阶段的交接点。

#### 其他集成流水线

| 流水线 | 覆盖范围 | URL |
|--------|---------|-----|
| **GATK Best Practices（WDL）** | 比对 → Mutect2 → 过滤 | https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl |
| **bcbio-nextgen** | 完整变异检测（体细胞 + 胚系） | https://github.com/bcbio/bcbio-nextgen |
| **NVIDIA Clara Parabricks** | GPU 加速比对 + 变异检测 | https://www.nvidia.com/en-us/clara/parabricks/ |
| **OpenCRAVAT** | 变异注释 + 解读 | https://github.com/KarchinLab/open-cravat |

---

## 3. 计算资源需求与运行时间估算

### 3.1 数据规模背景

对于 5 GB + 5 GB FASTQ（肿瘤 + 正常）：
- 若为 **WES**（全外显子组测序）：约 100–150× 深度——覆盖度较高，常见于临床 Panel
- 若为 **WGS**：约 3–5× 深度——覆盖度较低，不适合体细胞检测（最低建议 30×）
- 最可能的场景：**WES ~100×** 或 **靶向 Panel**

### 3.2 运行时间估算（WES ~100×，基于 CPU）

| 步骤 | 工具 | CPU 核数 | 内存 | 耗时 | 存储 |
|------|------|---------|------|------|------|
| 序列比对 | BWA-MEM2 | 16 | 16 GB | ~30–45 分钟/样本 | 2× ~15 GB BAM |
| 排序 + 去重 | SAMtools | 8 | 8 GB | ~15–20 分钟/样本 | 同 BAM，重写 |
| Mutect2 | GATK 4.6.2 | 4–8 | 16 GB | ~2–4 小时 | VCF（~10 MB） |
| FilterMutectCalls | GATK | 1 | 4 GB | ~5 分钟 | 过滤后 VCF |
| 注释 | VEP | 1–4 | 8 GB | ~10–20 分钟 | 注释 VCF |
| **上游合计** | | **16 核** | **16 GB** | **~4–6 小时** | **~50 GB 临时空间** |

### 3.3 运行时间估算（WGS 30×，基于 CPU——参考）

| 步骤 | 工具 | 耗时 | 说明 |
|------|------|------|------|
| 序列比对 | BWA-MEM2 | ~2–4 小时/样本 | 推荐 32 核 |
| 排序 + 去重 | SAMtools | ~1–2 小时/样本 | |
| Mutect2 | GATK | ~12–24 小时 | 可用 `--intervals` 按染色体并行 |
| 注释 | VEP | ~30–60 分钟 | |
| **合计** | | **~20–36 小时** | 单节点，32 核 |

### 3.4 GPU 加速（NVIDIA Clara Parabricks）

在单台 g4dn.12xlarge（4× T4 GPU）上使用 Parabricks：
- WGS 30× 比对 + Mutect2：**~2–3 小时合计**（比 CPU 快约 10×）
- 费用：AWS spot 约 $15–20/样本对

### 3.5 下游（结构预测 + 对接）——T2Lead 已有估算

| 步骤 | 工具 | 耗时 | 硬件 |
|------|------|------|------|
| ESMFold（单蛋白） | ESM-2 | ~30 秒–2 分钟 | 1× A100/4090 |
| 对接（每个化合物） | Vina/DiffDock | ~1–5 分钟 | CPU/GPU |
| MD（每个化合物，10 ns） | OpenMM | ~1–4 小时 | GPU |

---

## 4. 商业竞争格局

### 4.1 临床基因组学 → 靶点发现平台

| 公司 | 平台 | 功能 | 与 T2Lead 的关系 |
|------|------|------|----------------|
| **Tempus** | Tempus Loop | AI 驱动的靶点发现（RWD + PDO + CRISPR 筛选）；与阿斯利康达成 $2 亿合作 | 将基因组学与药物靶点验证相连接，但不做基于结构的对接 |
| **Foundation Medicine**（罗氏） | FoundationOne CDx | 全面基因组检测（324 基因），报告可干预突变 + 匹配疗法 | 仅临床报告——止步于"这是突变和已获批药物" |
| **Guardant Health** | Guardant360 | 液体活检 ctDNA 变异检测 | 仅上游（变异检测），不涉及药物设计 |
| **Illumina** | DRAGEN | FPGA 加速比对 + 变异检测 | 商业最快的体细胞流水线，无药物发现集成 |

### 4.2 计算药物发现平台

| 公司 | 平台 | 方法 | 相关性 |
|------|------|------|--------|
| **Schrödinger** | Maestro/Glide/FEP+ | 基于物理的对接 + 自由能微扰；基于结构设计的黄金标准 | 可替换 T2Lead 阶段四的对接；商业许可约 $5 万+/年 |
| **Relay Therapeutics** | Dynamo® | 基于运动的药物设计；cryo-EM + 长时间尺度 MD + AI/ML | 其 PI3Kα 抑制剂 zovegalisib（RLY-2608）专为突变特异性癌症设计——与 T2Lead 目标高度相似 |
| **Recursion Pharmaceuticals** | Recursion OS | 表型筛选 + ML；PB 级生物数据集 | 更偏表型驱动而非基因型驱动 |
| **Insilico Medicine** | TargetPro + Chemistry42 | AI 靶点识别（71.6% 临床靶点检出率）+ 生成式化学 | 端到端从靶点到分子，但为专有系统 |
| **insitro** | Virtual Human™ | 多模态细胞数据上的 ML 靶点识别；与 BMS 合作 | 靶点识别为主，不涉及基于结构的设计 |

### 4.3 T2Lead 在竞争格局中的定位

**目前没有任何平台同时覆盖以下全部环节：**
1. 从原始 FASTQ 进行体细胞变异检测
2. 变异注释到蛋白质层面的后果
3. 突变结构预测
4. 基于对接/MD 的化合物筛选（针对突变结构）
5. 先导化合物优化（REINVENT4 生成式化学）

Tempus/Foundation Medicine 覆盖（1–2），Schrödinger/Relay 覆盖（3–4）。T2Lead 致力于以**开源端到端流水线**跨越全链路——这正是其独特价值主张。

---

## 5. 集成模式：变异检测 → 基于结构的药物发现

### 5.1 模式 A："VCF → 突变序列 → 全从头折叠 → 对接"

```
VCF（已过滤并注释）
  │
  ├── 解析 VEP/SnpEff：提取基因、转录本、HGVS.p
  │     例：PIK3CA p.His1047Arg
  │
  ├── 获取野生型蛋白序列（UniProt/Ensembl REST API）
  │
  ├── 应用氨基酸替换 → 突变序列
  │
  ├── ESMFold/ColabFold → 突变 PDB
  │
  └── 对接 + MD（T2Lead 阶段四）
```

**优点：** 全自动，无需人工干预。  
**缺点：** 对于研究透彻的癌基因，全从头预测结构精度可能低于实验结构 + 点突变建模。

### 5.2 模式 B："VCF → 已知结构 + FoldX 突变"（T2Lead 推荐）

```
VCF（已过滤并注释）
  │
  ├── 解析 VEP：提取基因、UniProt ID、残基变化
  │
  ├── 查询 RCSB PDB 获取已有晶体结构
  │     （或从 AlphaFold DB 获取预测结构）
  │
  ├── 若结构存在：
  │     → FoldX BuildModel 或 Rosetta relax + 点突变
  │     → 精度更高，保留配体结合背景
  │
  ├── 若无结构：
  │     → ESMFold 输入完整突变序列（当前降级方案）
  │
  └── 对接 + MD（T2Lead 阶段四）
```

**优点：** 利用实验数据；FoldX/Rosetta 点突变建模对稳定性预测已有充分验证。  
**缺点：** 需要检查 PDB 结构可用性；逻辑分支更复杂。

### 5.3 模式 C："集成方案"（研究级）

同时运行模式 A 和模式 B，对两个结构分别进行对接，采用共识打分。已在学术论文中应用（例：KRAS 别构稳定剂研究，AlphaFold + GROMACS + AutoDock）。

### 5.4 已发表集成案例

1. **KRAS 突变通用型药物发现**（2025）：AlphaFold 结构 → P2Rank 口袋检测 → AutoDock Vina → GROMACS 100 ns MD → ADMET 过滤。识别了 Switch-I/II 沟槽的别构稳定剂。

2. **Relay Therapeutics PI3Kα**（临床）：PIK3CA 突变体全长 cryo-EM → 长时间尺度 MD → 计算设计 → zovegalisib（RLY-2608），现已进入 PIK3CA 突变型癌症临床试验。

3. **AnadrosPilotMD 流水线**（GitHub：https://github.com/tdextermorse/anadrospilotmd）：P2Rank + AutoDock + AmberTools + GROMACS 集成工作流，支持并行执行。

---

## 6. T2Lead 实践建议

### 6.1 近期（下一个冲刺——2026-03-30 当周）

1. **使用 nf-core/sarek 3.8.1** 作为上游变异检测引擎，而非手动编写 BWA-MEM2 → SAMtools → Mutect2 脚本。优势：
   - 经过测试、可重现、容器化
   - 全面实现 GATK Mutect2 Best Practices 流程
   - 内置 VEP 115.0 注释
   - `--input samplesheet.csv` 直接接受肿瘤/正常配对

2. **编写 VCF 到突变序列的解析器**，作为 T2Lead 新模块（`src/drugpipe/variant_analysis/vcf_to_mutant_seq.py`）：
   - 输入：VEP 注释 VCF
   - 提取：驱动基因的错义变异（按 IMPACT=HIGH 或 MODERATE 过滤）
   - 输出：供 ESMFold 使用的突变 FASTA 序列

3. **与现有流水线的集成点：**
   - 新增阶段零或阶段一前置模块
   - 若提供了体细胞 VCF，绕过阶段一的靶点选择（突变本身即是靶点）
   - 将突变结构直接输入阶段四

### 6.2 中期（未来 2–4 周）

1. **实现模式 B**（已知结构 + FoldX 降级为 ESMFold）：
   - 通过 RCSB REST API 查询 PDB（`https://data.rcsb.org/rest/v1/core/uniprot/{uniprot_id}`）
   - 集成 FoldX 5（学术许可，非商业免费）用于已有结构的点突变建模

2. **突变结构质量基准测试：**
   - 对 PIK3CA H1047R、EGFR L858R、BRAF V600E 比较 ESMFold 突变体 vs FoldX 突变体 vs 实验突变体结构
   - 以 RMSD（对比实验结构）和对接姿势一致性作为评估指标

3. **增加 GPU 加速变异检测选项：**
   - 对有 GPU 资源的用户支持 NVIDIA Clara Parabricks
   - 10× 加速使交互式/批量处理成为可能

### 6.3 架构：拟议模块布局

```
src/drugpipe/
├── variant_analysis/           # 新增——体细胞变异检测桥接
│   ├── __init__.py
│   ├── sarek_runner.py         # 通过 Nextflow 启动 nf-core/sarek
│   ├── vcf_parser.py           # 解析 VEP 注释 VCF → 驱动突变
│   ├── mutant_sequence.py      # 生成突变蛋白序列
│   └── structure_bridge.py     # 路由到 FoldX（若有 PDB）或 ESMFold
├── target_identification/      # 现有阶段一
├── virtual_screening/          # 现有阶段二
├── hit_to_lead/                # 现有阶段三
└── lead_optimization/          # 现有阶段四
```

### 6.4 最低硬件需求

| 配置 | 适用场景 | 估算费用 |
|------|---------|---------|
| **16 核，32 GB RAM，500 GB SSD** | WES 变异检测 + 下游分析 | 自建或云端约 $2/小时 |
| **32 核，64 GB RAM，1 TB SSD，1× A100** | WGS + ESMFold + MD | 云端约 $5–8/小时 |
| **GPU 节点（4× T4/A10）** | Parabricks 加速 WGS | 云端约 $4–6/小时（竞价实例） |

---

## 7. 总结

| 维度 | 评估 |
|------|------|
| **技术可行性** | 高——所有工具均为生产级且有完善文档 |
| **创新贡献** | 变异→结构→对接的桥接是独特价值；尚无开源流水线覆盖全链路 |
| **最大风险** | 用于对接的突变结构精度——需对比已知突变晶体结构进行验证 |
| **推荐方案** | 上游使用 nf-core/sarek + 定制 VCF 解析器 + FoldX/ESMFold 混合结构预测 |
| **WES 算力需求** | 16 核机器约 4–6 小时；符合 T2Lead 现有工作节奏 |
| **商业差异化** | T2Lead 将成为首个连接肿瘤原始 FASTQ 与先导化合物的开源流水线 |
