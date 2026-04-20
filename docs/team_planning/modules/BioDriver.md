# BioDriver 模块任务书（生物驱动 / 上游中台）

**负责人：** 新人（计算生物学 / 生信工程方向）  
**模块定位：** 唯一对接 BGI 的上游中台；把原始组学数据转成 4 个下游模块可直接消费的"知识产品"。  
**一句话目标：** 给定一个病人的 BGI 数据包，输出 4 份标准化交付包，分别给 DrugReflector / RNAi / ADC / ImmunoGen。

---

## 0) 一分钟导读

> **这个模块在干什么？** 做病人数据的"第一工序"。拿 BGI 原始 WES + scRNA-seq，洗干净、算签名、找致病靶点、做变异注释 / HLA 分型 / 新抗原预测，最后按 4 条下游疗法线的"口味"打包分发。
>
> - **上游**：BGI（唯一对接点）
> - **下游**：DrugReflector（单细胞签名）· RNAi（致病靶点 + transcript）· ADC（表面靶点 + 胞外域结构）· ImmunoGen（HLA + 新抗原）
> - **第一阶段 smoke 目标**：跑通 1 个病人，只交付 2 份下游包（DrugReflector 必做 + 明天新人选定的另一条线）。
> - **关键风险**：数据质量依赖 BGI；BioDriver 自己要把 QC 关做扎实，否则下游全歪。

---

## 1) 输入数据（从 BGI 来）

**详细清单参见**：[`../BGI_DATA_REQUEST.md`](../BGI_DATA_REQUEST.md)

简要版：
- scRNA-seq：`.h5ad/.rds` + 细胞元数据 + 基因注释
- WES（肿瘤+正常配对 FASTQ）
- 样本临床元数据（脱敏）

---

## 2) 中间过程（工具链，允许自主选型）

### 2.1 体细胞变异检测（Somatic Variant Calling）
- 比对：BWA-MEM2
- 排序/去重：SAMtools / Picard
- 变异调用：GATK Mutect2（肿瘤+正常配对模式）
- 过滤：FilterMutectCalls + 污染估计
- 注释：VEP（首选，支持 `--hgvsp` 蛋白层级变异）或 SnpEff

### 2.2 HLA 分型
- OptiType（HLA-I，主流）/ HLA-HD / arcasHLA
- 输入：WES 的肿瘤或正常 BAM
- 输出：`HLA-A*02:01` 格式的 allele 列表

### 2.3 单细胞处理
- Scanpy 或 Seurat
- QC → 归一化 → 批次校正 → 聚类 → 细胞类型注释
- 肿瘤 vs 正常 DEA（Wilcoxon / MAST）
- 通路富集（gseapy / fgsea）作为合理性检查

### 2.4 靶点提名（Target Nomination）
- 对接 OpenTargets / DepMap / COSMIC 做证据整合
- （可选）Phosformer 做 PTM（磷酸化）风险评分
- 输出：一份**带分层的靶点清单**（按 Surface / Intracellular / Kinase / TF 分类），供 4 下游各取所需

### 2.5 结构准备（给 ADC 和 Simulation Hub 用）
- UniProt 查询 → AlphaFold DB 或 ESMFold 获取结构
- 对膜蛋白：提取胞外域（Ectodomain）
- 保持输出结构标准化（链 ID、残基编号与 UniProt 一致）

---

## 3) 输出数据（给 4 下游）

> **关键要求**：每类输出都要有 `version`, `run_id`, `patient_id`, `timestamp` 字段，并放到独立子目录。

### 3.1 → DrugReflector

```
deliveries/<run_id>/to_drugreflector/
├── preprocessed.h5ad          # QC/归一化后的 AnnData（细胞注释已加）
├── disease_signature.csv      # gene_symbol, v_score, lfc, adj_p
└── meta.json                  # 样本分组、DEA 方法、参数
```

### 3.2 → RNAi-Therapy

```
deliveries/<run_id>/to_rnai/
├── pathogenic_targets.csv     # gene_symbol, ensembl_gene_id, canonical_transcript_id, priority_score, reason
├── transcripts.fasta          # 上述 transcript 的 mRNA 序列
└── meta.json
```

### 3.3 → PrecisionDelivery (ADC)

```
deliveries/<run_id>/to_adc/
├── surface_targets.csv        # gene_symbol, uniprot_id, surface_confidence, tumor_expr, normal_expr, lfc
├── ectodomains/<uniprot>.pdb  # 每个靶点的胞外域结构
└── meta.json
```

### 3.4 → ImmunoGen (mRNA 疫苗)

```
deliveries/<run_id>/to_immunogen/
├── neoantigen_candidates.csv  # mutation, mut_peptide, wt_peptide, transcript_id, variant_vaf
├── hla_typing.json            # { "HLA-A": ["HLA-A*02:01", ...], "HLA-B": [...], "HLA-C": [...] }
└── meta.json
```

---

## 4) 第一阶段交付范围（重要！按下游实际选择来裁剪）

项目第一阶段**只启动 2 个下游模块**：`DrugReflector` + 明天新人选中的 1 个下游。  
BioDriver 第一阶段**不用一次性交付 4 份包**，按下游实际选择裁剪：

| 第一阶段下游 | BioDriver 必须做 | 可延后（第二阶段） |
|-------------|-----------------|-------------------|
| **DrugReflector**（已定） | 单细胞处理 + disease signature | — |
| **+ ImmunoGen** | 变异调用 + **HLA 分型 + Neoantigen 预测** | 表面靶点 + 胞外域 PDB / transcripts FASTA |
| **+ PrecisionDelivery** | 变异调用 + **表面靶点筛选 + 胞外域 PDB** | HLA / Neoantigen / transcripts FASTA |
| **+ RNAi-Therapy** | 变异调用 + **靶点 + transcripts FASTA** | HLA / Neoantigen / 胞外域 PDB |

> **结论**：第一阶段 BioDriver 基本就是"单细胞 + 变异调用 + 1 个下游特化输出"。三条特化线的复杂度：RNAi（轻）< ADC（中）< ImmunoGen（重）。

## 5) 建议执行步骤（第一阶段 smoke run）

1. 用**一个病人**的 BGI 数据跑通单细胞线 + 变异线。
2. 根据第一阶段选中的下游，跑对应特化线（HLA / 结构 / transcripts 三选一或两）。
3. 产出**对应的 1-2 份交付包**，和下游负责人对齐字段。
4. 发布 **v0.1 接口契约文档**（自己仓库 `INTERFACES.md`），下游按此实现。
5. 跑通后，把 smoke run 的交付包放到共享存储，供下游做 smoke test。

---

## 5) 交付物（你要给 Lead 看的）

- 独立 GitHub 仓库 `biodriver/`（README + 环境配置 + 运行命令）
- 一份 `INTERFACES.md`：**字段级**接口定义 + 示例文件
- 一份 `REPORT.md`：
  - 样本 / 数据来源
  - 工具链版本与参数
  - QC 结果（变异数、过滤前后、DEGs 数、HLA 分型置信度）
  - 4 份交付包的路径
- smoke run 的 4 份交付包放到共享存储（autodl-fs 或 S3）

---

## 6) 风险与边界

- **HLA 分型置信度**：低覆盖度 WES 可能给出低置信 allele；要在交付包里暴露置信度字段给 ImmunoGen。
- **单细胞分组**：肿瘤 vs 正常的"正常"是匹配正常组织还是只用非肿瘤细胞？要在 `meta.json` 写明策略。
- **不做的事**：不做药物设计、不做分子动力学、不做实验方案；只做"数据 → 知识产品"的转换。

---

## 7) 可探索方向（Roadmap，第二阶段）

- 多模态整合：把 WES 变异 + scRNA-seq 表达 + ATAC-seq 做联合证据评分
- AlphaGenome / Phosformer 深度整合做非编码区与 PTM 功能推断
- 自动化多病人批处理 + 全流程质控仪表盘
