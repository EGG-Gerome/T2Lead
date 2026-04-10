# T2Lead 周报（2026-04-07 ~ 2026-04-10，汇报版）

[T2Lead 项目仓库](https://github.com/EGG-Gerome/T2Lead)

> 技术细节版见 `weekly_progress_20260407_20260410.md`。

## 一、本周结论

本周的主线，是把 T2Lead 从“标准路径能稳定产出结果”推进到“主路径默认不覆盖、可接真实病人数据、结果交付更规范”。

本周已经取得的四个阶段性成果如下：

1. **Stage 4 可靠性继续增强**  
   完成 AM1-BCC 运行链路恢复、benchmark checkpoint 污染修复和 strict MD-reliable 策略收敛，Stage 4 的输出现在更可信、更容易解释。

2. **Breast `PIK3CA` 结果进一步稳定化**  
   当前 breast 标准路径结果已稳定维持 10 个优化先导、4 个 benchmark comparator，并完成 dashboard / payload 一致性校验，适合继续作为当前最强展示样例。

3. **主路径默认改成“病人级隔离输出”**  
   现在变异驱动路径不再和标准 disease 路径共用同一结果目录，默认会按 `sample_id + 时间戳` 生成唯一运行目录，避免覆盖已有结果。

4. **真实数据入口方向已明确**  
   真实数据第一阶段优先走 GDC 的 BRCA `VCF` 路线，而不是直接上 FASTQ；当前已完成 synthetic VCF smoke，下一步是获取真实 `mutect2` somatic VCF、补 VEP 注释并进入主路径正式验证。

---

## 二、本周主要成果

### 2.1 Stage 4 与 benchmark

- 修复了 AM1-BCC 运行时链路，避免 MD 因环境暴露不稳定而悄悄掉分子。  
- 修复 benchmark 与 lead 共用 `md_progress.csv` 导致的 checkpoint 复用污染。  
- 保持 `md_reliable` / strict gating 作为默认决策边界，让结构项更可信。  

### 2.2 结果交付与图像资产整理

- `dashboard.html`、`report_payload.json` 继续作为稳定交付物。  
- 2D 结构图不再散落在 `stage4_optimization/` 根目录，而是整理为：
  - `visual_assets/leads/...`
  - `visual_assets/benchmark/...`
- 同时补了导出与整理脚本，减少后续人工打包成本。  

### 2.3 主路径（variant / WES 路线）工程化推进

本周对主路径的推进，不只是“又写了一些代码”，而是把真实落地需要的几个关键问题先解决了：

- 默认 Sarek profile 改为 `singularity`，更适配 AutoDL / HPC；  
- 主路径默认输出改为病人级隔离目录，避免覆盖标准路径；  
- 新增主路径一键脚本，减少长命令和手工配置负担；  
- 准备了最小 VCF、tiny FASTQ samplesheet 和公开 WES 模板，便于逐层验证。  

### 2.4 主路径现在做到哪一步

目前最准确的说法是：

- **已经完成**：  
  - VCF 入口 smoke；  
  - mutant FASTA / mutant structure 生成；  
  - 输出目录唯一化；  
  - GDC BRCA / WXS / VCF 真实数据入口筛选。  

- **还没完成**：  
  - 用真实 BRCA / PIK3CA VCF 跑到 `stage2_hits/` / `stage3_leads/`；  
  - 用真实数据完整跑到 `optimized_leads.csv`。  

也就是说，主路径已经进入“可正式接真实数据验证”的阶段，但还不能写成“真实病人端到端闭环已完成”。

---

## 三、当前重点问题

1. **真实数据不是“看见就能直接用”**  
   GDC 当前可见的 BRCA `raw_somatic_mutation.vcf.gz` 不一定带 VEP 注释，因此真实主路径验证前还需要检查 header，并在必要时补做 VEP。

2. **主路径完整自动化仍依赖 `PIK3CA -> ChEMBL target` 成功映射**  
   这条在 `PIK3CA` 场景下问题不大，但对其他不易映射的基因，后续仍需设计 fallback。

3. **历史结果虽已开始归档，但仍需进一步标准化**  
   目前已经把历史 variant 遗留目录从 breast 标准路径根目录中迁出，但未来最好补 `run_manifest` / run index，让多轮结果更容易管理。

---

## 四、下周计划

1. 选取一个真实 BRCA `mutect2_raw_somatic_mutation.vcf.gz`，完成是否带 VEP 的检查。  
2. 若无 `CSQ`，补做 VEP 注释，并把它正式接入主路径。  
3. 完成 BRCA / `PIK3CA` 的真实主路径首轮验证，至少跑到 `stage2_hits/` / `stage3_leads/`。  
4. 继续把主路径结果整理为与标准路径相同等级的可汇报交付物。  

---

## 五、一句话版

本周 T2Lead 的进展重点，是把主路径从“概念上支持变异 / WES”推进到“默认不覆盖、默认适配 AutoDL、可以开始接真实 BRCA / `PIK3CA` 数据验证”；同时继续提升 Stage 4 和 benchmark 的可信度，并把图像与报告输出整理为更规范的交付结构。

