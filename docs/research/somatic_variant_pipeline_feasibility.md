# Somatic Variant Detection → Drug Discovery Pipeline: Feasibility & Landscape Analysis

**Date:** 2026-03-30  
**Context:** T2Lead project — integrating tumor/normal paired variant calling into the existing 4-stage drug discovery pipeline (Section 5.3 of weekly progress 2026-03-27)

---

## 1. Pipeline Feasibility Assessment

### 1.1 Proposed Pipeline

```
Tumor FASTQ ──┐                                          ┌─→ Mutant PDB (ESMFold/Boltz-2)
              ├─→ BWA-MEM2 → SAMtools → GATK Mutect2 → VEP/SnpEff ──┤
Normal FASTQ ─┘         (align)   (sort/dedup)  (somatic call) (annotate)  └─→ Variant Report
                                                                                    │
                                                                          T2Lead Stage 4
                                                                        (Docking + MD on
                                                                         mutant structure)
```

### 1.2 Verdict: Fully Feasible

Each step uses production-grade, well-maintained tools with established interoperability. The GATK Best Practices workflow for somatic SNVs/indels is the de facto standard at major cancer centers. The novel contribution in T2Lead is the **downstream connection**: feeding annotated missense mutations into structure prediction and then into docking/MD — a step most clinical pipelines stop short of.

**Key risk:** The variant-to-structure bridge (VEP annotation → mutant protein sequence → ESMFold → docking) requires custom glue code. No off-the-shelf tool covers this span end-to-end.

---

## 2. Open Source Tools: Versions, URLs, and Recommendations

### 2.1 Alignment — BWA-MEM2

| Item | Detail |
|------|--------|
| **Latest version** | v2.3 (2025-06-30) |
| **GitHub** | https://github.com/bwa-mem2/bwa-mem2 |
| **License** | MIT |
| **Install** | BioConda: `conda install -c bioconda bwa-mem2` |
| **Reference genome** | GRCh38 (hg38) — use the GATK resource bundle version with ALT-aware contigs |

**Why BWA-MEM2 over original BWA-MEM:** 3–5× faster on x86 via SIMD vectorization (AVX-512/AVX2/SSE4.1). Produces identical SAM output, so all downstream tools are fully compatible.

**Alternative to watch:** `mm2-fast` (https://github.com/bwa-mem2/mm2-fast) — minimap2 with BWA-MEM2's SIMD acceleration, still experimental for short reads.

### 2.2 BAM Processing — SAMtools / Sambamba

| Item | Detail |
|------|--------|
| **SAMtools latest** | v1.21 (2024-09) |
| **GitHub** | https://github.com/samtools/samtools |
| **License** | MIT/Expat |
| **Key operations** | `samtools sort`, `samtools markdup` (or Picard/GATK MarkDuplicatesSpark) |

**Practical note:** For T2Lead's scale (5 GB + 5 GB FASTQ ≈ 30× WES or ~5× WGS), `samtools markdup` is faster than Picard. For full WGS at 30×+, consider `sambamba markdup` for parallelism or GATK's Spark-based dedup.

### 2.3 Somatic Variant Calling — GATK Mutect2

| Item | Detail |
|------|--------|
| **Latest version** | GATK 4.6.2.0 (2025-04-14) |
| **GitHub** | https://github.com/broadinstitute/gatk |
| **Docker** | `broadinstitute/gatk:4.6.2.0` |
| **License** | BSD-3-Clause |
| **Required resources** | `af-only-gnomad.hg38.vcf.gz`, `small_exac_common_3.hg38.vcf.gz`, panel of normals (optional but recommended) |

**Full somatic calling workflow:**

```bash
# 1. Call somatic variants (tumor/normal paired mode)
gatk Mutect2 \
  -R hg38.fa \
  -I tumor.bam \
  -I normal.bam \
  -normal normal_sample_name \
  --germline-resource af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals pon.vcf.gz \
  --f1r2-tar-gz f1r2.tar.gz \
  -O unfiltered.vcf.gz

# 2. Learn read orientation model (cross-linked artifact filter)
gatk LearnReadOrientationModel -I f1r2.tar.gz -O artifact-priors.tar.gz

# 3. Estimate contamination
gatk GetPileupSummaries -I tumor.bam -V small_exac_common.vcf.gz \
  -L small_exac_common.vcf.gz -O tumor_pileups.table
gatk CalculateContamination -I tumor_pileups.table -O contamination.table

# 4. Filter
gatk FilterMutectCalls \
  -R hg38.fa \
  -V unfiltered.vcf.gz \
  --contamination-table contamination.table \
  --ob-priors artifact-priors.tar.gz \
  -O filtered.vcf.gz
```

### 2.4 Variant Annotation — VEP & SnpEff

#### Ensembl VEP

| Item | Detail |
|------|--------|
| **Latest version** | release/115.2 (2025-09-26) |
| **GitHub** | https://github.com/Ensembl/ensembl-vep |
| **Docker** | `ensemblorg/ensembl-vep:release_115.2` |
| **License** | Apache 2.0 |
| **Key plugins** | CADD, REVEL, SpliceAI, AlphaMissense, ClinVar (somatic classifications as of 115.0) |

**Critical for T2Lead integration:** VEP's `--hgvsp` output gives protein-level consequence (e.g., `ENSP00000275493.2:p.Leu858Arg`) which can be parsed to generate the mutant amino acid sequence for ESMFold input.

#### SnpEff / SnpSift

| Item | Detail |
|------|--------|
| **Latest version** | v5.4c (2026-02-23) |
| **Website** | https://pcingola.github.io/SnpEff |
| **License** | MIT |
| **Requires** | Java 21+ |

**VEP vs SnpEff decision for T2Lead:**

| Criterion | VEP | SnpEff |
|-----------|-----|--------|
| Protein-level HGVS | Detailed, with transcript selection | Good, HGVS.p included |
| Plugin ecosystem | Richer (CADD, AlphaMissense, REVEL) | SnpSift for database filtering |
| Speed (WES VCF) | ~5–15 min | ~1–3 min |
| Ease of install | Heavier (Perl + cache download) | Simpler (single JAR) |
| **Recommendation** | **Primary** — better transcript-level detail for mutant sequence extraction | Secondary / validation |

### 2.5 Structure Prediction — ESMFold & Alternatives

| Tool | Version | GitHub | Key Trait |
|------|---------|--------|-----------|
| **ESMFold** | esm v2.0.0 | https://github.com/facebookresearch/esm | 10–30× faster than AF2; no MSA needed; single-sequence input ideal for quick mutant modeling |
| **ColabFold** | v1.6.1 (2026-03) | https://github.com/sokrypton/ColabFold | AF2 + MMseqs2 MSA; higher accuracy than ESMFold for distant homologs |
| **OpenFold3** | v0.4.0 (2026-03) | https://github.com/aqlaboratory/openfold-3 | Open AF3 reimplementation; protein+ligand+nucleic acid complexes; Apache 2.0 |
| **Boltz-2** | See project repo | MIT license claimed | Structure + affinity; most relevant for T2Lead's hit-to-lead stage |

**For mutant structure prediction specifically:**
- ESMFold is the pragmatic choice — feed it the mutant full-length sequence (wildtype + point substitution from VEP), get a PDB in seconds on a single A100.
- AlphaFold2/ColabFold produces higher-confidence structures but requires MSA computation (~10–30 min per sequence).
- For the T2Lead use case (single missense mutations in known cancer genes like PIK3CA, EGFR, BRAF), the wildtype PDB from RCSB + FoldX/Rosetta point-mutation modeling may actually be more accurate than full de novo prediction with ESMFold, since crystal structures exist for most common oncogene targets.

**Practical recommendation for T2Lead:**

```
if wildtype PDB exists in RCSB/AlphaFold DB:
    → FoldX BuildModel (fast point-mutation, ~seconds)
    → or Rosetta ddg_monomer
else:
    → ESMFold with mutant sequence (current T2Lead default)
    → ColabFold for higher accuracy if compute allows
```

### 2.6 Docking & MD — Existing T2Lead Stage 4 Tools + Additions

| Tool | Purpose | GitHub/URL |
|------|---------|------------|
| **AutoDock Vina** | Classical docking | https://github.com/ccsb-scripps/AutoDock-Vina |
| **DiffDock** | ML-based blind docking (38% top-1 vs 23% Vina) | https://github.com/gcorso/DiffDock |
| **GROMACS** | MD simulations | https://github.com/gromacs/gromacs |
| **OpenMM** | MD simulations (already in T2Lead) | https://github.com/openmm/openmm |
| **P2Rank** | Binding pocket prediction | https://github.com/rdk/p2rank |

### 2.7 Integrated Pipelines

#### nf-core/sarek (Variant Calling)

| Item | Detail |
|------|--------|
| **Latest version** | 3.8.1 (2026-02-12) |
| **GitHub** | https://github.com/nf-core/sarek |
| **Framework** | Nextflow ≥ 25.10.2 |
| **Somatic callers included** | Mutect2, Strelka2, Manta, FreeBayes, TIDDIT, CNVkit |
| **Annotation** | VEP (115.0), SnpEff, SnpSift |

Sarek covers steps 1–4 of the proposed pipeline (alignment through annotation) as a single `nextflow run nf-core/sarek` command. It handles tumor/normal pairing, produces filtered VCFs with VEP annotation.

**T2Lead integration strategy:** Use Sarek for the upstream variant calling, then parse its annotated VCF output as the handoff point into T2Lead's structure prediction and drug discovery stages.

#### Other Integrated Pipelines

| Pipeline | Scope | URL |
|----------|-------|-----|
| **GATK Best Practices (WDL)** | Alignment → Mutect2 → filtering | https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl |
| **bcbio-nextgen** | Full variant calling (somatic + germline) | https://github.com/bcbio/bcbio-nextgen |
| **NVIDIA Clara Parabricks** | GPU-accelerated alignment + variant calling | https://www.nvidia.com/en-us/clara/parabricks/ |
| **OpenCRAVAT** | Variant annotation + interpretation | https://github.com/KarchinLab/open-cravat |

---

## 3. Compute Requirements & Runtime Estimates

### 3.1 Data Scale Context

For 5 GB + 5 GB FASTQ (tumor + normal):
- If **WES** (whole exome): this is ~100–150× depth — very high coverage, typical for clinical panels
- If **WGS**: this is ~3–5× depth — low coverage, atypical for somatic calling (30× minimum recommended)
- Most likely interpretation: **WES at ~100×** or **targeted panel**

### 3.2 Runtime Estimates (WES ~100×, CPU-based)

| Step | Tool | CPU Cores | RAM | Wall Time | Storage |
|------|------|-----------|-----|-----------|---------|
| Alignment | BWA-MEM2 | 16 | 16 GB | ~30–45 min (per sample) | 2× ~15 GB BAM |
| Sort + Dedup | SAMtools | 8 | 8 GB | ~15–20 min (per sample) | Same BAMs, rewritten |
| Mutect2 | GATK 4.6.2 | 4–8 | 16 GB | ~2–4 hours | VCF (~10 MB) |
| FilterMutectCalls | GATK | 1 | 4 GB | ~5 min | Filtered VCF |
| Annotation | VEP | 1–4 | 8 GB | ~10–20 min | Annotated VCF |
| **Total upstream** | | **16 cores** | **16 GB** | **~4–6 hours** | **~50 GB scratch** |

### 3.3 Runtime Estimates (WGS 30×, CPU-based — for reference)

| Step | Tool | Wall Time | Notes |
|------|------|-----------|-------|
| Alignment | BWA-MEM2 | ~2–4 hours per sample | 32 cores recommended |
| Sort + Dedup | SAMtools | ~1–2 hours per sample | |
| Mutect2 | GATK | ~12–24 hours | Parallelizable by chromosome with `--intervals` |
| Annotation | VEP | ~30–60 min | |
| **Total** | | **~20–36 hours** | Single node, 32 cores |

### 3.4 GPU Acceleration (NVIDIA Clara Parabricks)

Parabricks on a single g4dn.12xlarge (4× T4 GPUs):
- WGS 30× alignment + Mutect2: **~2–3 hours total** (10× faster than CPU)
- Cost: ~$15–20 per sample pair on AWS spot

### 3.5 Downstream (Structure + Docking) — Already Known from T2Lead

| Step | Tool | Time | Hardware |
|------|------|------|----------|
| ESMFold (single protein) | ESM-2 | ~30s–2 min | 1× A100/4090 |
| Docking (per compound) | Vina/DiffDock | ~1–5 min | CPU/GPU |
| MD (per compound, 10 ns) | OpenMM | ~1–4 hours | GPU |

---

## 4. Commercial Landscape

### 4.1 Clinical Genomics → Target ID Platforms

| Company | Platform | What It Does | Relevance to T2Lead |
|---------|----------|-------------|---------------------|
| **Tempus** | Tempus Loop | AI-powered target discovery using RWD + PDOs + CRISPR screening; AstraZeneca partnership ($200M) | Connects genomic profiling to drug target validation; does not do structure-based docking |
| **Foundation Medicine** (Roche) | FoundationOne CDx | Comprehensive genomic profiling (324 genes); reports actionable mutations + matched therapies | Clinical reporting only — stops at "here's the mutation and approved drug" |
| **Guardant Health** | Guardant360 | Liquid biopsy ctDNA variant calling | Upstream only (variant detection), no drug design |
| **Illumina** | DRAGEN | FPGA-accelerated alignment + variant calling | Fastest somatic pipeline commercially; no drug discovery integration |

### 4.2 Computational Drug Discovery Platforms

| Company | Platform | Approach | Relevance |
|---------|----------|----------|-----------|
| **Schrödinger** | Maestro/Glide/FEP+ | Physics-based docking + free energy perturbation; gold standard for structure-based design | Could replace T2Lead Stage 4 docking; commercial license ~$50K+/year |
| **Relay Therapeutics** | Dynamo® | Motion-based drug design; cryo-EM + long-timescale MD + AI/ML | Their PI3Kα inhibitor zovegalisib (RLY-2608) was designed for mutation-specific cancers — directly analogous to T2Lead's goal |
| **Recursion Pharmaceuticals** | Recursion OS | Phenotypic screening + ML; petabyte-scale biological datasets | More phenotype-driven than genotype-driven |
| **Insilico Medicine** | TargetPro + Chemistry42 | AI target identification (71.6% clinical target retrieval) + generative chemistry | End-to-end from target to molecule, but proprietary |
| **insitro** | Virtual Human™ | ML on multi-modal cellular data for target ID; BMS partnership | Target identification focus, not structure-based |

### 4.3 Where T2Lead Fits in the Landscape

**No existing platform combines all of:**
1. Somatic variant calling from raw FASTQ
2. Variant annotation to protein-level consequences
3. Mutant structure prediction
4. Compound screening via docking/MD against the mutant structure
5. Lead optimization (REINVENT4 generative chemistry)

Tempus/Foundation Medicine cover (1–2). Schrödinger/Relay cover (3–4). T2Lead aims to bridge the entire span as an **open-source, end-to-end pipeline** — this is the unique value proposition.

---

## 5. Integration Patterns: Variant Calling → Structure-Based Drug Discovery

### 5.1 Pattern A: "VCF → Mutant Sequence → De Novo Fold → Dock"

```
VCF (filtered, annotated)
  │
  ├── Parse VEP/SnpEff: extract gene, transcript, HGVS.p
  │     e.g., PIK3CA p.His1047Arg
  │
  ├── Fetch wildtype protein sequence (UniProt/Ensembl REST API)
  │
  ├── Apply substitution → mutant sequence
  │
  ├── ESMFold/ColabFold → mutant PDB
  │
  └── Docking + MD (T2Lead Stage 4)
```

**Pros:** Fully automated, no manual intervention.  
**Cons:** De novo fold may be less accurate than experimental structure + point-mutation modeling for well-studied oncogenes.

### 5.2 Pattern B: "VCF → Known Structure + FoldX Mutation" (Recommended for T2Lead)

```
VCF (filtered, annotated)
  │
  ├── Parse VEP: extract gene, UniProt ID, residue change
  │
  ├── Query RCSB PDB for existing crystal structure
  │     (or AlphaFold DB for predicted structure)
  │
  ├── If structure exists:
  │     → FoldX BuildModel or Rosetta relax + point mutation
  │     → Higher accuracy, retains ligand binding context
  │
  ├── If no structure:
  │     → ESMFold with full mutant sequence (current fallback)
  │
  └── Docking + MD (T2Lead Stage 4)
```

**Pros:** Leverages experimental data; FoldX/Rosetta point mutations are well-validated for stability predictions.  
**Cons:** Requires PDB structure availability check; more complex branching logic.

### 5.3 Pattern C: "Ensemble Approach" (Research Grade)

Run both Pattern A and Pattern B, dock against both structures, use consensus scoring. Used in academic papers (e.g., the KRAS allosteric stabilizer study using AlphaFold + GROMACS + AutoDock).

### 5.4 Published Integration Examples

1. **KRAS mutation-agnostic drug discovery** (2025): AlphaFold structure → P2Rank pocket detection → AutoDock Vina → GROMACS 100 ns MD → ADMET filtering. Identified allosteric stabilizers for the Switch-I/II groove.

2. **Relay Therapeutics PI3Kα** (clinical): Full-length cryo-EM of PIK3CA mutant → long-timescale MD → computational design → zovegalisib (RLY-2608), now in clinical trials for PIK3CA-mutant cancers.

3. **AnadrosPilotMD pipeline** (GitHub: https://github.com/tdextermorse/anadrospilotmd): P2Rank + AutoDock + AmberTools + GROMACS integrated workflow with parallel execution.

---

## 6. Practical Recommendations for T2Lead

### 6.1 Short-Term (Next Sprint — Week of 2026-03-30)

1. **Use nf-core/sarek 3.8.1** as the upstream variant calling engine rather than scripting BWA-MEM2 → SAMtools → Mutect2 manually. Benefits:
   - Tested, reproducible, containerized
   - Handles the full GATK Mutect2 best practices workflow
   - Built-in VEP 115.0 annotation
   - `--input samplesheet.csv` accepts tumor/normal pairs directly

2. **Write a VCF-to-mutant-sequence parser** as a new T2Lead module (`src/drugpipe/variant_analysis/vcf_to_mutant_seq.py`):
   - Input: VEP-annotated VCF
   - Extract: missense variants in driver genes (filter by IMPACT=HIGH or MODERATE)
   - Output: mutant FASTA sequences ready for ESMFold

3. **Integration point with existing pipeline:**
   - New Stage 0 or pre-Stage 1 module
   - If a somatic VCF is provided, bypass Stage 1's target selection (the mutation IS the target)
   - Feed mutant structure directly into Stage 4

### 6.2 Medium-Term (Next 2–4 Weeks)

1. **Implement Pattern B** (known structure + FoldX fallback to ESMFold):
   - Add PDB lookup via RCSB REST API (`https://data.rcsb.org/rest/v1/core/uniprot/{uniprot_id}`)
   - Integrate FoldX 5 (academic license, free for non-commercial) for point mutations on existing structures

2. **Benchmark mutant structure quality:**
   - Compare ESMFold mutant vs. FoldX mutant vs. experimental mutant structures for PIK3CA H1047R, EGFR L858R, BRAF V600E
   - Use RMSD to experimental structure + docking pose consistency as metrics

3. **Add GPU-accelerated variant calling option:**
   - NVIDIA Clara Parabricks for users with GPU access
   - 10× speedup makes interactive/batch processing practical

### 6.3 Architecture: Proposed Module Layout

```
src/drugpipe/
├── variant_analysis/           # NEW — somatic variant calling bridge
│   ├── __init__.py
│   ├── sarek_runner.py         # Launch nf-core/sarek via Nextflow
│   ├── vcf_parser.py           # Parse VEP-annotated VCF → driver mutations
│   ├── mutant_sequence.py      # Generate mutant protein sequences
│   └── structure_bridge.py     # Route to FoldX (if PDB exists) or ESMFold
├── target_identification/      # Existing Stage 1
├── virtual_screening/          # Existing Stage 2
├── hit_to_lead/                # Existing Stage 3
└── lead_optimization/          # Existing Stage 4
```

### 6.4 Minimum Hardware Requirements

| Configuration | Suitable For | Estimated Cost |
|--------------|-------------|----------------|
| **16 cores, 32 GB RAM, 500 GB SSD** | WES variant calling + downstream | On-prem or ~$2/hr cloud |
| **32 cores, 64 GB RAM, 1 TB SSD, 1× A100** | WGS + ESMFold + MD | ~$5–8/hr cloud |
| **GPU node (4× T4/A10)** | Parabricks-accelerated WGS | ~$4–6/hr cloud (spot) |

---

## 7. Summary

| Dimension | Assessment |
|-----------|------------|
| **Technical feasibility** | High — all tools are production-grade and well-documented |
| **Novel contribution** | The variant→structure→docking bridge is unique; no open-source pipeline covers this span |
| **Biggest risk** | Mutant structure accuracy for docking — need validation against known mutant crystal structures |
| **Recommended approach** | nf-core/sarek for upstream + custom VCF parser + FoldX/ESMFold hybrid for structures |
| **Compute for WES** | ~4–6 hours on a 16-core machine; fits within T2Lead's existing workflow cadence |
| **Commercial differentiation** | T2Lead would be the first open-source pipeline connecting raw tumor FASTQ to lead compounds |
