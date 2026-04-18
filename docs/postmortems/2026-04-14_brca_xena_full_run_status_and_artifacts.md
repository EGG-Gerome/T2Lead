# BRCA Xena Full-Mutation Run Status and Artifacts

Status: Interrupted (not completed)
Date: 2026-04-14
Run ID: `20260413_191441`
Primary log: `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/logs/20260413_191442_breast_cancer_full.log`

## Executive summary

This run used the full-mutation Xena TSV path with single-GPU settings (`parallel_workers=1`) and the "experimental structure + point mutation + local optimization, else ESMFold fallback" workflow. The job started normally and produced substantial variant-analysis artifacts, but it stopped before pipeline completion. No active `run_mainpath` process was found during follow-up checks.

Confidence note: the "why process stopped" part is a high-probability inference, not a proven root cause, because no terminal-level kill message or Python traceback was captured in available logs.

## What succeeded

- Xena TSV auto-conversion worked and produced a VEP-style input VCF.
- Large-scale mutant FASTA generation worked (thousands of entries).
- Many structure-resolution tasks succeeded via experimental PDB download.
- Point-mutation modeling + local restrained optimization succeeded for multiple variants (mutant local-opt PDBs generated).
- Shared-cache strategy was correctly in place (ChEMBL/report cache symlinked to shared path in the mainpath workspace).

## Main problems observed, with likely causes

1. **Run not completed / process no longer running**
   - Symptom: log stopped updating at `2026-04-13 22:13:35 +0800`, no live pipeline process.
   - Likely cause: foreground session interruption (SSH disconnect, terminal closure, manual interruption, or system-level termination).

2. **RCSB search/query occasional failures (HTTP 400)**
   - Symptom: intermittent `POST /rcsbsearch/v2/query` `400`.
   - Likely cause: query-specific server rejection or transient API-side validation mismatch for some entries.
   - Impact: fallback path used; some variants still succeeded through alternative PDB candidates or later attempts.

3. **ESMFold rejection for long sequences (HTTP 413)**
   - Symptom: `POST /foldSequence/v1/pdb/` returned `413 Request Entity Too Large`.
   - Likely cause: sequence-length/request-size limits in external ESMFold API for long proteins.
   - Impact: those variants required fallback; if fallback also unavailable, structure unresolved.

4. **AlphaFold DB missing entries (HTTP 404) for some UniProt IDs**
   - Symptom: `AF-<UniProt>-F1-model_v4.pdb` not found.
   - Likely cause: no public AlphaFold DB model available for that accession/version.
   - Impact: "No structure resolved" for affected variants unless another source succeeded.

5. **Some candidate PDB IDs unavailable (HTTP 404)**
   - Symptom: example `8C6J.pdb` not found, then switched to `6ID1.pdb` successfully.
   - Likely cause: candidate entry unavailable in PDB download endpoint (obsolete/superseded/path inconsistency).
   - Impact: minor delay from retries; not fatal when alternate entry exists.

## Effectiveness of this run attempt

- The new structure bridge strategy is functioning in real workload conditions (not just unit/smoke scale), confirmed by generated `*_mut_localopt.pdb` files.
- Failure modes were mostly external-data/API coverage constraints rather than core workflow logic failure.
- However, this run did not reach downstream full completion, so end-to-end Stage2/3/4 outputs for the full 4k+ variant set were not produced yet.

## Key generated files and what they are

### Pipeline-run root
- `/root/autodl-fs/T2Lead_mainpath/breast_cancer/variant_runs/xena_full_max_variants/20260413_191441/`
  - Run root directory for this attempt.

### Logs
- `.../logs/20260413_191442_breast_cancer_full.log`
  - Detailed runtime log (debug/info/warning).
- `.../logs/20260413_191442_breast_cancer_summary.log`
  - High-level summary log.

### Variant input conversion
- `.../variant_analysis/converted_inputs/xena_converted.TCGA-AN-A046-01A.vep.vcf`
  - Xena TSV converted into VEP-style VCF for downstream variant parsing.

### Mutant sequence artifacts
- `.../variant_analysis/mutant_fastas/`
  - Per-variant mutant FASTA files (4638 files generated in this attempt).

### Structure artifacts
- `.../variant_analysis/structures/`
  - Experimental template PDBs (e.g., `6e1f.pdb`, `9ejy.pdb`).
  - Mutant local-optimization structures (e.g., `6e1f_R251Q_mut_localopt.pdb`, `9ejy_G189R_mut_localopt.pdb`).
  - ESMFold outputs for some variants (e.g., `AGMAT_E320K_esmfold.pdb`).

## Recommendations for next run

1. Launch long run under `tmux`/`screen`/`nohup` to prevent SSH-disconnect interruption.
2. Keep `parallel_workers=1` for current 1-GPU machine.
3. For throughput realism, run in batches (gene-priority subsets or chunked variant lists) rather than one huge full-mutation set.
4. Keep this postmortem linked to a per-run report under `docs/experiments/runs/<run-id>/`.
