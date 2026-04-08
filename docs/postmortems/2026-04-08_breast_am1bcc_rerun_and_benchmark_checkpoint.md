# Breast Stage 4 AM1-BCC Rerun and Benchmark Checkpoint Postmortem

Status: Completed
Date: 2026-04-08
Scope: breast cancer standard pipeline, selected target `PIK3CA (CHEMBL4005)`, Stage 4 + benchmark/report regeneration

## Summary

The breast Stage 4 shortlist could not be treated as final because many MD jobs were being skipped for environment reasons, not scientific reasons. The immediate cause was that the real runtime did not always expose the AmberTools executables required by OpenFF AM1-BCC charge assignment. During validation of the rerun, a second issue was found: benchmark MD could silently reuse lead checkpoint rows from `md_progress.csv`, contaminating benchmark comparator values.

This postmortem covers both incidents because they affected the same delivery: the breast `PIK3CA` rerun, the exported lead shortlist, and the benchmark comparison shown in the dashboard.

## User-visible impact

- The earlier 7-lead shortlist was incomplete as a fair "full MD compared" result because some molecules were skipped before real MD due to charge-backend failure.
- Benchmark comparator MD values were temporarily untrustworthy because some values came from the lead checkpoint file rather than from fresh benchmark simulations.
- The dashboard at the user-visible path could lag behind the actual regenerated file because the pipeline wrote under `/autodl-fs/data/...` while the user checked `/root/autodl-fs/...`.

## Root causes

### 1. AM1-BCC backend not guaranteed in the real runtime

OpenFF small-molecule parameterization ultimately calls `assign_partial_charges("am1bcc")`. That path requires AmberTools binaries such as `antechamber` and `sqm`. The environment contained the Python packages needed for MD, but the actual execution mode did not reliably guarantee those executables were on `PATH`, so AM1-BCC parameterization could fail at runtime.

### 2. MD startup did not fail loudly enough

Before the fix, MD setup could reach the charge-assignment stage and then fail per molecule, which looked like a molecule problem even when the actual issue was a broken environment. This made the run appear to progress while silently reducing coverage.

### 3. Benchmark MD reused the lead checkpoint namespace

`MDSimulator` stores resumable results in `stage4_optimization/md_trajectories/md_progress.csv`, keyed by the DataFrame row index (`idx`). Benchmark runs also use small integer row indices (`0..N-1`). Because the benchmark path shared the same checkpoint file, benchmark rows could collide with lead rows and inherit lead MD results.

### 4. Output-path aliasing made dashboard verification confusing

The generated dashboard was written under `/autodl-fs/data/T2Lead/...`, while the user inspected `/root/autodl-fs/T2Lead/...`. Without an explicit sync, the visible HTML could stay stale even when the real report payload had already been regenerated.

## Changes made

### Runtime and code fixes

- `src/drugpipe/lead_optimization/md_simulation.py`
  - added an AM1-BCC preflight probe for AmberTools binaries and charge assignment
  - added a fail-loud requirement so MD stops early when the charge backend is unavailable
- `Dockerfile`
  - added `ambertools` to the reproducible conda install set
- `Makefile`
  - added `ambertools` to the reproducible conda install target
- `src/drugpipe/benchmark/scorer.py`
  - forced benchmark MD to run with `resume_from_checkpoint = False` so it cannot reuse the lead checkpoint file

### Operational cleanup and rerun

- removed stale breast Stage 4 artifact directories and old manual benchmark logs
- deleted the old lead MD checkpoint before rerunning breast Stage 4
- reran breast Stage 4 from scratch for `CHEMBL4005` with `device=cuda`
- regenerated benchmark/report outputs after the checkpoint-isolation fix
- synced the visible dashboard path after regeneration

## Validation evidence

### AM1-BCC path

- The rerun completed without the earlier AM1-BCC/OpenFF charge-chain failures that had been skipping molecules for environment reasons.
- The new runtime now checks AM1-BCC availability before spending hours in Stage 4.

### Benchmark contamination fix

- A fresh breast benchmark rerun produced 4 comparator rows.
- The regenerated comparator MD values no longer match the lead checkpoint rows from `md_progress.csv`.
- This confirmed that the benchmark numbers were no longer coming from lead checkpoint reuse.

### Final outputs

- `optimized_leads.csv`: 10 exported rows
- `strict_reliable_count`: 13 across 45 screened candidates
- dashboard payload verified at the user-visible path
- top generated lead `rank_score`: `0.843944`
- top benchmark drug: `ALPELISIB` with `rank_score = 0.815535`

## Timing

- Breast Stage 4 clean rerun: about `10h 04m 48s`
- Breast benchmark/report clean regeneration after checkpoint fix: about `44m 04s`

## What worked well

- Conservative cleanup kept the rerun targeted to the breast `PIK3CA` path instead of broad destructive deletion.
- The fail-loud AM1-BCC preflight turns an ambiguous late failure into an immediate actionable error.
- The benchmark checkpoint collision was caught before packaging the run as a final deliverable.

## What did not work well

- Sharing a checkpoint file keyed only by local row index was too fragile once multiple MD call sites existed in the same Stage 4 directory.
- Report regeneration across `/autodl-fs/data` and `/root/autodl-fs` remained operationally confusing.
- Structure-package assets were not regenerated automatically as part of the Stage 4/report flow.

## Follow-ups

1. Replace checkpoint identity based on local DataFrame index with a stable molecule key, such as canonical-SMILES hash plus run context.
2. Separate lead and benchmark MD checkpoint locations so collisions are impossible by construction.
3. Add a small post-run packaging/export helper for the doctor package so 2D assets and tables stay in sync with the latest CSVs.
4. Keep the current MM-GBSA overflow gating behavior (`DeltaG > 50` treated as failed), but document it clearly as a scientific reliability gate rather than an environment problem.
