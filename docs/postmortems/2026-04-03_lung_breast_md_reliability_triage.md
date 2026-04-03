# Lung/Breast MD Reliability Triage (2026-04-03)

Status: Draft (live run + mini-validation pending implementation)

## Scope and assumptions

- Scope includes:
  - Lung run `20260403_105734` (Stage 4 + benchmark already completed).
  - Breast run started at the same time (`/tmp/breast_stage4.log`), still incomplete at the time of this draft.
- This note focuses on MD reliability and downstream scoring behavior (`md_binding_energy`, `md_rmsd_mean`, `opt_score`).
- We treat docking and ADMET as upstream signals; we do not claim medicinal efficacy from this triage.

## Timeline and observed state

### Lung (complete)

- Stage 4 and benchmark finished; dashboard generated.
- Repeated pattern in Stage 4 and benchmark:
  - many `DeltaG` values very close to 0,
  - many `md_rmsd_mean` values high (often > 20 A),
  - resulting in weak structural discrimination in MD-derived score components.

### Breast (incomplete at observation time)

- Process is no longer running, and log has no `Stage 4 complete` / `PIPELINE COMPLETE`.
- Parsed from `/tmp/breast_stage4.log`:
  - attempted molecules: 8 (`[1, 2, 6, 7, 17, 28, 39, 42]`)
  - completed MM-GBSA molecules: 5 (`[1, 2, 6, 39, 42]`)
  - overflow failures: 2 (`[7, 28]`)
  - started but not finished in log: 1 (`[17]`)
  - completed subset stats:
    - `DeltaG`: min `-23.90`, median `-18.27`, max `-17.49` kcal/mol
    - `md_rmsd_mean`: min `7.46`, median `9.16`, max `12.29` A
    - count with `md_rmsd_mean > 15 A`: `0`

Interpretation: breast currently shows mostly healthy completed MD samples but is not end-to-end complete yet; conclusions remain provisional until completion.

## Problem statement

For Lung, MD-derived terms appear low-signal for ranking in this batch. The root-cause hypothesis is simulation reliability (e.g., ligand drift / conditioning), not merely score normalization.

## Pre-registered mini-validation protocol

Goal: decide whether to adopt restrained implicit MD and/or explicit-solvent fallback for Lung-like cases.

### Molecules

- Use 2 Lung molecules from current top-10:
  - one high-drift case (worst `md_rmsd_mean`),
  - one mid-drift case (near median `md_rmsd_mean`).

### Conditions (A/B/C per molecule)

- A: baseline implicit MD (current config).
- B: implicit MD + mild ligand heavy-atom positional restraint.
- C: explicit-solvent MD (TIP3P/PME existing path).

### Metrics

- Primary:
  - `md_binding_energy`
  - `md_binding_energy_std`
  - `md_rmsd_mean`
  - overflow/failure status
- Secondary:
  - runtime per molecule

### Pre-defined acceptance thresholds (fixed before running)

- Reliability gate:
  - no overflow in test molecule result, and
  - finite `md_binding_energy`, finite `md_rmsd_mean`.
- Structural quality target:
  - `md_rmsd_mean <= 15 A` for at least one previously high-drift test molecule, or
  - at least 25% RMSD reduction vs baseline on both test molecules.
- Energy plausibility target:
  - `md_binding_energy <= -1.0 kcal/mol` on both molecules (avoid near-zero collapse).
- Runtime target for default method:
  - median runtime <= 3x baseline implicit.
  - If better than 3x but clearly more stable, keep as top-k fallback (not default).

### Decision rule

- Choose restrained implicit as default for affected targets if it passes reliability + quality targets with acceptable runtime.
- Choose explicit-solvent as fallback/top-k rescue when implicit remains unstable.
- Do not select a fix based only on cosmetic score scaling.

## Planned score guardrails

- Add `md_reliable` boolean to outputs and dashboard payload.
- When `md_reliable` is false, keep `opt_score` visible but use a non-MD fallback ranking signal (`fast_score`) for decision support.
- Add explicit dashboard annotation when quality gating is active.

## Open items

- Complete mini-validation runs (A/B/C) and append measured results.
- Complete breast run (or re-run with `nohup`) and finalize breast quality summary.
- Implement and verify code changes for reliability gating and fallback ranking.

## Mini-validation results (completed)

Artifacts:
- `/autodl-fs/data/T2Lead/lung_cancer/stage4_optimization/mini_validation/current_v3/summary.csv`
- `/tmp/lung_md_mini_validation_v3.log`

Setup:
- Molecules: `48` (high drift), `7` (mid drift)
- Conditions:
  - implicit baseline
  - implicit + restraint (`k=1.0 kcal/mol/A^2`, heavy atoms)
  - explicit solvent (shortened mini run)

Observed metrics:
- Implicit baseline:
  - `DeltaG` median: `0.0185` kcal/mol
  - `md_rmsd_mean` median: `12.72 A`
  - elapsed: `246.7 s`
- Implicit restrained:
  - `DeltaG` median: `0.0325` kcal/mol
  - `md_rmsd_mean` median: `0.48 A`
  - elapsed: `253.0 s`
- Explicit solvent:
  - `explicit_rmsd_mean` median: `6.39 A`
  - elapsed: `1471.9 s` (~6x slower than implicit)
  - `explicit_binding_energy` is very large negative total-system potential energy, not directly comparable to implicit `DeltaG`.

Interpretation:
- Restraints drastically reduce RMSD with little runtime cost, but do not fix near-zero implicit `DeltaG` collapse.
- Explicit solvent improves RMSD without over-constraining, but runtime is much higher and energy term is not directly comparable to implicit MM-GBSA.

## Fix-path decision (selected)

1. Immediate production path:
   - keep implicit MD as default,
   - add MD reliability gate (`md_reliable`) and fallback ranking to `fast_score` when MD is unreliable,
   - expose reliability status in CSV/dashboard to avoid over-trusting weak MD terms.
2. Restraints:
   - keep available as optional target-specific control,
   - do not force-enable globally yet.
3. Explicit solvent:
   - keep as optional top-k refinement/fallback due runtime,
   - avoid treating current `explicit_binding_energy` as directly comparable to implicit `DeltaG`.
