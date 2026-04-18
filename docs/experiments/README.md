# Experiment Records Convention

This directory is the canonical location for per-run experiment records.

## Directory layout

Use this pattern for every run:

`docs/experiments/runs/<date>_<run-name>_<run-id>/`

Example:

`docs/experiments/runs/2026-04-13_xena-full-max-variants_20260413_191441/`

## Required files per run

1. `REPORT.md`
   - Final run report for humans.
   - Must include: objective, config, runtime status, successes, failures, root causes, outputs, next actions.

2. `ARTIFACTS.md`
   - Machine/location-focused artifact index.
   - Must include: absolute paths, file purpose, whether complete/partial.

3. `COMMANDS.md`
   - Launch command, environment variables, and monitoring commands actually used.

4. `METADATA.yaml`
   - Structured metadata for automation (run_id, disease, sample_id, input source, start/end time, status, owner, notes).

## Status values (recommended)

- `completed`
- `interrupted`
- `failed`
- `partial-success`

## Writing rules

- Always use absolute paths for generated artifacts.
- Separate "observed facts" from "hypotheses/causes".
- List unresolved blockers clearly.
- If run is interrupted, explicitly state last log timestamp and last confirmed successful step.
