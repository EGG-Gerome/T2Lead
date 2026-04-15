---
name: experiment-reporting
description: Write standardized experiment run reports for T2Lead, including status, problems, root causes, successes, effects, generated artifacts, and next actions. Use when the user asks for experiment summary, run复盘, 问题原因分析, 结果汇报, 产物路径梳理, or postmortem documentation.
---

# Experiment Reporting

Use this skill to generate consistent run reports under `docs/experiments/runs/<date>_<run-name>_<run-id>/`.

## Output files to produce

1. `REPORT.md`
2. `ARTIFACTS.md`
3. `COMMANDS.md`
4. `METADATA.yaml`

If needed, also produce a higher-level postmortem in `docs/postmortems/`.
For Chinese-first teams, also produce `REPORT_zh.md`, `ARTIFACTS_zh.md`, and `COMMANDS_zh.md` unless the user explicitly wants English-only.

## Required report content

- Run identity: run_id, disease, input source, sample strategy/sample id.
- Runtime state: running/completed/interrupted/failed and evidence.
- Successes and effective outcomes.
- Problems encountered and likely root causes.
- Impact assessment (what was affected and what still usable).
- Generated file inventory with absolute paths and meanings.
- Explicit "next actions" checklist.

## Evidence rules

- Prefer observed log/process facts over assumptions.
- Separate facts from hypotheses.
- Always include last known log timestamp for non-completed runs.
- Use absolute paths for all key artifacts.

## Writing style

- Concise, technical, operator-friendly.
- Use bullet points for scanability.
- Keep one section for "what succeeded" and one for "what failed/why".
- Language default: Chinese for internal reports; English copies only when requested for collaborators.

## Path convention

Always store per-run records at:

`docs/experiments/runs/<date>_<run-name>_<run-id>/`

When missing, create the directory and all four files.

## Template reference

Use [report_template.md](report_template.md) for default structure.
