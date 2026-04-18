---
name: t2lead-run-guardrails
description: Enforce T2Lead runtime guardrails for environment, testing, and safety-focused reporting. Use when running pipeline commands, validating code changes, checking run health, or summarizing biomedical experiment outcomes.
---

# T2Lead Run Guardrails

## Mandatory execution rules

1. Use `t2lead` environment for project execution commands.
2. Prefer `conda run -n t2lead <command>` over implicit shell state.
3. Before long runs, confirm environment and key runtime settings.
4. For jobs expected to exceed 10 minutes, launch in `tmux`/`screen`/`nohup` by default unless the user explicitly requests foreground.
5. For large variant sets, do a smoke subset first, then run in chunks/batches.

## Mandatory verification rules

1. After non-trivial code edits, run targeted tests or smoke checks.
2. Do not report "done" without command/log evidence.
3. Separate observed facts from inferred causes in all run summaries.
4. Work in small increments: complete one small change, test it, then continue.
5. After meaningful edits, run `git status` and `git diff` to verify scope.

## Biomedical safety emphasis

- Favor conservative statements.
- Highlight uncertainty and residual risk explicitly.
- Prefer reproducible commands and absolute paths.

## Long-run operations checklist

- Use `tmux`/`screen`/`nohup` for jobs expected to outlive SSH sessions.
- Provide process-check and log-tail commands.
- If interrupted, record last successful step and last log timestamp.
- Auto-create/update run docs at `docs/experiments/runs/<date>_<run-name>_<run-id>/`.
