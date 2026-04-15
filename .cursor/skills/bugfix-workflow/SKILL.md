---
name: bugfix-workflow
description: Fix bugs with evidence-first debugging and regression-safe validation. Use when diagnosing failures, log anomalies, or pipeline regressions.
---

# Bugfix Workflow

## Steps

1. Reproduce and capture exact symptom (log, traceback, wrong output).
2. Isolate likely cause and confidence level.
3. Implement minimal fix.
4. Re-run targeted validation.
5. Add/adjust test to prevent regression.
6. Check `git status` and `git diff` before handoff.

## Reporting

- Separate observed facts from inferred causes.
- Include last known good step and failure boundary.
