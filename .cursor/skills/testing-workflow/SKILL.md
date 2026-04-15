---
name: testing-workflow
description: Run layered tests for T2Lead changes: unit, smoke, and run-level checks. Use when validating code edits or pipeline run readiness.
---

# Testing Workflow

## Validation layers

1. Unit tests for changed module.
2. Smoke command for affected pipeline path.
3. Output-path and log sanity checks.

## Minimum evidence to report pass

- Command used
- Exit code
- Key output/log line

## Failure handling

- Stop, capture evidence, and propose smallest next fix.
