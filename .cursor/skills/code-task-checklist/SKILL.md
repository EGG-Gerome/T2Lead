---
name: code-task-checklist
description: Implement code tasks with small-step delivery, immediate testing, and git scope checks. Use when writing or refactoring code in T2Lead.
---

# Code Task Checklist

## Workflow

1. Clarify the smallest implementable change.
2. Implement that change only.
3. Run targeted test/smoke for the changed scope.
4. Check `git status` and `git diff`.
5. Repeat for next small change.

## Safety notes

- Prefer minimal, verifiable edits over large rewrites.
- Report facts from logs/tests; mark assumptions explicitly.
