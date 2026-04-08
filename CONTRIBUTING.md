# Contributing to T2Lead

## Workflow

- Branch from the latest `main`; do not push directly to `main`.
- Preferred branch names: `feat/*`, `fix/*`, `docs/*`, `chore/*`.
- Run local checks before pushing:
  - `make test`
  - or `python -m pytest tests/ -v`
- Push your branch and open a pull request to `main`.
- Wait for the GitHub Actions `CI` workflow to pass before merging.

## CI

- The `CI` workflow runs on every pull request and every branch push.
- It includes:
  - `quality-checks`: dependency sanity check and Python syntax compilation
  - `unit-tests`: pytest on Python `3.9`, `3.11`, and `3.12`
- The local equivalent of the main CI gate is `make test`.

## Branch Protection

- When configuring required status checks, use checks from the `CI` workflow.
- Typical required checks are `CI / quality-checks` and the `CI / unit-tests (...)` jobs you want to enforce.
