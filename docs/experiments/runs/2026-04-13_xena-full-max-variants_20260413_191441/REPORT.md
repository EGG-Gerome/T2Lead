# Run Report: BRCA Xena Full Max Variants

- status: `interrupted`
- run_id: `20260413_191441`
- sample_id: `TCGA-AN-A046-01A`
- strategy: `xena-sample-strategy=max_variants`
- runtime mode: `parallel_workers=1`, `device=cuda`

## Objective

Run mutation-driven pipeline on the selected Xena sample using best-practice structure handling:

1. experimental PDB first,
2. point mutation modeling,
3. local restrained optimization,
4. ESMFold fallback when template-based mutation is unavailable.

## Outcome snapshot

- Input conversion succeeded.
- Variant analysis produced large-scale FASTA and structure artifacts.
- Pipeline did not reach final completion; process was not active during follow-up check.

## Major successes

- Xena TSV -> VEP-VCF conversion worked.
- Multiple variants successfully produced `*_mut_localopt.pdb` outputs.
- Shared cache usage was in place and avoided redundant downloads for key cache directories.

## Major issues and causes

1. **Run interruption before completion**
   - Last log write: `2026-04-13 22:13:35 +0800`.
   - No active `run_mainpath` process found later.
   - Most likely a session/process interruption rather than deterministic code crash.

2. **External API coverage/limit issues**
   - RCSB search occasional `400`.
   - ESMFold API `413` for long proteins.
   - AlphaFold DB `404` for some accessions.
   - Effects: per-variant structure-resolution gaps for a subset of mutations.

## Last confirmed successful step

Downloaded alternate experimental template `6ID1.pdb` for variant `SYF2 K71T` after `8C6J` download failures.

## Next actions

1. Relaunch under `tmux`/`screen`/`nohup`.
2. Keep `parallel_workers=1` on current single-GPU host.
3. Prefer chunked/batched variant execution for better operational control and resumability.
