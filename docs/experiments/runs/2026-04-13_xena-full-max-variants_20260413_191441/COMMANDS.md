# Commands Used

## Launch intent (single GPU)

Environment:

- `DP_PIPELINE__DEVICE=cuda`
- `PARALLEL_WORKERS=1`

Core run arguments:

- `--vcf-path /root/autodl-fs/T2Lead_mainpath/data/test_inputs/xena/TCGA-BRCA.somaticmutation_wxs.tsv`
- `--xena-sample-strategy max_variants`
- `--sample-id xena_full_max_variants`
- shared-cache mode enabled (default shared ChEMBL behavior in current script logic)

## Monitoring commands

- Process check:
  - `pgrep -af "run_mainpath|python.*drugpipe"`
- Log check:
  - inspect `.../logs/*_full.log`
- Freshness check:
  - `stat .../logs/*_full.log`

## Recommended next launch mode (for long run)

Use `tmux` or `nohup` to avoid interruption when SSH disconnects.
