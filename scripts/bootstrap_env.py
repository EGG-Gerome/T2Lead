#!/usr/bin/env python3
"""Update .env with DP_* paths after `make install` (idempotent)."""
from __future__ import annotations

import argparse
from pathlib import Path


MANAGED_KEYS = (
    "DP_HIT_TO_LEAD__REINVENT4__REINVENT_PATH",
    "DP_HIT_TO_LEAD__REINVENT4__PRIOR_PATH",
    "DP_HIT_TO_LEAD__ANALOG_GEN__CREM_DB_PATH",
)


def _strip_managed(lines: list[str]) -> list[str]:
    out: list[str] = []
    skip_block = False
    for line in lines:
        stripped = line.strip()
        if stripped == "# >>> t2lead make install >>>":
            skip_block = True
            continue
        if stripped == "# <<< t2lead make install <<<":
            skip_block = False
            continue
        if skip_block:
            continue
        if any(stripped.startswith(k + "=") for k in MANAGED_KEYS):
            continue
        out.append(line.rstrip("\n"))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--env-file", type=Path, required=True)
    p.add_argument("--conda-prefix", required=True)
    p.add_argument("--prior-path", required=True)
    p.add_argument("--crem-db", required=True)
    args = p.parse_args()

    reinvent_bin = Path(args.conda_prefix) / "bin" / "reinvent"
    block_lines = [
        "",
        "# >>> t2lead make install >>>",
        "# Paths below are updated by `make install` — safe to re-run.",
        f"{MANAGED_KEYS[0]}={reinvent_bin}",
        f"{MANAGED_KEYS[1]}={args.prior_path}",
        f"{MANAGED_KEYS[2]}={args.crem_db}",
        "# <<< t2lead make install <<<",
        "",
    ]

    text = args.env_file.read_text(encoding="utf-8")
    lines = text.splitlines()
    kept = _strip_managed(lines)
    while kept and kept[-1] == "":
        kept.pop()
    new_text = "\n".join(kept) + "\n" + "\n".join(block_lines)
    args.env_file.write_text(new_text, encoding="utf-8")
    print(f"Updated {args.env_file} with REINVENT4 and CReM paths.")


if __name__ == "__main__":
    main()
