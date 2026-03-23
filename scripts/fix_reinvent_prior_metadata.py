#!/usr/bin/env python3
"""Remove custom/invalid metadata from a REINVENT4 .prior pickle.

REINVENT4 logs ERROR with a large dict when ``metadata.hash_id`` does not match
the file content.  Stripping ``metadata`` lets REINVENT inject a fresh default
on load (same as many upstream priors).  A backup ``*.prior.bak`` is created.

Usage:
  conda run -n t2lead python scripts/fix_reinvent_prior_metadata.py /path/to/reinvent.prior
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import torch


def main() -> int:
    p = argparse.ArgumentParser(description="Strip metadata from REINVENT4 prior pickle.")
    p.add_argument("prior", type=Path, help="Path to reinvent.prior")
    p.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not write *.prior.bak (dangerous)",
    )
    args = p.parse_args()
    path: Path = args.prior.resolve()
    if not path.is_file():
        print(f"Not a file: {path}", file=sys.stderr)
        return 1

    backup = path.with_suffix(path.suffix + ".bak")
    if not args.no_backup:
        shutil.copy2(path, backup)
        print(f"Backup: {backup}")

    data = torch.load(path, map_location="cpu", weights_only=False)
    if isinstance(data, dict) and "metadata" in data:
        del data["metadata"]
        torch.save(data, path)
        print(f"Removed metadata from {path} — REINVENT4 will assign defaults on load.")
    else:
        print("No top-level 'metadata' key; nothing to do.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
