"""HTML dashboard report generation from pipeline outputs."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


def generate_report(
    cfg: Dict[str, Any],
    run_root: Path,
    layout: Dict[str, Path],
    stage4_dir: Optional[Path] = None,
) -> Optional[Path]:
    """Write ``dashboard.html`` (and payload JSON) under *run_root*.

    For per-variant Stage 4 (when *stage4_dir* is a subdirectory of the stage-4
    root), files are written under that folder instead so the disease-level
    ``dashboard.html`` is left unchanged.
    """
    from drugpipe.report.html_builder import write_dashboard

    path = write_dashboard(cfg, Path(run_root), layout, stage4_dir=stage4_dir)
    logger.info("Dashboard written: %s", path)
    return path
