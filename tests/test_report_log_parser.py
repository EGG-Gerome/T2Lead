"""Log parser for HTML report."""
from __future__ import annotations

from pathlib import Path

from drugpipe.report.log_parser import Severity, issues_to_json_serializable, parse_full_log


def test_parse_full_log_groups_benign(tmp_path: Path):
    log = tmp_path / "t_full.log"
    log.write_text(
        "2026-01-01 00:00:00 [WARNING] openmm: Did not recognize residue UNL\n"
        "2026-01-01 00:00:01 [WARNING] openmm: Did not recognize residue UNL\n"
        "2026-01-01 00:00:02 [WARNING] pint.util: Redefining '[electric_potential]'\n"
        "2026-01-01 00:00:03 [ERROR] pipeline: failed badly\n",
        encoding="utf-8",
    )
    issues = parse_full_log(log)
    assert any(i.severity == Severity.ERROR for i in issues)
    assert any(i.severity == Severity.INFO and i.count >= 2 for i in issues)
    ser = issues_to_json_serializable(issues)
    assert any(x["collapsible"] for x in ser)
