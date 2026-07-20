#!/usr/bin/env python3
"""Compare a completed run using its declared workflow policy."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

from comparison.reporting import print_adaptive_summary, print_fixed_summary
from comparison.workflow import compare_completed_run
from support.errors import HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", metavar="RUN_DIRECTORY", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument(
        "--tolerance-profile",
        dest="tolerance_profile",
    )
    parser.add_argument("--report", type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--tolerances", required=True, type=Path, help=argparse.SUPPRESS
    )
    args = parser.parse_args(argv)

    try:
        policy, report_path, report = compare_completed_run(
            args.run_directory,
            args.cases,
            args.tolerances,
            candidate_override=args.candidate,
            reference_override=args.reference,
            tolerance_profile_override=args.tolerance_profile,
            report_override=args.report,
        )
    except HarnessError as exc:
        print(f"comparison failed: {exc}", file=sys.stderr)
        return 1

    print_comparison(policy, report_path, report)
    return 0 if report["status"] == "passed" else 1


def print_comparison(
    policy: str,
    report_path: Path,
    report: dict[str, Any],
) -> None:
    """Print the appropriate concise summary for a comparison report."""
    print(f"comparison policy: {policy}")
    print(f"comparison {report['status']}: {report_path}")
    if policy == "fixed_hdf5":
        _print_inputs(report)
        print_fixed_summary(report)
    elif policy == "mesh_independent":
        _print_inputs(report, report["files"])
        print_adaptive_summary(report)
    else:
        print(
            f"stages checked: {report['checked_stage_count']}/"
            f"{report['total_stage_count']}"
        )
        for failure in report["failures"]:
            print(f"FAIL: {failure}")


def _print_inputs(
    report: dict[str, Any],
    files: dict[str, str] | None = None,
) -> None:
    files = files or report
    print(f"candidate: {files['candidate']}")
    print(f"reference: {files['reference']}")
    print(f"tolerance profile: {report['tolerance_profile']['id']}")


if __name__ == "__main__":
    raise SystemExit(main())
