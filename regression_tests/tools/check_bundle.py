#!/usr/bin/env python3
"""Read-only validation for an external MHDG regression data bundle."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bundle import validation
from bundle.models import ValidationSummary
from support.errors import HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--settings",
        required=True,
        type=Path,
        help="local regression settings file",
    )
    parser.add_argument(
        "--cases",
        required=True,
        type=Path,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        summary = validation.validate_bundle(args.settings, args.cases)
    except HarnessError as exc:
        print(f"bundle validation failed: {exc}", file=sys.stderr)
        return 1

    _print_summary(summary)
    return 0


def _print_summary(summary: ValidationSummary) -> None:
    print(f"bundle valid: {summary.bundle_id} version {summary.bundle_version}")
    print(
        f"verified {summary.verified_artifact_count} of "
        f"{summary.artifact_count} artifacts ({summary.verified_bytes} bytes)"
    )
    print(f"required roles checked for case: {summary.case_id}")
    for warning in summary.warnings:
        print(f"warning: {warning}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
