#!/usr/bin/env python3
"""Verify previously saved workflow-by-layout suite results."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from suite.reporting import print_verification_summary
from suite.verification import verify_suite
from support.errors import HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_summary", metavar="SUITE_SUMMARY", type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--tolerances",
        required=True,
        type=Path,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        path, summary = verify_suite(
            args.suite_summary,
            args.cases,
            args.tolerances,
        )
    except HarnessError as exc:
        print(f"suite verification failed: {exc}", file=sys.stderr)
        return 1

    print_verification_summary(summary, path)
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
