#!/usr/bin/env python3
"""Promote accepted regression results into a complete golden bundle."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bundle import promotion
from support.errors import HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_summary", metavar="SUITE_SUMMARY", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bundle-version", required=True, metavar="VERSION")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        result = promotion.promote_bundle(
            args.settings,
            args.suite_summary,
            args.output,
            args.bundle_version,
            args.cases,
        )
    except HarnessError as exc:
        print(f"golden bundle promotion failed: {exc}", file=sys.stderr)
        return 1

    print(f"golden bundle created: {args.output.expanduser().resolve()}")
    print(f"bundle: {result.bundle_id} version {result.bundle_version}")
    print(
        f"verified {result.artifact_count} artifacts "
        f"({result.verified_bytes} bytes)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
