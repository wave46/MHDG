#!/usr/bin/env python3
"""Create a validated MHDG regression bundle from prepared case files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bundle import creation
from support.errors import HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        required=True,
        dest="case_id",
        metavar="CASE",
        help="tracked case identifier",
    )
    parser.add_argument(
        "--source",
        required=True,
        type=Path,
        help="prepared case directory",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="new bundle path",
    )
    parser.add_argument(
        "--bundle-version",
        default="1.0.0",
        metavar="VERSION",
        help="bundle version label (default: %(default)s)",
    )
    parser.add_argument(
        "--cases",
        required=True,
        type=Path,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        summary = creation.create_bundle(
            args.case_id,
            args.source,
            args.output,
            args.cases,
            args.bundle_version,
        )
    except HarnessError as exc:
        print(f"bundle creation failed: {exc}", file=sys.stderr)
        return 1

    print(f"bundle created: {args.output}")
    print(
        f"recorded {summary.artifact_count} artifacts "
        f"({summary.verified_bytes} bytes)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
