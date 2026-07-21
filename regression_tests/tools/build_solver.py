#!/usr/bin/env python3
"""Build the serial and parallel regression executables reproducibly."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from build import configuration, workflow
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
        "--jobs",
        type=configuration.parse_build_jobs,
        metavar="N",
        help="parallel jobs for each make invocation",
    )
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        result = workflow.build_solver(args.settings, args.repository_root, args.jobs)
    except HarnessError as exc:
        print(f"build failed: {exc}", file=sys.stderr)
        return 1

    print(f"build completed: {result.path}")
    print(f"generated settings: {result.settings_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
