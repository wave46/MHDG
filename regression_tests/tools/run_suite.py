#!/usr/bin/env python3
"""Run a tracked workflow-by-layout regression suite."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from build_solver import build_solver, parse_build_jobs
from suite.reporting import print_run_summary
from suite.runner import run_suite
from support.errors import BundleError, HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_id", metavar="SUITE")
    parser.add_argument("--run-id")
    parser.add_argument(
        "--run-only",
        action="store_true",
        help="save solver results without comparing them",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="continue an existing suite run without repeating recorded cells",
    )
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--layouts", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--suites", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--tolerances",
        required=True,
        type=Path,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--require-bundle-class",
        choices=("candidate", "golden"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="build clean serial and parallel executables before running",
    )
    parser.add_argument("--build-jobs", type=parse_build_jobs, metavar="N")
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        settings_path = args.settings
        if args.build:
            build = build_solver(settings_path, args.repository_root, args.build_jobs)
            settings_path = build.settings_path
            print(f"using build: {build.path}")
        elif args.build_jobs is not None:
            raise BundleError("--build-jobs requires --build")
        summary_path, summary = run_suite(
            settings_path,
            args.suite_id,
            args.cases,
            args.layouts,
            args.suites,
            args.tolerances,
            args.run_id,
            args.require_bundle_class,
            compare=not args.run_only,
            resume=args.resume,
        )
    except HarnessError as exc:
        print(f"suite failed: {exc}", file=sys.stderr)
        return 1

    print_run_summary(summary, summary_path)
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
