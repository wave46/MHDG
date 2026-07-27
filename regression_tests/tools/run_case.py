#!/usr/bin/env python3
"""Prepare and execute one isolated MHDG regression run."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

from bundle.settings import read_settings
from execution.models import RunResult
from execution.single import execute_run
from execution.staged import execute_staged_run
from preparation.models import (
    PreparedExecution,
    PreparedRun,
    PreparedStagedRun,
)
from prepare_run import prepare_run
from support.errors import HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id", metavar="CASE")
    parser.add_argument("workflow_ids", metavar="WORKFLOW", nargs="+")
    parser.add_argument("--layout", required=True, dest="layout_id")
    parser.add_argument("--run-id")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--layouts", required=True, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        settings = read_settings(args.settings)
    except HarnessError as exc:
        print(f"run failed: {exc}", file=sys.stderr)
        return 1

    completed = True
    for workflow_id in args.workflow_ids:
        print(f"workflow: {workflow_id}")
        try:
            prepared = prepare_run(
                args.settings,
                args.case_id,
                workflow_id,
                args.layout_id,
                args.cases,
                args.layouts,
                args.run_id,
            )
            _print_prepared(prepared)
            result = execute_prepared(prepared, settings)
        except HarnessError as exc:
            print(f"run failed: {exc}", file=sys.stderr)
            completed = False
            continue
        _print_result(result)
        completed = completed and result.status == "completed"
    return 0 if completed else 1


def execute_prepared(
    prepared: PreparedExecution,
    settings: dict[str, str],
) -> RunResult:
    """Execute either one warm run or a sequential staged workflow."""
    if isinstance(prepared, PreparedStagedRun):
        return execute_staged_run(prepared, settings)
    return execute_run(prepared, settings)


def _print_prepared(prepared: PreparedExecution) -> None:
    print(f"run directory: {prepared.path}")
    if isinstance(prepared, PreparedRun):
        print(f"command: {shlex.join(prepared.command)}")
        print("solver output: stdout.log and stderr.log", flush=True)
    else:
        print(f"stages: {len(prepared.stages)}", flush=True)


def _print_result(result: RunResult) -> None:
    print(f"status: {result.status}")
    print(f"solver exit code: {result.exit_code}")
    print(f"runtime: {result.duration_seconds:.3f} s")
    print(f"HDF5 outputs: {len(result.hdf5_outputs)}")


if __name__ == "__main__":
    raise SystemExit(main())
