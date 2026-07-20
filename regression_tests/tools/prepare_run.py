#!/usr/bin/env python3
"""Prepare an isolated MHDG regression run without executing the solver."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

from preparation.configuration import load_preparation_inputs
from preparation.models import PreparedExecution, PreparedRun
from preparation.workflows import prepare_staged_run, prepare_warm_run
from support.errors import BundleError, HarnessError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id", metavar="CASE")
    parser.add_argument("workflow_id", metavar="WORKFLOW")
    parser.add_argument("--layout", required=True, dest="layout_id")
    parser.add_argument("--run-id")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--layouts", required=True, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        prepared = prepare_run(
            args.settings,
            args.case_id,
            args.workflow_id,
            args.layout_id,
            args.cases,
            args.layouts,
            args.run_id,
        )
    except HarnessError as exc:
        print(f"run preparation failed: {exc}", file=sys.stderr)
        return 1

    print(f"run prepared: {prepared.path}")
    if isinstance(prepared, PreparedRun):
        print(f"OMP_NUM_THREADS={prepared.omp_threads}")
        print(f"command: {shlex.join(prepared.command)}")
    else:
        for stage in prepared.stages:
            print(f"stage {stage.stage_id}: {shlex.join(stage.run.command)}")
    return 0


def prepare_run(
    settings_path: Path,
    case_id: str,
    workflow_id: str,
    layout_id: str,
    case_dir: Path,
    layouts_path: Path,
    run_id: str | None = None,
    validate_bundle: bool = True,
) -> PreparedExecution:
    """Create one validated, isolated run or staged workflow directory."""
    inputs = load_preparation_inputs(
        settings_path,
        case_id,
        workflow_id,
        layout_id,
        case_dir,
        layouts_path,
        run_id,
        validate_bundle,
    )
    workflow_kind = inputs.workflow["kind"]
    if workflow_kind == "warm_same_state":
        return prepare_warm_run(inputs)
    if workflow_kind in {"staged_fixed_mesh", "staged_adaptive_mesh"}:
        return prepare_staged_run(inputs)
    raise BundleError(
        f"run preparation does not support workflow kind {workflow_kind}"
    )


if __name__ == "__main__":
    raise SystemExit(main())
