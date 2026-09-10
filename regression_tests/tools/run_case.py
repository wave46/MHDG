"""Prepare and execute one isolated MHDG regression run."""

from __future__ import annotations

from execution.models import RunResult
from execution.single import execute_run
from execution.staged import execute_staged_run
from regression_tests.prepare import (
    PreparedExecution,
    PreparedStagedRun,
)


def execute_prepared(
    prepared: PreparedExecution,
    settings: dict[str, str],
) -> RunResult:
    """Execute either one warm run or a sequential staged workflow."""
    if isinstance(prepared, PreparedStagedRun):
        return execute_staged_run(prepared, settings)
    return execute_run(prepared, settings)
