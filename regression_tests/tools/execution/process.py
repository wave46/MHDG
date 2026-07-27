"""Launch one prepared solver process with captured logs."""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

from execution.models import RuntimeEnvironment, SolverProcessResult
from preparation.commands import openmp_environment
from preparation.models import PreparedRun
from support.environments import source_environment
from support.errors import BundleError
from support.time import utc_now


def runtime_environment(
    settings: dict[str, str],
    omp_threads: int,
) -> RuntimeEnvironment:
    """Load the configured shell environment and deterministic OpenMP values."""
    configured_path = settings.get("MHDG_ENVIRONMENT_SCRIPT")
    if configured_path:
        script = Path(configured_path).expanduser()
        if not script.is_absolute():
            raise BundleError("MHDG_ENVIRONMENT_SCRIPT must be an absolute path")
        script, values = source_environment(script)
    else:
        script = None
        values = dict(os.environ)

    openmp = openmp_environment(omp_threads)
    values.update(openmp)
    return RuntimeEnvironment(values, openmp, script)


def launch_solver(
    prepared: PreparedRun,
    environment: dict[str, str],
) -> SolverProcessResult:
    """Run the solver once and capture launch status and elapsed time."""
    started_utc = utc_now()
    started_clock = time.monotonic()
    exit_code: int | None = None
    launch_error: str | None = None
    try:
        with (prepared.path / "stdout.log").open(
            "w", encoding="utf-8"
        ) as stdout, (prepared.path / "stderr.log").open(
            "w", encoding="utf-8"
        ) as stderr:
            completed = subprocess.run(
                prepared.command,
                cwd=prepared.path,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        exit_code = completed.returncode
    except OSError as exc:
        launch_error = str(exc)
    return SolverProcessResult(
        started_utc,
        time.monotonic() - started_clock,
        exit_code,
        launch_error,
    )
