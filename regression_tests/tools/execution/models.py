"""Typed results produced while executing prepared solver runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class RunResult:
    path: Path
    status: str
    exit_code: int | None
    duration_seconds: float
    hdf5_outputs: list[str]


@dataclass(frozen=True)
class RuntimeEnvironment:
    values: dict[str, str]
    openmp: dict[str, str]
    setup_script: Path | None


@dataclass(frozen=True)
class SolverProcessResult:
    started_utc: str
    duration_seconds: float
    exit_code: int | None
    launch_error: str | None


@dataclass(frozen=True)
class RunObservations:
    output_files: list[dict[str, Any]]
    hdf5_outputs: list[str]
    fatal_log_messages: list[str]
    status: str
