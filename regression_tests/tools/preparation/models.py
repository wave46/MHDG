"""Typed data passed through run preparation and execution."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PreparedRun:
    path: Path
    command: list[str]
    omp_threads: int
    executable: Path
    runtime_files: dict[str, Path]


@dataclass(frozen=True)
class PreparedStage:
    stage_id: str
    restart_from: str
    run: PreparedRun


@dataclass(frozen=True)
class PreparedStagedRun:
    path: Path
    stages: list[PreparedStage]


PreparedExecution = PreparedRun | PreparedStagedRun


@dataclass(frozen=True)
class PreparationInputs:
    run_directory: Path
    executable: Path
    launcher: Path | None
    layout: dict[str, Any]
    artifacts: dict[str, Path]
    runtime_files: dict[str, Path]
    case: dict[str, Any]
    workflow_id: str
    workflow: dict[str, Any]
    layout_id: str
    bundle_root: Path
    manifest: dict[str, Any]
    requested_overrides: dict[str, bool | float | int | str]
