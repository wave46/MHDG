"""Typed inputs shared while executing suite cells."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SuiteRunInputs:
    settings_path: Path
    settings: dict[str, str]
    case_id: str
    run_id: str
    case_directory: Path
    layouts_path: Path
    tolerances_path: Path
    compare: bool
    parameter_overrides: dict[str, bool | float | int | str]
