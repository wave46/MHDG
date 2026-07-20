"""Typed inputs shared by comparison workflow orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ComparisonInputs:
    run_directory: Path
    case_directory: Path
    tolerances_path: Path
    plan: dict[str, Any]
    case: dict[str, Any]
    workflow: dict[str, Any]


@dataclass(frozen=True)
class ComparisonOverrides:
    candidate: Path | None = None
    reference: Path | None = None
    tolerance_profile: str | None = None
