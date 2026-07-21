"""Typed inputs shared by comparison workflow orchestration."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition
from support.documents import load_json
from support.errors import ComparisonError
from support.paths import require_directory


@dataclass(frozen=True)
class ComparisonInputs:
    run_directory: Path
    case_directory: Path
    tolerances_path: Path
    plan: dict[str, Any]
    metadata: dict[str, Any]
    case: dict[str, Any]
    workflow: dict[str, Any]


@dataclass(frozen=True)
class ComparisonOverrides:
    candidate: Path | None = None
    reference: Path | None = None
    tolerance_profile: str | None = None


def load_comparison_inputs(
    run_directory: Path,
    case_directory: Path,
    tolerances_path: Path,
) -> ComparisonInputs:
    """Load and validate the documents needed to compare a completed run."""
    run_directory = require_directory(run_directory, "run")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    metadata = _load_completed_metadata(run_directory)
    case = load_case_definition(plan["case_id"], case_directory)
    workflow = case["workflows"].get(plan.get("workflow_id"))
    if workflow is None:
        raise ComparisonError("run plan refers to an unknown workflow")
    return ComparisonInputs(
        run_directory=run_directory,
        case_directory=case_directory,
        tolerances_path=tolerances_path,
        plan=plan,
        metadata=metadata,
        case=case,
        workflow=workflow,
    )


def load_stage_inputs(
    parent: ComparisonInputs,
    run_directory: Path,
) -> ComparisonInputs:
    """Use parent workflow inputs with one completed stage's run data."""
    run_directory = require_directory(run_directory, "stage run")
    return replace(
        parent,
        run_directory=run_directory,
        metadata=_load_completed_metadata(run_directory),
    )


def _load_completed_metadata(run_directory: Path) -> dict[str, Any]:
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")
    return metadata
