"""Prepare an isolated MHDG regression run without executing the solver."""

from __future__ import annotations

from pathlib import Path

from preparation.configuration import load_preparation_inputs
from preparation.models import PreparedExecution
from preparation.workflows import prepare_staged_run, prepare_warm_run
from support.errors import BundleError


def prepare_run(
    settings_path: Path,
    case_id: str,
    workflow_id: str,
    layout_id: str,
    case_dir: Path,
    layouts_path: Path,
    run_id: str | None = None,
    validate_bundle: bool = True,
    requested_overrides: dict[str, bool | float | int | str] | None = None,
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
        requested_overrides,
    )
    workflow_kind = inputs.workflow["kind"]
    if workflow_kind == "warm_same_state":
        return prepare_warm_run(inputs)
    if workflow_kind in {"staged_fixed_mesh", "staged_adaptive_mesh"}:
        return prepare_staged_run(inputs)
    raise BundleError(
        f"run preparation does not support workflow kind {workflow_kind}"
    )
