"""Create, resume, and finalize persistent suite summaries."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from support.documents import load_json
from support.errors import BundleError
from support.time import utc_now


def new_summary(
    suite_id: str,
    run_id: str,
    suite: dict[str, Any],
    workflow_ids: list[str],
    comparison_mode: str,
) -> dict[str, Any]:
    """Create the initial running summary for a new suite."""
    summary = {
        "schema_version": 1,
        "started_utc": utc_now(),
        "finished_utc": None,
        "duration_seconds": 0.0,
        "status": "running",
        "suite_id": suite_id,
        "run_id": run_id,
        "description": suite["description"],
        "case_id": suite["case_id"],
        "workflow_ids": workflow_ids,
        "layout_ids": suite["layouts"],
        "comparison_mode": comparison_mode,
        "results": [],
    }
    if len(workflow_ids) == 1:
        summary["workflow_id"] = workflow_ids[0]
    return summary


def resume_summary(
    path: Path,
    suite_id: str,
    run_id: str,
    workflow_ids: list[str],
    layout_ids: list[str],
    comparison_mode: str,
) -> dict[str, Any]:
    """Load an existing compatible summary and mark it running again."""
    summary = load_json(path, "suite summary")
    expected = {
        "suite_id": suite_id,
        "run_id": run_id,
        "workflow_ids": workflow_ids,
        "layout_ids": layout_ids,
        "comparison_mode": comparison_mode,
    }
    mismatched = [key for key, value in expected.items() if summary.get(key) != value]
    if mismatched:
        raise BundleError(f"suite summary does not match: {', '.join(mismatched)}")
    if not isinstance(summary.get("results"), list):
        raise BundleError("suite summary has invalid results")
    _validate_recorded_cells(summary["results"], workflow_ids, layout_ids)
    summary["status"] = "running"
    summary["finished_utc"] = None
    return summary


def completed_cells(summary: dict[str, Any]) -> set[tuple[str, str]]:
    """Return workflow-layout cells already recorded by a summary."""
    return {
        (result["workflow_id"], result["layout_id"])
        for result in summary["results"]
    }


def finalize_summary(summary: dict[str, Any]) -> None:
    """Mark a completed summary passed or failed from its cell results."""
    summary["finished_utc"] = utc_now()
    summary["status"] = (
        "passed"
        if all(result["status"] == "passed" for result in summary["results"])
        else "failed"
    )


def _validate_recorded_cells(
    results: list[Any],
    workflow_ids: list[str],
    layout_ids: list[str],
) -> None:
    planned = {
        (workflow, layout)
        for workflow in workflow_ids
        for layout in layout_ids
    }
    recorded = [
        (result.get("workflow_id"), result.get("layout_id"))
        for result in results
        if isinstance(result, dict)
    ]
    if len(recorded) != len(results) or any(cell not in planned for cell in recorded):
        raise BundleError("suite summary contains an invalid workflow/layout cell")
    if len(set(recorded)) != len(recorded):
        raise BundleError("suite summary contains duplicate workflow/layout cells")
