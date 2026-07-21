"""Compare completed suite runs directly across execution layouts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.shared.outputs import select_candidate
from comparison.workflow import compare_completed_run
from support.documents import load_json
from support.errors import ComparisonError, HarnessError
from support.paths import recorded_directory


def compare_layout_pairs(
    summary: dict[str, Any],
    case_directory: Path,
    tolerances_path: Path,
) -> list[dict[str, Any]]:
    """Compare every declared candidate with its same-workflow baseline."""
    return [
        _compare_pair(
            summary,
            workflow_id,
            pair,
            case_directory,
            tolerances_path,
        )
        for workflow_id in summary["workflow_ids"]
        for pair in summary["layout_comparisons"]
    ]


def _compare_pair(
    summary: dict[str, Any],
    workflow_id: str,
    pair: dict[str, str],
    case_directory: Path,
    tolerances_path: Path,
) -> dict[str, Any]:
    baseline_layout = pair["baseline"]
    candidate_layout = pair["candidate"]
    result = {
        "workflow_id": workflow_id,
        "baseline_layout_id": baseline_layout,
        "candidate_layout_id": candidate_layout,
        "baseline_run_directory": None,
        "candidate_run_directory": None,
        "baseline_output": None,
        "comparison_policy": None,
        "comparison_report": None,
        "status": "failed",
        "failures": [],
    }
    try:
        baseline = _completed_run(summary, workflow_id, baseline_layout)
        candidate = _completed_run(summary, workflow_id, candidate_layout)
        result["baseline_run_directory"] = str(baseline)
        result["candidate_run_directory"] = str(candidate)
        reference = _selected_output(baseline)
        result["baseline_output"] = str(reference)
        policy, report_path, report = compare_completed_run(
            candidate,
            case_directory,
            tolerances_path,
            reference_override=reference,
            tolerance_profile_override=summary["tolerance_profile"],
        )
        result["comparison_policy"] = policy
        result["comparison_report"] = str(report_path)
        result["status"] = report["status"]
        result["failures"] = report["failures"]
    except HarnessError as exc:
        result["failures"] = [str(exc)]
    return result


def _completed_run(
    summary: dict[str, Any],
    workflow_id: str,
    layout_id: str,
) -> Path:
    matches = [
        result
        for result in summary["results"]
        if result.get("workflow_id") == workflow_id
        and result.get("layout_id") == layout_id
    ]
    if len(matches) != 1:
        raise ComparisonError(
            f"suite has no unique result for {workflow_id}/{layout_id}"
        )
    result = matches[0]
    if result.get("run_status") != "completed":
        raise ComparisonError(
            f"suite run did not complete: {workflow_id}/{layout_id}"
        )
    return recorded_directory(result.get("run_directory"), "suite run")


def _selected_output(run_directory: Path) -> Path:
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    return select_candidate(run_directory, metadata)
