"""Execute and optionally compare one workflow-layout suite cell."""

from __future__ import annotations

from typing import Any

from comparison.workflow import compare_completed_run
from prepare_run import prepare_run
from run_case import execute_prepared
from suite.models import SuiteRunInputs
from support.errors import HarnessError


def run_cell(
    inputs: SuiteRunInputs,
    workflow_id: str,
    layout_id: str,
) -> dict[str, Any]:
    """Prepare, execute, and optionally compare one suite cell."""
    result = _empty_result(workflow_id, layout_id)
    try:
        prepared = prepare_run(
            inputs.settings_path,
            inputs.case_id,
            workflow_id,
            layout_id,
            inputs.case_directory,
            inputs.layouts_path,
            inputs.run_id,
            validate_bundle=False,
        )
        result["run_directory"] = str(prepared.path)
        run = execute_prepared(prepared, inputs.settings)
        result["run_status"] = run.status
        result["duration_seconds"] = run.duration_seconds
        result["status"] = run.status
        if run.status != "completed":
            result["failures"] = [f"solver run status is {run.status}"]
            return result
        if not inputs.compare:
            result["comparison_status"] = "not_run"
            result["status"] = "passed"
            return result

        policy, report_path, comparison = compare_completed_run(
            run.path,
            inputs.case_directory,
            inputs.tolerances_path,
        )
        result["comparison_policy"] = policy
        result["comparison_status"] = comparison["status"]
        result["comparison_report"] = str(report_path)
        result["failures"] = comparison["failures"]
        result["status"] = (
            "passed"
            if comparison["status"] == "passed"
            else "comparison_failed"
        )
    except HarnessError as exc:
        result["failures"] = [str(exc)]
    return result


def _empty_result(workflow_id: str, layout_id: str) -> dict[str, Any]:
    return {
        "workflow_id": workflow_id,
        "layout_id": layout_id,
        "status": "error",
        "run_status": None,
        "comparison_policy": None,
        "comparison_status": None,
        "duration_seconds": None,
        "run_directory": None,
        "comparison_report": None,
        "failures": [],
    }
