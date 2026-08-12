"""Verify recorded suite results without rerunning the solver."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition
from comparison.shared.convergence import (
    effective_newton_maximum,
    read_newton_convergence,
)
from comparison.shared.tolerances import (
    load_adaptive_tolerances,
    load_fixed_tolerances,
)
from comparison.workflow import compare_completed_run
from suite.pairs import compare_layout_pairs
from support.documents import load_json, write_json_atomic
from support.errors import BundleError, HarnessError
from support.paths import require_file
from support.time import utc_now


def verify_suite(
    suite_summary_path: Path,
    case_directory: Path,
    tolerances_path: Path,
    *,
    include_layout_pairs: bool = True,
) -> tuple[Path, dict[str, Any]]:
    """Compare all recorded suite runs without executing the solver."""
    suite_summary_path = require_file(suite_summary_path, "suite summary")
    source = load_json(suite_summary_path, "suite summary")
    _validate_source_summary(source)
    case = load_case_definition(source["case_id"], case_directory)

    output_path = suite_summary_path.parent / "verification_summary.json"
    summary = _new_verification_summary(source, suite_summary_path)
    write_json_atomic(output_path, summary, "verification summary")
    reference_comparisons = source.get(
        "reference_comparisons",
        not source.get("layout_comparisons"),
    )
    if reference_comparisons:
        for source_result in source["results"]:
            summary["results"].append(
                _verify_result(
                    source_result,
                    case,
                    case_directory,
                    tolerances_path,
                )
            )
            write_json_atomic(output_path, summary, "verification summary")

    if include_layout_pairs and source.get("layout_comparisons"):
        summary["comparisons"] = compare_layout_pairs(
            source,
            case_directory,
            tolerances_path,
        )

    checks = [*summary["results"], *summary.get("comparisons", [])]
    summary["status"] = (
        "passed"
        if all(result["status"] == "passed" for result in checks)
        else "failed"
    )
    summary["finished_utc"] = utc_now()
    write_json_atomic(output_path, summary, "verification summary")
    return output_path, summary


def _new_verification_summary(
    source: dict[str, Any],
    source_path: Path,
) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "created_utc": utc_now(),
        "status": "running",
        "source_summary": str(source_path),
        "suite_id": source["suite_id"],
        "run_id": source["run_id"],
        "case_id": source["case_id"],
        "results": [],
    }


def _verify_result(
    source: dict[str, Any],
    case: dict[str, Any],
    case_directory: Path,
    tolerances_path: Path,
) -> dict[str, Any]:
    workflow_id = source.get("workflow_id")
    result = {
        "workflow_id": workflow_id,
        "layout_id": source.get("layout_id"),
        "run_directory": source.get("run_directory"),
        "run_status": source.get("run_status"),
        "comparison_policy": None,
        "comparison_report": None,
        "convergence_status": None,
        "status": "failed",
        "failures": [],
    }
    if source.get("run_status") != "completed":
        result["failures"] = ["solver run did not complete"]
        return result
    if workflow_id not in case["workflows"]:
        result["failures"] = [f"unknown workflow: {workflow_id}"]
        return result

    try:
        run_directory = Path(source["run_directory"])
        policy, report_path, report = compare_completed_run(
            run_directory,
            case_directory,
            tolerances_path,
        )
        result["comparison_policy"] = policy
        result["comparison_report"] = str(report_path)
        result["convergence_status"] = (
            "passed"
            if _producer_converged(
                source,
                case["workflows"][workflow_id],
                policy,
                report,
                tolerances_path,
            )
            else "failed"
        )
        result["failures"] = report["failures"]
        result["status"] = report["status"]
    except (HarnessError, TypeError) as exc:
        result["failures"] = [str(exc)]
    return result


def _producer_converged(
    source: dict[str, Any],
    workflow: dict[str, Any],
    policy: str,
    report: dict[str, Any],
    tolerances_path: Path,
) -> bool:
    """Check convergence independently of old-reference field agreement."""
    convergence = report.get("convergence")
    if isinstance(convergence, dict):
        return convergence.get("passed") is True
    if policy != "reference_matrix" or not workflow.get("stages"):
        return False

    run_directory = Path(source["run_directory"])
    metadata = load_json(
        run_directory / "run_metadata.json",
        "staged run metadata",
    )
    records = metadata.get("stages")
    definitions = workflow["stages"]
    if not isinstance(records, list) or len(records) != len(definitions):
        return False
    if any(
        record.get("stage_id") != definition["stage_id"]
        or record.get("status") != "completed"
        for record, definition in zip(records, definitions)
    ):
        return False

    maximum = _stage_newton_maximum(
        workflow,
        source["layout_id"],
        tolerances_path,
    )
    for record, definition in zip(records, definitions):
        effective_maximum = effective_newton_maximum(
            maximum,
            definition["newton_check"],
        )
        convergence = read_newton_convergence(
            Path(record["run_directory"]) / "stdout.log",
            effective_maximum,
        )
        if not convergence.passed:
            return False
    return True


def _stage_newton_maximum(
    workflow: dict[str, Any],
    layout_id: str,
    tolerances_path: Path,
) -> float | None:
    profile = workflow.get("stage_tolerance_profile")
    if workflow.get("comparison_policy") == "fixed_hdf5":
        _, tolerances = load_fixed_tolerances(
            tolerances_path,
            workflow,
            layout_id,
            profile,
        )
    else:
        _, tolerances = load_adaptive_tolerances(
            tolerances_path,
            workflow,
            profile,
        )
    return tolerances["newton_error_max"]


def _validate_source_summary(summary: dict[str, Any]) -> None:
    required = {"schema_version", "suite_id", "run_id", "case_id", "results"}
    missing = sorted(required - summary.keys())
    if missing:
        raise BundleError(f"suite summary is missing: {', '.join(missing)}")
    if summary["schema_version"] != 2:
        raise BundleError("suite summary has unsupported schema version")
    if not isinstance(summary["results"], list) or not summary["results"]:
        raise BundleError("suite summary contains no results")
    result_fields = {"workflow_id", "layout_id", "run_directory", "run_status"}
    for result in summary["results"]:
        if not isinstance(result, dict) or not result_fields <= result.keys():
            raise BundleError("suite summary contains an invalid result")
    if "layout_comparisons" in summary:
        pair_fields = {"workflow_ids", "tolerance_profile"}
        if not pair_fields <= summary.keys():
            raise BundleError("paired suite summary is incomplete")
        pairs = summary["layout_comparisons"]
        if not isinstance(pairs, list) or not pairs:
            raise BundleError("paired suite summary has no layout comparisons")
        if any(
            not isinstance(pair, dict)
            or set(pair) != {"baseline", "candidate"}
            for pair in pairs
        ):
            raise BundleError("paired suite summary has an invalid comparison")
