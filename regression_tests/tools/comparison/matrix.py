"""Staged reference-matrix comparison and report assembly."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.adaptive.run import compare_adaptive_run
from comparison.fixed.run import compare_fixed_run
from comparison.inputs import (
    ComparisonInputs,
    ComparisonOverrides,
    load_stage_inputs,
)
from reference_matrix import ReferenceMatrix, load_reference_matrix
from support.documents import write_json_atomic
from support.errors import BundleError, ComparisonError
from support.paths import recorded_directory, recorded_file
from support.time import utc_now


def compare_reference_matrix(
    context: ComparisonInputs,
    matrix: ReferenceMatrix,
    report_override: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Compare recorded workflow stages until completion or first failure."""
    stages = _validated_stage_records(
        context.plan,
        context.workflow,
        context.metadata,
    )

    stage_reports = []
    failures = []
    for stage in stages:
        stage_report = _compare_stage(stage, context, matrix)
        stage_reports.append(stage_report)
        failures = _stage_failures(stage_report)
        if failures:
            break

    report = _matrix_report(
        context,
        len(stages),
        stage_reports,
        failures,
    )
    output = _save_matrix_report(
        context.run_directory, report_override, stage_reports, report
    )
    return output, report


def load_plan_matrix(
    plan: dict[str, Any],
    case: dict[str, Any],
    schema_dir: Path,
) -> ReferenceMatrix | None:
    """Load the run's reference matrix and confirm its bundle identity."""
    bundle = plan.get("bundle")
    if not isinstance(bundle, dict) or not isinstance(bundle.get("root"), str):
        raise BundleError("run plan has no bundle identity")
    matrix = load_reference_matrix(Path(bundle["root"]), case, schema_dir)
    if matrix is not None and (
        bundle.get("bundle_id") != matrix.bundle_id
        or bundle.get("bundle_version") != matrix.bundle_version
    ):
        raise BundleError("run and golden-matrix bundle identities differ")
    return matrix


def _validated_stage_records(
    plan: dict[str, Any],
    workflow: dict[str, Any],
    metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")

    expected = [stage["stage_id"] for stage in workflow["stages"]]
    plan_stages = [stage.get("stage_id") for stage in plan.get("stages", [])]
    recorded = metadata.get("stages", [])
    recorded_stages = [stage.get("stage_id") for stage in recorded]
    if plan_stages != expected or recorded_stages != expected:
        raise ComparisonError("recorded stages differ from the tracked workflow")
    return recorded


def _compare_stage(
    stage: dict[str, Any],
    context: ComparisonInputs,
    matrix: ReferenceMatrix,
) -> dict[str, Any]:
    stage_id = stage["stage_id"]
    if stage.get("status") != "completed":
        raise ComparisonError(f"stage {stage_id} did not complete")

    stage_directory = recorded_directory(stage.get("run_directory"), "stage run")
    candidate = recorded_file(stage.get("selected_hdf5"), "stage result")
    reference = matrix.reference_for(
        context.plan["workflow_id"], context.plan["layout_id"], stage_id
    )
    stage_inputs = load_stage_inputs(context, stage_directory)
    overrides = ComparisonOverrides(
        candidate=candidate,
        reference=reference,
        tolerance_profile=context.workflow.get("stage_tolerance_profile"),
    )
    report_path, report = _run_stage_comparison(
        context.workflow["comparison_policy"],
        stage_inputs,
        overrides,
    )
    return {
        "stage_id": stage_id,
        "status": report["status"],
        "comparison_report": str(report_path),
        "reference": str(reference),
        "candidate": str(candidate),
        "failures": report["failures"],
    }


def _run_stage_comparison(
    policy: str,
    inputs: ComparisonInputs,
    overrides: ComparisonOverrides,
) -> tuple[Path, dict[str, Any]]:
    report_path = inputs.run_directory / "comparison.json"
    if policy == "fixed_hdf5":
        return compare_fixed_run(
            inputs,
            overrides=overrides,
            report_path=report_path,
        )
    if policy == "mesh_independent":
        return compare_adaptive_run(
            inputs,
            overrides=overrides,
            report_path=report_path,
        )
    raise ComparisonError(f"unsupported comparison policy: {policy}")


def _stage_failures(stage_report: dict[str, Any]) -> list[str]:
    if stage_report["status"] == "passed":
        return []
    failures = stage_report["failures"] or ["comparison failed"]
    return [f"{stage_report['stage_id']}: {failure}" for failure in failures]


def _matrix_report(
    context: ComparisonInputs,
    total_stage_count: int,
    stage_reports: list[dict[str, Any]],
    failures: list[str],
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": (
            "passed"
            if not failures and len(stage_reports) == total_stage_count
            else "failed"
        ),
        "run_directory": str(context.run_directory),
        "case_id": context.case["case_id"],
        "workflow_id": context.plan["workflow_id"],
        "layout_id": context.plan["layout_id"],
        "comparison_policy": "reference_matrix",
        "checked_stage_count": len(stage_reports),
        "total_stage_count": total_stage_count,
        "first_failed_stage": stage_reports[-1]["stage_id"] if failures else None,
        "stages": stage_reports,
        "failures": failures,
    }


def _save_matrix_report(
    run_directory: Path,
    report_override: Path | None,
    stage_reports: list[dict[str, Any]],
    report: dict[str, Any],
) -> Path:
    output = (report_override or run_directory / "matrix_comparison.json").resolve()
    reserved = {
        Path(stage[name]).resolve()
        for stage in stage_reports
        for name in ("candidate", "reference", "comparison_report")
    }
    if output in reserved:
        raise ComparisonError("matrix report cannot replace a stage artifact")
    write_json_atomic(output, report, "matrix report")
    return output
