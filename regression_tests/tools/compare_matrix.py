"""Dispatch completed runs to final-state or staged golden comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from check_bundle import load_case_definition
from compare_adaptive import compare_adaptive_run
from compare_run import compare_run
from reference_matrix import ReferenceMatrix, load_reference_matrix
from support.documents import load_json, write_json_atomic
from support.errors import BundleError, ComparisonError
from support.paths import recorded_directory, recorded_file, require_directory
from support.time import utc_now


def compare_completed_run(
    run_directory: Path,
    case_dir: Path,
    tolerances_path: Path,
    candidate_override: Path | None = None,
    reference_override: Path | None = None,
    profile_override: str | None = None,
    report_override: Path | None = None,
) -> tuple[str, Path, dict[str, Any]]:
    """Use a matching stage matrix when present, otherwise compare the final state."""
    run_directory = require_directory(run_directory, "run")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    case = load_case_definition(plan["case_id"], case_dir)
    workflow = case["workflows"].get(plan.get("workflow_id"))
    if workflow is None:
        raise ComparisonError("run plan refers to an unknown workflow")

    if workflow.get("stages"):
        matrix = _matrix_for_plan(plan, case, case_dir.parent / "schemas")
        if matrix is not None:
            if any(
                value is not None
                for value in (
                    candidate_override,
                    reference_override,
                    profile_override,
                )
            ):
                raise ComparisonError(
                    "candidate, reference, and profile overrides are not supported "
                    "for staged comparisons"
                )
            path, report = _compare_stages(
                run_directory,
                plan,
                case,
                workflow,
                matrix,
                case_dir,
                tolerances_path,
                report_override,
            )
            return "reference_matrix", path, report

    policy = workflow.get("comparison_policy")
    if policy == "fixed_hdf5":
        path, report = compare_run(
            run_directory,
            case_dir,
            tolerances_path,
            candidate_override=candidate_override,
            reference_override=reference_override,
            profile_override=profile_override,
            report_override=report_override,
        )
    elif policy == "mesh_independent":
        path, report = compare_adaptive_run(
            run_directory,
            case_dir,
            tolerances_path,
            report_path=report_override,
            candidate_override=candidate_override,
            reference_override=reference_override,
            profile_override=profile_override,
        )
    else:
        raise ComparisonError(f"unsupported comparison policy: {policy}")
    return policy, path, report


def _compare_stages(
    run_directory: Path,
    plan: dict[str, Any],
    case: dict[str, Any],
    workflow: dict[str, Any],
    matrix: ReferenceMatrix,
    case_dir: Path,
    tolerances_path: Path,
    report_override: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")

    expected = [stage["stage_id"] for stage in workflow["stages"]]
    plan_stages = [stage.get("stage_id") for stage in plan.get("stages", [])]
    recorded = metadata.get("stages", [])
    recorded_stages = [stage.get("stage_id") for stage in recorded]
    if plan_stages != expected or recorded_stages != expected:
        raise ComparisonError("recorded stages differ from the tracked workflow")

    profile = workflow.get("stage_tolerance_profile")
    policy = workflow["comparison_policy"]
    stage_reports = []
    failures = []
    for stage in recorded:
        stage_id = stage["stage_id"]
        if stage.get("status") != "completed":
            raise ComparisonError(f"stage {stage_id} did not complete")
        stage_directory = recorded_directory(
            stage.get("run_directory"), "stage run"
        )
        candidate = recorded_file(stage.get("selected_hdf5"), "stage result")
        reference = matrix.reference_for(
            plan["workflow_id"], plan["layout_id"], stage_id
        )
        report_path = stage_directory / "comparison.json"
        if policy == "fixed_hdf5":
            report_path, report = compare_run(
                stage_directory,
                case_dir,
                tolerances_path,
                candidate_override=candidate,
                reference_override=reference,
                profile_override=profile,
                report_override=report_path,
            )
        elif policy == "mesh_independent":
            report_path, report = compare_adaptive_run(
                stage_directory,
                case_dir,
                tolerances_path,
                report_path=report_path,
                candidate_override=candidate,
                reference_override=reference,
                profile_override=profile,
            )
        else:
            raise ComparisonError(f"unsupported comparison policy: {policy}")

        stage_reports.append(
            {
                "stage_id": stage_id,
                "status": report["status"],
                "comparison_report": str(report_path),
                "reference": str(reference),
                "candidate": str(candidate),
                "failures": report["failures"],
            }
        )
        if report["status"] != "passed":
            stage_failures = report["failures"] or ["comparison failed"]
            failures.extend(f"{stage_id}: {item}" for item in stage_failures)
            break

    total = len(expected)
    aggregate = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "passed" if not failures and len(stage_reports) == total else "failed",
        "run_directory": str(run_directory),
        "case_id": case["case_id"],
        "workflow_id": plan["workflow_id"],
        "layout_id": plan["layout_id"],
        "comparison_policy": "reference_matrix",
        "checked_stage_count": len(stage_reports),
        "total_stage_count": total,
        "first_failed_stage": stage_reports[-1]["stage_id"] if failures else None,
        "stages": stage_reports,
        "failures": failures,
    }
    output = (report_override or run_directory / "matrix_comparison.json").resolve()
    reserved = {
        Path(stage[name]).resolve()
        for stage in stage_reports
        for name in ("candidate", "reference", "comparison_report")
    }
    if output in reserved:
        raise ComparisonError("matrix report cannot replace a stage artifact")
    write_json_atomic(output, aggregate, "matrix report")
    return output, aggregate


def _matrix_for_plan(
    plan: dict[str, Any],
    case: dict[str, Any],
    schema_dir: Path,
) -> ReferenceMatrix | None:
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
