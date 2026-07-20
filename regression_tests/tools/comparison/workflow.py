"""Select the comparison policy declared by a completed workflow run."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from check_bundle import load_case_definition
from comparison.adaptive.run import compare_adaptive_run
from comparison.fixed.run import compare_run
from comparison.inputs import ComparisonInputs, ComparisonOverrides
from comparison.matrix import compare_reference_matrix, load_plan_matrix
from support.documents import load_json
from support.errors import ComparisonError
from support.paths import require_directory


def compare_completed_run(
    run_directory: Path,
    case_dir: Path,
    tolerances_path: Path,
    candidate_override: Path | None = None,
    reference_override: Path | None = None,
    tolerance_profile_override: str | None = None,
    report_override: Path | None = None,
) -> tuple[str, Path, dict[str, Any]]:
    """Dispatch one completed run to its final-state or staged comparison."""
    context = _load_context(run_directory, case_dir, tolerances_path)
    overrides = ComparisonOverrides(
        candidate=candidate_override,
        reference=reference_override,
        tolerance_profile=tolerance_profile_override,
    )

    if context.workflow.get("stages"):
        matrix = load_plan_matrix(
            context.plan,
            context.case,
            context.case_directory.parent / "schemas",
        )
        if matrix is not None:
            _reject_staged_overrides(overrides)
            path, report = compare_reference_matrix(
                context,
                matrix,
                report_override=report_override,
            )
            return "reference_matrix", path, report

    return _compare_final_state(context, overrides, report_override)


def _load_context(
    run_directory: Path,
    case_dir: Path,
    tolerances_path: Path,
) -> ComparisonInputs:
    run_directory = require_directory(run_directory, "run")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    case = load_case_definition(plan["case_id"], case_dir)
    workflow = case["workflows"].get(plan.get("workflow_id"))
    if workflow is None:
        raise ComparisonError("run plan refers to an unknown workflow")
    return ComparisonInputs(
        run_directory=run_directory,
        case_directory=case_dir,
        tolerances_path=tolerances_path,
        plan=plan,
        case=case,
        workflow=workflow,
    )


def _compare_final_state(
    context: ComparisonInputs,
    overrides: ComparisonOverrides,
    report_path: Path | None,
) -> tuple[str, Path, dict[str, Any]]:
    policy = context.workflow.get("comparison_policy")
    if policy == "fixed_hdf5":
        path, report = compare_run(
            context.run_directory,
            context.case_directory,
            context.tolerances_path,
            candidate_override=overrides.candidate,
            reference_override=overrides.reference,
            tolerance_profile_override=overrides.tolerance_profile,
            report_override=report_path,
        )
    elif policy == "mesh_independent":
        path, report = compare_adaptive_run(
            context.run_directory,
            context.case_directory,
            context.tolerances_path,
            report_path=report_path,
            candidate_override=overrides.candidate,
            reference_override=overrides.reference,
            tolerance_profile_override=overrides.tolerance_profile,
        )
    else:
        raise ComparisonError(f"unsupported comparison policy: {policy}")
    return policy, path, report


def _reject_staged_overrides(overrides: ComparisonOverrides) -> None:
    if any(
        value is not None
        for value in (
            overrides.candidate,
            overrides.reference,
            overrides.tolerance_profile,
        )
    ):
        raise ComparisonError(
            "candidate, reference, and tolerance-profile overrides are not supported "
            "for staged comparisons"
        )
