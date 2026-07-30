"""Select the comparison policy declared by a completed workflow run."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.adaptive.run import compare_adaptive_run
from comparison.fixed.run import compare_fixed_run
from comparison.inputs import (
    ComparisonInputs,
    ComparisonOverrides,
    load_comparison_inputs,
)
from comparison.matrix import compare_reference_matrix, load_plan_matrix
from support.errors import ComparisonError


def compare_completed_run(
    run_directory: Path,
    case_dir: Path,
    tolerances_path: Path,
    candidate_override: Path | None = None,
    reference_override: Path | None = None,
    tolerance_profile_override: str | None = None,
    comparison_policy_override: str | None = None,
    report_override: Path | None = None,
) -> tuple[str, Path, dict[str, Any]]:
    """Dispatch one completed run to its final-state or staged comparison."""
    inputs = load_comparison_inputs(run_directory, case_dir, tolerances_path)
    overrides = ComparisonOverrides(
        candidate=candidate_override,
        reference=reference_override,
        tolerance_profile=tolerance_profile_override,
    )

    explicit_final_state = any(
        value is not None
        for value in (
            overrides.candidate,
            overrides.reference,
            overrides.tolerance_profile,
            comparison_policy_override,
        )
    )
    if inputs.workflow.get("stages") and not explicit_final_state:
        matrix = load_plan_matrix(
            inputs.plan,
            inputs.case,
            inputs.case_directory.parent / "schemas",
        )
        if matrix is not None:
            path, report = compare_reference_matrix(
                inputs,
                matrix,
                report_override=report_override,
            )
            return "reference_matrix", path, report

    return _compare_final_state(
        inputs,
        overrides,
        report_override,
        comparison_policy_override,
    )


def _compare_final_state(
    inputs: ComparisonInputs,
    overrides: ComparisonOverrides,
    report_path: Path | None,
    policy_override: str | None = None,
) -> tuple[str, Path, dict[str, Any]]:
    policy = policy_override or inputs.workflow.get("comparison_policy")
    if policy == "fixed_hdf5":
        path, report = compare_fixed_run(
            inputs,
            overrides,
            report_path,
        )
    elif policy == "mesh_independent":
        path, report = compare_adaptive_run(
            inputs,
            overrides,
            report_path,
        )
    else:
        raise ComparisonError(f"unsupported comparison policy: {policy}")
    return policy, path, report
