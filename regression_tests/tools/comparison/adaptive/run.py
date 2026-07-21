"""Compare adaptive MHDG solutions on different meshes."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.adaptive.fields import compare_sampled_fields
from comparison.adaptive.sampling import (
    SampledFields,
    reference_sample_points,
    sample_solution,
)
from comparison.inputs import ComparisonInputs, ComparisonOverrides
from comparison.shared.convergence import NewtonConvergence, read_newton_convergence
from comparison.shared.outputs import resolve_run_file, select_candidate
from comparison.shared.report import (
    comparison_report_fields,
    merge_comparison_failures,
    save_comparison_report,
)
from comparison.shared.tolerances import load_adaptive_tolerances
from support.errors import ComparisonError
from support.paths import recorded_file, require_file


def compare_adaptive_files(
    reference_path: Path,
    candidate_path: Path,
    fekete_path: Path,
    samples_per_element: int = 4,
    tolerances: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Interpolate two files at common points and compare their fields."""
    reference_path = require_file(reference_path, "reference")
    candidate_path = require_file(candidate_path, "candidate")
    fekete_path = require_file(fekete_path, "Fekete-node")
    points = reference_sample_points(reference_path, samples_per_element)
    reference = sample_solution(reference_path, points, fekete_path)
    candidate = sample_solution(candidate_path, points, fekete_path)
    report = compare_sampled_fields(
        reference,
        candidate,
        samples_per_element,
        tolerances,
    )
    report["files"] = {
        "reference": str(reference_path),
        "candidate": str(candidate_path),
        "fekete_nodes": str(fekete_path),
    }
    return report


def compare_adaptive_run(
    inputs: ComparisonInputs,
    overrides: ComparisonOverrides,
    report_path: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Compare one completed adaptive run with its bundled reference."""
    if inputs.workflow.get("comparison_policy") != "mesh_independent":
        raise ComparisonError("run does not define mesh-independent comparison")
    profile_id, tolerances = load_adaptive_tolerances(
        inputs.tolerances_path,
        inputs.workflow,
        overrides.tolerance_profile,
    )
    candidate = select_candidate(
        inputs.run_directory,
        inputs.metadata,
        overrides.candidate,
    )
    reference = resolve_run_file(
        inputs.run_directory,
        overrides.reference,
        "inputs/reference.h5",
        "reference",
    )
    fekete = _fekete_nodes(inputs.metadata)

    field_report = compare_adaptive_files(
        reference,
        candidate,
        fekete,
        tolerances["samples_per_element"],
        tolerances,
    )
    convergence = read_newton_convergence(
        inputs.run_directory / "stdout.log",
        tolerances["newton_error_max"],
    )
    failures = merge_comparison_failures(field_report["failures"], convergence)
    report = _comparison_report(
        inputs,
        profile_id,
        tolerances,
        convergence,
        field_report,
        failures,
    )
    output = save_comparison_report(
        inputs.run_directory,
        report_path,
        report,
        (candidate, reference, fekete),
    )
    return output, report


def _comparison_report(
    inputs: ComparisonInputs,
    tolerance_profile_id: str,
    tolerances: dict[str, Any],
    convergence: NewtonConvergence,
    field_report: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    report = {
        **field_report,
        **comparison_report_fields(
            inputs,
            tolerance_profile_id,
            tolerances,
            convergence,
            failures,
        ),
    }
    report.pop("tolerances", None)
    return report


def _fekete_nodes(metadata: dict[str, Any]) -> Path:
    runtime = metadata.get("runtime_files", {}).get(
        "positionFeketeNodesTri2D.h5", {}
    )
    return recorded_file(runtime.get("path"), "Fekete-node")
