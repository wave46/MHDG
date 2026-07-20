"""Compare fields already sampled at common physical points."""

from __future__ import annotations

from typing import Any

import numpy as np

from comparison.metrics import calculate_error_norms
from comparison.sampling import SampledFields
from support.errors import ComparisonError


def compare_sampled_fields(
    reference: SampledFields,
    candidate: SampledFields,
    samples_per_element: int,
    tolerances: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare sampled conservative state and gradient arrays."""
    _validate_sample_shapes(reference, candidate)
    failures = []
    common = reference.inside & candidate.inside
    point_count = int(common.size)
    common_count = int(np.count_nonzero(common))
    coverage = common_count / point_count if point_count else 0.0
    if common_count == 0:
        failures.append("the meshes have no common sampling points")

    if reference.equation_names != candidate.equation_names:
        failures.append("conservative variable names differ")
    names = reference.equation_names
    datasets = {}
    arrays = {
        "solution": (reference.solution, candidate.solution),
        "gradient_x": (reference.gradient[:, :, 0], candidate.gradient[:, :, 0]),
        "gradient_y": (reference.gradient[:, :, 1], candidate.gradient[:, :, 1]),
    }
    for dataset_name, (first, second) in arrays.items():
        equations = {}
        for index, name in enumerate(names):
            limits = None
            if tolerances is not None:
                limits = tolerances[
                    "solution" if dataset_name == "solution" else "gradient"
                ]
            metrics = _numeric_metrics(
                first[common, index],
                second[common, index],
                limits,
            )
            equations[name] = metrics
            if not metrics["finite"]:
                failures.append(f"{dataset_name}/{name} contains non-finite values")
            elif metrics["passed"] is False:
                failures.append(f"{dataset_name}/{name} exceeds tolerance")
        datasets[dataset_name] = {"equations": equations}

    if tolerances is not None and coverage < tolerances["minimum_point_coverage"]:
        failures.append("common point coverage is below tolerance")
    status = "failed" if failures else "passed" if tolerances else "characterized"
    return {
        "schema_version": 1,
        "status": status,
        "sampling": {
            "method": "reference_triangle_interior",
            "samples_per_element": samples_per_element,
            "point_count": point_count,
            "reference_points": int(np.count_nonzero(reference.inside)),
            "candidate_points": int(np.count_nonzero(candidate.inside)),
            "common_points": common_count,
            "common_coverage": coverage,
        },
        "equation_count": len(names),
        "equation_names": names,
        "datasets": datasets,
        "tolerances": tolerances,
        "failures": failures,
    }


def _numeric_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    limits: dict[str, float] | None,
) -> dict[str, Any]:
    norms = calculate_error_norms(reference, candidate)
    passed = None
    if limits is not None:
        passed = bool(
            norms.available
            and norms.relative_l2 <= limits["relative_l2_max"]
            and norms.normalized_linf <= limits["normalized_linf_max"]
        )
    return {
        "finite": norms.finite,
        "relative_l2": norms.relative_l2,
        "normalized_linf": norms.normalized_linf,
        "passed": passed,
    }


def _validate_sample_shapes(
    reference: SampledFields,
    candidate: SampledFields,
) -> None:
    expected = reference.solution.shape
    if expected != candidate.solution.shape or len(expected) != 2:
        raise ComparisonError("sampled solution shapes differ")
    gradient_shape = (*expected, 2)
    if (
        reference.gradient.shape != gradient_shape
        or candidate.gradient.shape != gradient_shape
    ):
        raise ComparisonError("sampled gradient shapes differ")
    if reference.inside.shape != (expected[0],) or candidate.inside.shape != (
        expected[0],
    ):
        raise ComparisonError("sampling masks have incompatible shapes")
    if len(reference.equation_names) != expected[1] or len(
        candidate.equation_names
    ) != expected[1]:
        raise ComparisonError("equation names do not match sampled fields")
