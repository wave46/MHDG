"""Compare conservative fields and gradients at common interior mesh points."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from support.errors import ComparisonError
from support.paths import require_file
from .compare_common import calculate_error_norms


@dataclass(frozen=True)
class SampledFields:
    equation_names: list[str]
    inside: np.ndarray
    solution: np.ndarray
    gradient: np.ndarray


def reference_sample_points(
    path: Path,
    samples_per_element: int = 4,
) -> np.ndarray:
    """Return deterministic points strictly inside each reference triangle."""
    barycentric = _barycentric_samples(samples_per_element)
    try:
        with h5py.File(path, "r") as handle:
            mesh = handle["mesh"] if "mesh" in handle else handle
            coordinates = _coordinate_rows(np.asarray(mesh["X"]))
            triangles = _triangle_rows(np.asarray(mesh["Tlin"])) - 1
    except (OSError, KeyError) as exc:
        raise ComparisonError(f"cannot read reference mesh: {exc}") from exc

    if triangles.min(initial=0) < 0 or triangles.max(initial=-1) >= len(coordinates):
        raise ComparisonError("reference mesh/Tlin contains invalid node indices")
    vertices = coordinates[triangles]
    return np.einsum("sc,ecd->esd", barycentric, vertices).reshape(-1, 2)


def sample_solution(
    path: Path,
    points: np.ndarray,
    fekete_path: Path,
) -> SampledFields:
    """Interpolate one MHDG solution at the requested points."""
    try:
        from hdg_postprocess.api import load_solution
    except ImportError as exc:
        raise ComparisonError(
            "HDG_postprocess is required for adaptive interpolation"
        ) from exc

    try:
        solution = load_solution(f"{path.parent}{os.sep}", path.stem)
        order = int(solution.mesh.metadata.p_order)
        with h5py.File(fekete_path, "r") as handle:
            nodes = np.asarray(handle[f"P{order}"]).T
        solution.mesh.metadata.reference_element = {"NodesCoord": nodes}
        solution.sample.define_interpolators()
        locator = solution.mesh.geometry.element_locator
        inside = np.fromiter(
            (int(locator(x, y)) >= 0 for x, y in points),
            dtype=bool,
            count=len(points),
        )
        equation_count = int(solution.neq)
        values = np.full((len(points), equation_count), np.nan)
        gradients = np.full((len(points), equation_count, 2), np.nan)
        x_values, y_values = points[inside].T
        for index in range(equation_count):
            values[inside, index] = solution.interpolators.solution[
                index
            ].evaluate_many(x_values, y_values)
            for component in range(2):
                gradients[inside, index, component] = solution.interpolators.gradient[
                    index
                ][component].evaluate_many(x_values, y_values)
        names = _equation_names(solution, equation_count)
    except (
        AttributeError,
        IndexError,
        KeyError,
        OSError,
        RuntimeError,
        TypeError,
        ValueError,
    ) as exc:
        raise ComparisonError(f"cannot sample {path}: {exc}") from exc
    return SampledFields(names, inside, values, gradients)


def _equation_names(solution: Any, equation_count: int) -> list[str]:
    names = [f"equation_{index + 1}" for index in range(equation_count)]
    for raw_name, raw_index in getattr(solution, "_cons_idx", {}).items():
        index = int(raw_index)
        if 0 <= index < equation_count:
            names[index] = (
                raw_name.decode(errors="replace")
                if isinstance(raw_name, bytes)
                else str(raw_name)
            )
    return names


def _barycentric_samples(count: int) -> np.ndarray:
    if count == 1:
        return np.array([[1 / 3, 1 / 3, 1 / 3]])
    if count == 4:
        return np.array(
            [
                [1 / 3, 1 / 3, 1 / 3],
                [0.6, 0.2, 0.2],
                [0.2, 0.6, 0.2],
                [0.2, 0.2, 0.6],
            ]
        )
    raise ComparisonError("samples per element must be 1 or 4")


def _coordinate_rows(array: np.ndarray) -> np.ndarray:
    if array.ndim != 2 or 2 not in array.shape:
        raise ComparisonError(
            "reference mesh/X must be a two-dimensional coordinate array"
        )
    return array.T if array.shape[0] == 2 else array


def _triangle_rows(array: np.ndarray) -> np.ndarray:
    if array.ndim != 2 or 3 not in array.shape:
        raise ComparisonError(
            "reference mesh/Tlin must contain three-node triangles"
        )
    return (array.T if array.shape[0] == 3 else array).astype(np.int64)


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
        "schema_version": 2,
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
