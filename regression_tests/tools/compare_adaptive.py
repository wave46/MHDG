#!/usr/bin/env python3
"""Characterize two MHDG solutions on different meshes at common points."""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from check_bundle import load_case_definition
from comparison.convergence import (
    NEWTON_CONVERGENCE_FAILURE,
    read_newton_convergence,
)
from comparison.metrics import calculate_error_norms, format_metric, maximum_metric
from compare_run import select_candidate
from support.documents import load_json, write_json_atomic
from support.errors import ComparisonError, HarnessError
from support.time import utc_now


@dataclass(frozen=True)
class SampledFields:
    equation_names: list[str]
    inside: np.ndarray
    solution: np.ndarray
    gradient: np.ndarray


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--fekete-nodes", required=True, type=Path)
    parser.add_argument("--samples-per-element", choices=(1, 4), default=4, type=int)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)

    try:
        report = compare_adaptive_files(
            args.reference,
            args.candidate,
            args.fekete_nodes,
            args.samples_per_element,
        )
        if args.report:
            write_json_atomic(args.report, report, "adaptive report")
    except HarnessError as exc:
        print(f"adaptive comparison failed: {exc}", file=sys.stderr)
        return 1

    print(f"adaptive comparison {report['status']}")
    print(
        f"common points: {report['sampling']['common_points']}/"
        f"{report['sampling']['point_count']} "
        f"({report['sampling']['common_coverage']:.3%})"
    )
    for dataset_name, dataset in report["datasets"].items():
        metrics = list(dataset["equations"].values())
        print(
            f"{dataset_name}: max relL2="
            f"{format_metric(maximum_metric(metrics, 'relative_l2'))} "
            f"max nLinf={format_metric(maximum_metric(metrics, 'normalized_linf'))}"
        )
    if args.report:
        print(f"report: {args.report.expanduser().resolve()}")
    return 0 if report["status"] != "failed" else 1


def compare_adaptive_files(
    reference_path: Path,
    candidate_path: Path,
    fekete_path: Path,
    samples_per_element: int = 4,
    tolerances: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Interpolate both files at fixed interior points and compare their fields."""
    reference_path = _existing_file(reference_path, "reference")
    candidate_path = _existing_file(candidate_path, "candidate")
    fekete_path = _existing_file(fekete_path, "Fekete-node")
    points = reference_sample_points(reference_path, samples_per_element)
    reference = _sample_solution(reference_path, points, fekete_path)
    candidate = _sample_solution(candidate_path, points, fekete_path)
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
    run_directory: Path,
    case_dir: Path,
    tolerances_path: Path,
    report_path: Path | None = None,
    candidate_override: Path | None = None,
    reference_override: Path | None = None,
    profile_override: str | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Compare one completed adaptive run with its bundled reference."""
    run_directory = _existing_directory(run_directory, "run")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")

    case = load_case_definition(plan["case_id"], case_dir)
    workflow = case["workflows"].get(plan["workflow_id"])
    if workflow is None or workflow.get("comparison_policy") != "mesh_independent":
        raise ComparisonError("run does not define mesh-independent comparison")
    profile_id, tolerances = _adaptive_tolerance_profile(
        tolerances_path, workflow, profile_override
    )
    candidate = select_candidate(run_directory, metadata, candidate_override)
    reference = _run_input(
        run_directory, reference_override, "inputs/reference.h5", "reference"
    )
    fekete = _fekete_nodes(metadata)
    report = compare_adaptive_files(
        reference,
        candidate,
        fekete,
        tolerances["samples_per_element"],
        tolerances,
    )

    convergence = read_newton_convergence(
        run_directory / "stdout.log", tolerances["newton_error_max"]
    )
    if not convergence.passed:
        report["failures"].append(NEWTON_CONVERGENCE_FAILURE)
    report.update(
        {
            "created_utc": utc_now(),
            "status": "passed" if not report["failures"] else "failed",
            "run_directory": str(run_directory),
            "case_id": plan["case_id"],
            "workflow_id": plan["workflow_id"],
            "layout_id": plan["layout_id"],
            "tolerance_profile": {"id": profile_id, **tolerances},
            "convergence": convergence.as_report(),
        }
    )
    report.pop("tolerances", None)
    report_path = (report_path or run_directory / "comparison.json").expanduser().resolve()
    if report_path in {candidate, reference, fekete}:
        raise ComparisonError("adaptive report cannot replace a comparison input")
    write_json_atomic(report_path, report, "adaptive report")
    return report_path, report


def reference_sample_points(path: Path, samples_per_element: int = 4) -> np.ndarray:
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


def compare_sampled_fields(
    reference: SampledFields,
    candidate: SampledFields,
    samples_per_element: int,
    tolerances: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare already sampled conservative state and gradient arrays."""
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
            metrics = _numeric_metrics(first[common, index], second[common, index], limits)
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


def _sample_solution(path: Path, points: np.ndarray, fekete_path: Path) -> SampledFields:
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


def _numeric_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    limits: dict[str, float] | None,
) -> dict[str, Any]:
    norms = calculate_error_norms(reference, candidate)
    finite = norms.finite
    relative_l2 = norms.relative_l2
    normalized_linf = norms.normalized_linf
    passed = None
    if limits is not None:
        passed = bool(
            norms.available
            and relative_l2 <= limits["relative_l2_max"]
            and normalized_linf <= limits["normalized_linf_max"]
        )
    return {
        "finite": finite,
        "relative_l2": relative_l2,
        "normalized_linf": normalized_linf,
        "passed": passed,
    }


def _validate_sample_shapes(
    reference: SampledFields, candidate: SampledFields
) -> None:
    expected = reference.solution.shape
    if expected != candidate.solution.shape or len(expected) != 2:
        raise ComparisonError("sampled solution shapes differ")
    gradient_shape = (*expected, 2)
    if reference.gradient.shape != gradient_shape or candidate.gradient.shape != gradient_shape:
        raise ComparisonError("sampled gradient shapes differ")
    if reference.inside.shape != (expected[0],) or candidate.inside.shape != (
        expected[0],
    ):
        raise ComparisonError("sampling masks have incompatible shapes")
    if len(reference.equation_names) != expected[1] or len(
        candidate.equation_names
    ) != expected[1]:
        raise ComparisonError("equation names do not match sampled fields")


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
        raise ComparisonError("reference mesh/X must be a two-dimensional coordinate array")
    return array.T if array.shape[0] == 2 else array


def _triangle_rows(array: np.ndarray) -> np.ndarray:
    if array.ndim != 2 or 3 not in array.shape:
        raise ComparisonError("reference mesh/Tlin must contain three-node triangles")
    return (array.T if array.shape[0] == 3 else array).astype(np.int64)


def _existing_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise ComparisonError(f"{label} file does not exist: {path}")
    return path


def _existing_directory(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_dir():
        raise ComparisonError(f"{label} directory does not exist: {path}")
    return path


def _fekete_nodes(metadata: dict[str, Any]) -> Path:
    runtime = metadata.get("runtime_files", {}).get(
        "positionFeketeNodesTri2D.h5", {}
    )
    path = runtime.get("path")
    if not isinstance(path, str):
        raise ComparisonError("run metadata has no Fekete-node path")
    return _existing_file(Path(path), "Fekete-node")


def _adaptive_tolerance_profile(
    path: Path,
    workflow: dict[str, Any],
    profile_override: str | None = None,
) -> tuple[str, dict[str, Any]]:
    document = load_json(path, "tolerance definitions")
    profile_id = profile_override or workflow.get("tolerance_profile")
    profile = document.get("profiles", {}).get(profile_id)
    required = {
        "newton_error_max",
        "samples_per_element",
        "minimum_point_coverage",
        "solution",
        "gradient",
    }
    if not isinstance(profile, dict) or not required <= profile.keys():
        raise ComparisonError(f"invalid adaptive tolerance profile: {profile_id}")
    for dataset in ("solution", "gradient"):
        limits = profile[dataset]
        if not isinstance(limits, dict) or not {
            "relative_l2_max",
            "normalized_linf_max",
        } <= limits.keys():
            raise ComparisonError(
                f"adaptive tolerance profile has invalid {dataset} limits"
            )
    return profile_id, profile


def _run_input(
    run_directory: Path,
    override: Path | None,
    default: str,
    label: str,
) -> Path:
    path = override if override is not None else Path(default)
    path = path.expanduser()
    if not path.is_absolute():
        path = run_directory / path
    return _existing_file(path, label)


if __name__ == "__main__":
    raise SystemExit(main())
