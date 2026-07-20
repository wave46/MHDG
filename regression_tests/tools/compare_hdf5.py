"""Library-independent fixed-mesh comparison for MHDG HDF5 solutions."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from comparison.metrics import calculate_error_norms
from support.errors import ComparisonError


def compare_hdf5_files(
    reference_path: Path,
    candidate_path: Path,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    """Compare grouped or legacy-flat solutions on the same fixed mesh."""
    failures: list[str] = []
    report: dict[str, Any] = {"failures": failures}

    try:
        with h5py.File(reference_path, "r") as reference, h5py.File(
            candidate_path, "r"
        ) as candidate:
            report["formats"] = {
                "reference": _format_name(reference),
                "candidate": _format_name(candidate),
            }
            report["mesh"] = _compare_mesh(
                reference, candidate, tolerances, failures
            )
            report["solution"] = _compare_solution(
                reference, candidate, tolerances, failures
            )
            report["transport_1d"] = _compare_transport(
                reference, candidate, tolerances, failures
            )
    except (OSError, KeyError) as exc:
        raise ComparisonError(f"cannot read HDF5 comparison data: {exc}") from exc
    except ComparisonError as exc:
        failures.append(str(exc))

    report["status"] = "passed" if not failures else "failed"
    return report


def _compare_mesh(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    report: dict[str, Any] = {"connectivity": {}}
    if tolerances["mesh_connectivity"] != "exact":
        raise ComparisonError("only exact mesh connectivity is supported")

    for name in ("T", "Tlin", "Tb"):
        first = _optional_array(reference, "mesh", name)
        second = _optional_array(candidate, "mesh", name)
        if first is None and second is None and name == "Tb":
            continue
        if first is None or second is None:
            passed = False
            details = {"passed": False, "reason": "dataset missing"}
        else:
            passed = first.shape == second.shape and np.array_equal(first, second)
            details = {
                "passed": passed,
                "reference_shape": list(first.shape),
                "candidate_shape": list(second.shape),
            }
        report["connectivity"][name] = details
        if not passed:
            failures.append(f"mesh/{name} connectivity differs")

    first = _array(reference, "mesh", "X")
    second = _array(candidate, "mesh", "X")
    finite = _finite(first) and _finite(second)
    same_shape = first.shape == second.shape
    maximum_error = (
        float(np.max(np.abs(second - first))) if finite and same_shape else None
    )
    passed = (
        finite
        and same_shape
        and maximum_error is not None
        and maximum_error <= tolerances["mesh_coordinate_atol"]
    )
    report["coordinates"] = {
        "passed": passed,
        "finite": finite,
        "reference_shape": list(first.shape),
        "candidate_shape": list(second.shape),
        "maximum_absolute_error": maximum_error,
        "absolute_tolerance": tolerances["mesh_coordinate_atol"],
    }
    if not passed:
        failures.append("mesh/X coordinates exceed tolerance")
    return report


def _compare_solution(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    reference_neq = _equation_count(reference)
    candidate_neq = _equation_count(candidate)
    if reference_neq != candidate_neq:
        raise ComparisonError(
            f"equation counts differ: {reference_neq} != {candidate_neq}"
        )

    reference_names = _valid_equation_names(reference, reference_neq)
    candidate_names = _valid_equation_names(candidate, candidate_neq)
    names = reference_names or candidate_names or [
        f"equation_{index + 1}" for index in range(reference_neq)
    ]
    names_match = not (
        reference_names and candidate_names and reference_names != candidate_names
    )
    if not names_match:
        failures.append("conservative variable names differ")

    datasets = {}
    for dataset_name in ("u", "q", "u_tilde"):
        first = _array(reference, "solution", dataset_name).reshape(-1)
        second = _array(candidate, "solution", dataset_name).reshape(-1)
        dataset_report = {
            "reference_size": int(first.size),
            "candidate_size": int(second.size),
            "equations": {},
        }
        datasets[dataset_name] = dataset_report

        if first.size != second.size or first.size % reference_neq:
            dataset_report["passed"] = False
            failures.append(f"solution/{dataset_name} has incompatible size")
            continue

        first = first.reshape(-1, reference_neq)
        second = second.reshape(-1, reference_neq)
        for index, name in enumerate(names):
            metrics = _numeric_metrics(first[:, index], second[:, index], tolerances)
            dataset_report["equations"][name] = metrics
            if not metrics["passed"]:
                failures.append(f"solution/{dataset_name}/{name} exceeds tolerance")
        dataset_report["passed"] = all(
            result["passed"] for result in dataset_report["equations"].values()
        )

    return {
        "passed": all(dataset["passed"] for dataset in datasets.values())
        and names_match,
        "equation_count": reference_neq,
        "equation_names": names,
        "datasets": datasets,
    }


def _compare_transport(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    reference_group = _transport_group(reference)
    candidate_group = _transport_group(candidate)
    if reference_group is None:
        return {"present": False, "passed": True, "datasets": {}}
    if candidate_group is None:
        failures.append("candidate is missing transport_1d data")
        return {"present": True, "passed": False, "datasets": {}}

    datasets = {}
    for section in ("coefficients", "profiles"):
        if section not in reference_group:
            continue
        candidate_section = candidate_group.get(section)
        for name, dataset in reference_group[section].items():
            if not isinstance(dataset, h5py.Dataset):
                continue
            label = f"{section}/{name}"
            candidate_dataset = (
                candidate_section.get(name)
                if isinstance(candidate_section, h5py.Group)
                else None
            )
            if (
                not isinstance(candidate_section, h5py.Group)
                or not isinstance(candidate_dataset, h5py.Dataset)
            ):
                datasets[label] = {"passed": False, "reason": "dataset missing"}
                failures.append(f"transport_1d/{label} is missing")
                continue
            metrics = _numeric_metrics(
                np.asarray(dataset),
                np.asarray(candidate_dataset),
                tolerances,
            )
            datasets[label] = metrics
            if not metrics["passed"]:
                failures.append(f"transport_1d/{label} exceeds tolerance")

    return {
        "present": True,
        "passed": all(result["passed"] for result in datasets.values()),
        "datasets": datasets,
    }


def _numeric_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    norms = calculate_error_norms(reference, candidate)
    if not norms.compatible or not norms.finite:
        return {
            "passed": False,
            "finite": norms.finite,
            "reference_shape": list(reference.shape),
            "candidate_shape": list(candidate.shape),
            "relative_l2": None,
            "normalized_linf": None,
        }

    relative_l2 = norms.relative_l2
    normalized_linf = norms.normalized_linf
    assert relative_l2 is not None and normalized_linf is not None
    passed = (
        (norms.finite or not tolerances["require_finite"])
        and relative_l2 <= tolerances["relative_l2_max"]
        and normalized_linf <= tolerances["normalized_linf_max"]
    )
    return {
        "passed": passed,
        "finite": norms.finite,
        "reference_shape": list(reference.shape),
        "candidate_shape": list(candidate.shape),
        "relative_l2": relative_l2,
        "relative_l2_max": tolerances["relative_l2_max"],
        "normalized_linf": normalized_linf,
        "normalized_linf_max": tolerances["normalized_linf_max"],
    }


def _equation_count(handle: h5py.File) -> int:
    for path in ("simulation_parameters/Neq", "Neq"):
        if path in handle:
            return int(np.asarray(handle[path]).reshape(-1)[0])

    names = _optional_names(handle)
    if names:
        return len(names)

    values = _array(handle, "solution", "u").size
    elements = int(_array(handle, "mesh", "Nelems").reshape(-1)[0])
    nodes = int(_array(handle, "mesh", "Nnodesperelem").reshape(-1)[0])
    if elements <= 0 or nodes <= 0 or values % (elements * nodes):
        raise ComparisonError("cannot infer the number of equations")
    return values // (elements * nodes)


def _valid_equation_names(handle: h5py.File, count: int) -> list[str]:
    names = _optional_names(handle)
    return names if len(names) == count else []


def _optional_names(handle: h5py.File) -> list[str]:
    for path in (
        "simulation_parameters/physics/conservative_variable_names",
        "conservative_variable_names",
    ):
        if path not in handle:
            continue
        return [
            value.decode(errors="replace").strip(" \x00")
            if isinstance(value, bytes)
            else str(value).strip()
            for value in np.asarray(handle[path]).reshape(-1)
        ]
    return []


def _array(handle: h5py.File, group: str, name: str) -> np.ndarray:
    array = _optional_array(handle, group, name)
    if array is None:
        raise ComparisonError(f"required dataset is missing: {group}/{name}")
    return array


def _optional_array(
    handle: h5py.File, group: str, name: str
) -> np.ndarray | None:
    grouped_path = f"{group}/{name}"
    for path in (grouped_path, name):
        if path in handle and isinstance(handle[path], h5py.Dataset):
            return np.asarray(handle[path])
    return None


def _transport_group(handle: h5py.File) -> h5py.Group | None:
    for path in ("transport_1d", "solution/transport_1d"):
        if path in handle and isinstance(handle[path], h5py.Group):
            return handle[path]
    return None


def _format_name(handle: h5py.File) -> str:
    solution = "grouped_solution" if "solution" in handle else "flat_solution"
    mesh = "grouped_mesh" if "mesh" in handle else "flat_mesh"
    return f"{solution}+{mesh}"


def _finite(array: np.ndarray) -> bool:
    return bool(np.isfinite(np.asarray(array, dtype=float)).all())
