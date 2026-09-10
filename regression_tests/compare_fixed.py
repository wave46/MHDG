"""Direct comparison of mesh, conservative fields, transport and magnetic data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from support.errors import ComparisonError
from .compare_common import calculate_error_norms


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
                "reference": storage_format(reference),
                "candidate": storage_format(candidate),
            }
            report["mesh"] = compare_mesh(
                reference,
                candidate,
                tolerances,
                failures,
            )
            if not report["mesh"]["passed"]:
                report["solution"] = {
                    "passed": False,
                    "reason": "mesh comparison failed",
                }
            else:
                report["solution"] = compare_solution(
                    reference,
                    candidate,
                    tolerances,
                    failures,
                )
            report["transport_1d"] = compare_transport(
                reference,
                candidate,
                tolerances,
                failures,
            )
            report["magnetic"] = compare_magnetic(
                reference,
                candidate,
                tolerances,
                failures,
            )
    except (OSError, KeyError) as exc:
        raise ComparisonError(f"cannot read HDF5 comparison data: {exc}") from exc
    except ComparisonError as exc:
        failures.append(str(exc))

    report["status"] = "passed" if not failures else "failed"
    return report


def required_array(handle: h5py.File, group: str, name: str) -> np.ndarray:
    """Read a required dataset from grouped or legacy-flat storage."""
    array = optional_array(handle, group, name)
    if array is None:
        raise ComparisonError(f"required dataset is missing: {group}/{name}")
    return array


def required_scalar(handle: h5py.File, group: str, name: str) -> int | float:
    """Read the scalar value stored in a required MHDG dataset."""
    return required_array(handle, group, name).reshape(-1)[0].item()


def optional_array(
    handle: h5py.File,
    group: str,
    name: str,
) -> np.ndarray | None:
    """Read an optional dataset from grouped or legacy-flat storage."""
    for path in (f"{group}/{name}", name):
        if path in handle and isinstance(handle[path], h5py.Dataset):
            return np.asarray(handle[path])
    return None


def optional_group(handle: h5py.File, *paths: str) -> h5py.Group | None:
    """Return the first matching HDF5 group."""
    for path in paths:
        if path in handle and isinstance(handle[path], h5py.Group):
            return handle[path]
    return None


def storage_format(handle: h5py.File) -> str:
    """Describe whether solution and mesh datasets are grouped or flat."""
    solution = "grouped_solution" if "solution" in handle else "flat_solution"
    mesh = "grouped_mesh" if "mesh" in handle else "flat_mesh"
    return f"{solution}+{mesh}"


def numeric_comparison(
    reference: np.ndarray,
    candidate: np.ndarray,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    """Build the standard fixed-mesh numeric tolerance report."""
    norms = calculate_error_norms(reference, candidate)
    if not norms.available:
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
        relative_l2 <= tolerances["relative_l2_max"]
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


def compare_mesh(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare mesh coordinates and connectivity in storage order."""
    connectivity = {}
    for name in ("T", "Tlin", "Tb"):
        details = _compare_connectivity(reference, candidate, name)
        if details is None:
            continue
        connectivity[name] = details

    coordinates = _compare_coordinates(reference, candidate, tolerances)
    report = {
        "mode": "exact",
        "connectivity": connectivity,
        "coordinates": coordinates,
    }
    report["passed"] = coordinates["passed"] and all(
        details["passed"] for details in connectivity.values()
    )
    _record_mesh_failures(report, failures)
    return report


def _record_mesh_failures(
    report: dict[str, Any], failures: list[str]
) -> None:
    coordinates = report.get("coordinates", {})
    if not coordinates.get("passed"):
        failures.append(
            f"mesh/X: {coordinates.get('reason', 'coordinates differ')}"
        )
        return
    connectivity = report.get("connectivity", {})
    for name, details in connectivity.items():
        if not details.get("passed"):
            reason = details.get("reason", "connectivity differs")
            failures.append(f"mesh/{name}: {reason}")


def _compare_connectivity(
    reference: h5py.File,
    candidate: h5py.File,
    name: str,
) -> dict[str, Any] | None:
    first = optional_array(reference, "mesh", name)
    second = optional_array(candidate, "mesh", name)
    if first is None and second is None and name == "Tb":
        return None
    if first is None or second is None:
        return {"passed": False, "reason": "dataset missing"}

    return {
        "passed": first.shape == second.shape and np.array_equal(first, second),
        "reference_shape": list(first.shape),
        "candidate_shape": list(second.shape),
    }


def _compare_coordinates(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
    first = required_array(reference, "mesh", "X")
    second = required_array(candidate, "mesh", "X")
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
    return {
        "passed": passed,
        "finite": finite,
        "reference_shape": list(first.shape),
        "candidate_shape": list(second.shape),
        "maximum_absolute_error": maximum_error,
        "absolute_tolerance": tolerances["mesh_coordinate_atol"],
    }


def _finite(array: np.ndarray) -> bool:
    return bool(np.isfinite(np.asarray(array, dtype=float)).all())


def compare_solution(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare each conservative field equation by equation."""
    equation_count = _matching_equation_count(reference, candidate)
    names, names_match = _equation_names(reference, candidate, equation_count)
    if not names_match:
        failures.append("conservative variable names differ")

    datasets = {
        name: _compare_dataset(
            reference,
            candidate,
            name,
            equation_count,
            names,
            tolerances,
            failures,
        )
        for name in ("u", "q", "u_tilde")
    }
    return {
        "passed": all(dataset["passed"] for dataset in datasets.values())
        and names_match,
        "equation_count": equation_count,
        "equation_names": names,
        "datasets": datasets,
    }


def _compare_dataset(
    reference: h5py.File,
    candidate: h5py.File,
    dataset_name: str,
    equation_count: int,
    equation_names: list[str],
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    reference_values = required_array(
        reference, "solution", dataset_name
    ).reshape(-1)
    candidate_values = required_array(
        candidate, "solution", dataset_name
    ).reshape(-1)
    report: dict[str, Any] = {
        "reference_size": int(reference_values.size),
        "candidate_size": int(candidate_values.size),
        "equations": {},
    }
    if (
        reference_values.size != candidate_values.size
        or reference_values.size % equation_count
    ):
        report["passed"] = False
        failures.append(f"solution/{dataset_name} has incompatible size")
        return report

    reference_values = _reshape_dataset(
        reference, dataset_name, reference_values, equation_count
    )
    candidate_values = _reshape_dataset(
        candidate, dataset_name, candidate_values, equation_count
    )

    reference_by_equation = np.moveaxis(reference_values, 2, -1).reshape(
        -1, equation_count
    )
    candidate_by_equation = np.moveaxis(candidate_values, 2, -1).reshape(
        -1, equation_count
    )
    for index, name in enumerate(equation_names):
        metrics = numeric_comparison(
            reference_by_equation[:, index],
            candidate_by_equation[:, index],
            tolerances,
        )
        report["equations"][name] = metrics
        if not metrics["passed"]:
            failures.append(f"solution/{dataset_name}/{name} exceeds tolerance")
    report["passed"] = all(
        result["passed"] for result in report["equations"].values()
    )
    return report


def _reshape_dataset(
    handle: h5py.File,
    name: str,
    values: np.ndarray,
    equation_count: int,
) -> np.ndarray:
    element_count = int(required_scalar(handle, "mesh", "Nelems"))
    nodes_per_element = int(
        required_scalar(handle, "mesh", "Nnodesperelem")
    )
    if name == "u":
        shape = (element_count, nodes_per_element, equation_count)
    elif name == "q":
        dimension = int(required_scalar(handle, "mesh", "Ndim"))
        shape = (element_count, nodes_per_element, equation_count, dimension)
    else:
        face_count = int(required_scalar(handle, "mesh", "Nfaces"))
        nodes_per_face = int(
            required_scalar(handle, "mesh", "Nnodesperface")
        )
        shape = (face_count, nodes_per_face, equation_count)
    return values.reshape(shape)


def _matching_equation_count(reference: h5py.File, candidate: h5py.File) -> int:
    reference_count = _equation_count(reference)
    candidate_count = _equation_count(candidate)
    if reference_count != candidate_count:
        raise ComparisonError(
            f"equation counts differ: {reference_count} != {candidate_count}"
        )
    return reference_count


def _equation_names(
    reference: h5py.File,
    candidate: h5py.File,
    count: int,
) -> tuple[list[str], bool]:
    reference_names = _valid_equation_names(reference, count)
    candidate_names = _valid_equation_names(candidate, count)
    names = reference_names or candidate_names or [
        f"equation_{index + 1}" for index in range(count)
    ]
    names_match = not (
        reference_names and candidate_names and reference_names != candidate_names
    )
    return names, names_match


def _equation_count(handle: h5py.File) -> int:
    for path in ("simulation_parameters/Neq", "Neq"):
        if path in handle:
            return int(np.asarray(handle[path]).reshape(-1)[0])

    names = _optional_names(handle)
    if names:
        return len(names)

    values = required_array(handle, "solution", "u").size
    elements = int(required_scalar(handle, "mesh", "Nelems"))
    nodes = int(required_scalar(handle, "mesh", "Nnodesperelem"))
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


def compare_transport(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare reference transport datasets when they are present."""
    reference_group = optional_group(
        reference,
        "transport_1d",
        "solution/transport_1d",
    )
    candidate_group = optional_group(
        candidate,
        "transport_1d",
        "solution/transport_1d",
    )
    if reference_group is None:
        return {"present": False, "passed": True, "datasets": {}}
    if candidate_group is None:
        failures.append("candidate is missing transport_1d data")
        return {"present": True, "passed": False, "datasets": {}}

    datasets = {}
    for section in ("coefficients", "profiles"):
        _compare_section(
            reference_group,
            candidate_group,
            section,
            tolerances,
            failures,
            datasets,
        )
    return {
        "present": True,
        "passed": all(result["passed"] for result in datasets.values()),
        "datasets": datasets,
    }


def _compare_section(
    reference_group: h5py.Group,
    candidate_group: h5py.Group,
    section: str,
    tolerances: dict[str, Any],
    failures: list[str],
    datasets: dict[str, dict[str, Any]],
) -> None:
    if section not in reference_group:
        return

    candidate_section = candidate_group.get(section)
    for name, reference_dataset in reference_group[section].items():
        if not isinstance(reference_dataset, h5py.Dataset):
            continue
        label = f"{section}/{name}"
        candidate_dataset = (
            candidate_section.get(name)
            if isinstance(candidate_section, h5py.Group)
            else None
        )
        if not isinstance(candidate_dataset, h5py.Dataset):
            datasets[label] = {"passed": False, "reason": "dataset missing"}
            failures.append(f"transport_1d/{label} is missing")
            continue

        metrics = numeric_comparison(
            np.asarray(reference_dataset),
            np.asarray(candidate_dataset),
            tolerances,
        )
        datasets[label] = metrics
        if not metrics["passed"]:
            failures.append(f"transport_1d/{label} exceeds tolerance")


def compare_magnetic(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    """Compare every direct dataset in the reference magnetic group."""
    reference_group = optional_group(reference, "magnetic")
    candidate_group = optional_group(candidate, "magnetic")
    if reference_group is None:
        return {"present": False, "passed": True, "datasets": {}}
    if candidate_group is None:
        failures.append("candidate is missing magnetic data")
        return {"present": True, "passed": False, "datasets": {}}

    datasets: dict[str, dict[str, Any]] = {}
    for name, reference_dataset in reference_group.items():
        if not isinstance(reference_dataset, h5py.Dataset):
            continue
        candidate_dataset = candidate_group.get(name)
        if not isinstance(candidate_dataset, h5py.Dataset):
            datasets[name] = {"passed": False, "reason": "dataset missing"}
            failures.append(f"magnetic/{name} is missing")
            continue

        reference_array = np.asarray(reference_dataset)
        candidate_array = np.asarray(candidate_dataset)
        if reference_array.dtype.kind in "biuOSU":
            metrics = _exact_comparison(reference_array, candidate_array)
        else:
            metrics = numeric_comparison(
                reference_array,
                candidate_array,
                tolerances,
            )
        datasets[name] = metrics
        if not metrics["passed"]:
            failures.append(f"magnetic/{name} differs")

    return {
        "present": True,
        "passed": all(result["passed"] for result in datasets.values()),
        "datasets": datasets,
    }


def _exact_comparison(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict[str, Any]:
    """Compare categorical and integer magnetic metadata exactly."""
    return {
        "passed": np.array_equal(reference, candidate),
        "mode": "exact",
        "reference_shape": list(reference.shape),
        "candidate_shape": list(candidate.shape),
    }
