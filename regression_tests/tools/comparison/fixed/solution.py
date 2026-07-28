"""Fixed-mesh conservative solution comparison."""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np

from comparison.fixed.alignment import MeshAlignment
from comparison.fixed.data import (
    numeric_comparison,
    required_array,
    required_scalar,
)
from support.errors import ComparisonError


def compare_solution(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
    alignment: MeshAlignment | None = None,
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
            alignment,
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
    alignment: MeshAlignment | None,
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
        "mesh_alignment_applied": alignment is not None,
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
    if alignment is not None:
        candidate_values = _align_dataset(
            candidate_values,
            dataset_name,
            alignment,
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


def _align_dataset(
    candidate_values: np.ndarray,
    name: str,
    alignment: MeshAlignment,
) -> np.ndarray:
    aligned_candidate = np.empty_like(candidate_values)
    if name in ("u", "q"):
        element_map = alignment.candidate_element_to_reference
        aligned_candidate[element_map] = candidate_values
        return aligned_candidate

    face_map = alignment.candidate_face_to_reference
    aligned_candidate[face_map] = candidate_values
    reversed_faces = np.flatnonzero(alignment.reverse_candidate_faces)
    aligned_candidate[face_map[reversed_faces]] = candidate_values[
        reversed_faces, ::-1
    ]
    return aligned_candidate


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
