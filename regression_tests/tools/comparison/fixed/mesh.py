"""Fixed-mesh connectivity and coordinate comparison."""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np

from comparison.fixed.alignment import (
    MeshAlignment,
    compare_numbering_invariant,
)
from comparison.fixed.data import optional_array, required_array


def compare_mesh(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
    failures: list[str],
) -> tuple[dict[str, Any], MeshAlignment | None]:
    """Compare mesh connectivity using the configured numbering policy."""
    if tolerances["mesh_connectivity"] == "numbering_invariant":
        report, alignment = compare_numbering_invariant(
            reference,
            candidate,
            tolerances["mesh_coordinate_atol"],
        )
    else:
        report = _compare_exact_mesh(reference, candidate, tolerances)
        alignment = None
    _record_mesh_failures(report, failures)
    return report, alignment


def _compare_exact_mesh(
    reference: h5py.File,
    candidate: h5py.File,
    tolerances: dict[str, Any],
) -> dict[str, Any]:
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
