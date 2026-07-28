"""Shared HDF5 dataset access and numeric comparison reports."""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np

from comparison.shared.metrics import calculate_error_norms
from support.errors import ComparisonError


def required_array(handle: h5py.File, group: str, name: str) -> np.ndarray:
    """Read a required dataset from grouped or legacy-flat storage."""
    array = optional_array(handle, group, name)
    if array is None:
        raise ComparisonError(f"required dataset is missing: {group}/{name}")
    return array


def required_scalar(handle: h5py.File, group: str, name: str) -> int | float:
    """Read the scalar value stored in a required MHDG dataset."""
    return required_array(handle, group, name).reshape(-1)[0].item()


def mesh_rows(
    handle: h5py.File,
    name: str,
    dtype: type[float] | type[np.int64],
) -> np.ndarray:
    """Read an MHDG mesh matrix with one entity per returned row."""
    return np.asarray(required_array(handle, "mesh", name), dtype=dtype).T


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
