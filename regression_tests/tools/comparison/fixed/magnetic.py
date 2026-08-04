"""Optional magnetic-geometry contract comparison."""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np

from comparison.fixed.data import numeric_comparison, optional_group


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
