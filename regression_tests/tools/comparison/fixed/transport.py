"""Optional one-dimensional transport-data comparison."""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np

from comparison.fixed.data import numeric_comparison, optional_group


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
