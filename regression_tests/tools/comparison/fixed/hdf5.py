"""Coordinate fixed-mesh comparison of two MHDG HDF5 solutions."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py

from comparison.fixed.data import storage_format
from comparison.fixed.mesh import compare_mesh
from comparison.fixed.solution import compare_solution
from comparison.fixed.transport import compare_transport
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
    except (OSError, KeyError) as exc:
        raise ComparisonError(f"cannot read HDF5 comparison data: {exc}") from exc
    except ComparisonError as exc:
        failures.append(str(exc))

    report["status"] = "passed" if not failures else "failed"
    return report
