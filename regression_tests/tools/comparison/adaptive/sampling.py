"""Sample MHDG solutions at deterministic points on a reference mesh."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from support.errors import ComparisonError


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
