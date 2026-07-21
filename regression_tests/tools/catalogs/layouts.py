"""Derive execution settings from declarative layout identifiers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from bundle.schemas import load_validated_json
from support.errors import BundleError


SERIAL_LAYOUT = re.compile(r"serial_omp([1-9][0-9]*)")
MPI_LAYOUT = re.compile(r"mpi([1-9][0-9]*)_omp([1-9][0-9]*)")


def load_layouts(path: Path) -> dict[str, dict[str, Any]]:
    """Load every declared layout and derive its execution parameters."""
    schema = path.parent / "schemas/layouts.schema.json"
    document = load_validated_json(path, schema, "layout definitions")
    return {layout_id: _layout(layout_id) for layout_id in document["layouts"]}


def load_layout(layout_id: str, path: Path) -> dict[str, Any]:
    """Return one layout or report the available identifiers."""
    layouts = load_layouts(path)
    try:
        return layouts[layout_id]
    except KeyError as exc:
        available = ", ".join(layouts)
        raise BundleError(
            f"unknown layout {layout_id}; available: {available}"
        ) from exc


def _layout(layout_id: str) -> dict[str, Any]:
    serial = SERIAL_LAYOUT.fullmatch(layout_id)
    if serial:
        threads = int(serial.group(1))
        return {
            "description": f"Serial solver with {threads} OpenMP thread(s)",
            "execution": "serial",
            "mpi_ranks": 1,
            "omp_threads": threads,
        }

    mpi = MPI_LAYOUT.fullmatch(layout_id)
    if mpi:
        ranks, threads = (int(value) for value in mpi.groups())
        return {
            "description": (
                f"Parallel solver with {ranks} MPI rank(s) and "
                f"{threads} OpenMP thread(s)"
            ),
            "execution": "mpi",
            "mpi_ranks": ranks,
            "omp_threads": threads,
        }
    raise BundleError(f"invalid layout identifier: {layout_id}")
