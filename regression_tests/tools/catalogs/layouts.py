"""Derive execution settings from declarative layout identifiers."""

from __future__ import annotations

import re
from itertools import combinations
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


def layout_pairs(
    layouts: dict[str, dict[str, Any]], relations: list[str],
) -> list[dict[str, str]]:
    """Generate requested comparisons, sharing each unordered pair once."""
    pairs = []
    seen = set()
    for relation in relations:
        if relation == "all_pairs":
            matches = list(combinations(layouts, 2))
        elif relation in {"openmp", "mpi", "hybrid"}:
            matches = []
            for baseline, first in layouts.items():
                if first["omp_threads"] != 1:
                    continue
                for candidate, second in layouts.items():
                    if relation == "openmp":
                        match = (
                            first["execution"] == second["execution"] == "serial"
                            and second["omp_threads"] > 1
                        )
                    elif relation == "mpi":
                        match = (
                            first["execution"] == "serial"
                            and second["execution"] == "mpi"
                            and second["omp_threads"] == 1
                        )
                    else:
                        match = (
                            first["execution"] == second["execution"] == "mpi"
                            and first["mpi_ranks"] == second["mpi_ranks"]
                            and second["omp_threads"] > 1
                        )
                    if match:
                        matches.append((baseline, candidate))
        else:
            raise BundleError(f"unknown layout relation: {relation}")
        if not matches:
            raise BundleError(f"layout relation {relation} has no matching layouts")
        for baseline, candidate in matches:
            identity = frozenset((baseline, candidate))
            if identity not in seen:
                pairs.append({"baseline": baseline, "candidate": candidate})
                seen.add(identity)
    return pairs


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
