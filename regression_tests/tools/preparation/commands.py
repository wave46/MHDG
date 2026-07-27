"""Construct deterministic solver commands and runtime environments."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def openmp_environment(threads: int) -> dict[str, str]:
    """Return deterministic OpenMP placement for one solver process."""
    return {
        "OMP_NUM_THREADS": str(threads),
        "OMP_PLACES": "cores",
        "OMP_PROC_BIND": "spread",
    }


def solver_command(
    run_directory: Path,
    executable: Path,
    launcher: Path | None,
    layout: dict[str, Any],
    restart: bool,
) -> list[str]:
    """Build the serial or MPI solver command for one prepared run."""
    arguments = [
        str(executable),
        str(run_directory / "inputs" / "mesh"),
    ]
    if restart:
        arguments.append(str(run_directory / "inputs" / "restart"))
    if launcher is None:
        return arguments
    return [
        str(launcher),
        "--bind-to",
        "core",
        "--map-by",
        f"slot:PE={layout['omp_threads']}",
        "-n",
        str(layout["mpi_ranks"]),
        *arguments,
    ]
