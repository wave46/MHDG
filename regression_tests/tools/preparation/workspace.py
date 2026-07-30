"""Populate isolated run directories from external bundle artifacts."""

from __future__ import annotations

import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from preparation.parameters import render_parameter_file
from support.errors import BundleError


WARM_INPUT_LINKS = {
    "mesh": "mesh.msh",
    "geometry": "geometry.geo",
    "equilibrium_magnetic_field": "equilibrium.h5",
    "equilibrium_current_density": "current_density.h5",
    "transport_configuration": "transport_model.nml",
}


@contextmanager
def temporary_run_directory(final_directory: Path) -> Iterator[Path]:
    """Build a run beside its destination and publish it atomically."""
    try:
        final_directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{final_directory.name}.",
            dir=final_directory.parent,
        ) as workspace:
            staging = Path(workspace) / "run"
            staging.mkdir()
            yield staging
            staging.rename(final_directory)
    except OSError as exc:
        raise BundleError(
            f"cannot prepare run {final_directory}: {exc}"
        ) from exc


def populate_warm_run(
    staging: Path,
    final_directory: Path,
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    workflow: dict[str, Any],
    parameter_overrides: dict[str, Any],
) -> None:
    """Link warm-run inputs and render its parameter file."""
    inputs = _create_run_directories(staging)
    for role, filename in WARM_INPUT_LINKS.items():
        (inputs / filename).symlink_to(artifacts[role])
    (inputs / "restart.h5").symlink_to(
        artifacts[workflow.get("restart_role", "warm_restart")]
    )
    (inputs / "reference.h5").symlink_to(
        artifacts[workflow.get("reference_role", "warm_reference")]
    )
    _link_runtime_files(staging, runtime_files)
    render_parameter_file(
        artifacts["warm_parameters"],
        staging / "param.txt",
        _parameter_replacements(final_directory),
        parameter_overrides,
    )


def populate_stage_run(
    staging: Path,
    final_directory: Path,
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    workflow: dict[str, Any],
    stage: dict[str, Any],
    parameter_overrides: dict[str, Any],
) -> None:
    """Link one stage's inputs and render its parameter file."""
    inputs = _create_run_directories(staging)
    input_sources = {
        "mesh.msh": artifacts[workflow["mesh_role"]],
        "geometry.geo": artifacts["geometry"],
        "equilibrium.h5": artifacts["equilibrium_magnetic_field"],
        "current_density.h5": artifacts["equilibrium_current_density"],
        "transport_model.nml": artifacts[stage["transport_configuration_role"]],
    }
    for filename, source in input_sources.items():
        (inputs / filename).symlink_to(source)
    _link_runtime_files(staging, runtime_files)
    render_parameter_file(
        artifacts[stage["parameter_role"]],
        staging / "param.txt",
        _parameter_replacements(final_directory),
        parameter_overrides,
    )


def _create_run_directories(staging: Path) -> Path:
    inputs = staging / "inputs"
    inputs.mkdir(parents=True)
    (staging / "outputs").mkdir()
    (staging / "res").mkdir()
    return inputs


def _link_runtime_files(
    staging: Path,
    runtime_files: dict[str, Path],
) -> None:
    for filename, source in runtime_files.items():
        (staging / filename).symlink_to(source)


def _parameter_replacements(final_directory: Path) -> dict[str, Path | str]:
    inputs = final_directory / "inputs"
    return {
        "transport_model_path": inputs / "transport_model.nml",
        "field_path": inputs / "equilibrium.h5",
        "jtor_path": inputs / "current_density.h5",
        "geometry_path": inputs / "geometry.geo",
        "save_folder": f"{final_directory / 'outputs'}/",
    }
