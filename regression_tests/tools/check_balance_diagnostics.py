#!/usr/bin/env python3
"""Validate detailed particle diagnostics in a completed staged suite."""

from __future__ import annotations

import argparse
import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import h5py

from support.documents import load_json, write_json_atomic
from support.errors import HarnessError
from support.paths import require_directory, require_file


IDENTITY_TOLERANCE = 1.0e-12
TERMINAL_TOLERANCE = 5.1e-4
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
VALUE_LINE = re.compile(rf"^\s*(?P<label>.*?)\s+(?P<value>{NUMBER})\s*$")
TIME_ITERATION = re.compile(r"Time iteration\s*=\s*(\d+)")
NEWTON_ITERATION = re.compile(r"NR iteration:\s*(\d+)")
BLOCK_START = "Particle diagnostics (detailed)"

SPECIES = ("plasma", "neutral", "total")
RATE_GROUPS = (
    "temporal",
    "volume",
    "boundary_physical_inward",
    "physical_imbalance",
    "hdg_tau_inward",
    "conservation",
)
PARTICLE_GROUPS = ("content", *RATE_GROUPS)
COMPONENT_GROUPS = {
    "particles/volume/plasma": "particles/components/plasma/volume",
    "particles/boundary_physical_inward/plasma": (
        "particles/components/plasma/boundary_inward"
    ),
    "particles/volume/neutral": "particles/components/neutral/volume",
    "particles/boundary_physical_inward/neutral": (
        "particles/components/neutral/boundary_inward"
    ),
}
REQUIRED_COMPONENTS = {
    "particles/components/plasma/volume": (
        "ionization",
        "recombination",
        "other_source",
    ),
    "particles/components/plasma/boundary_inward": (
        "parallel",
        "diffusion",
        "pinch",
    ),
    "particles/components/neutral/volume": (
        "ionization",
        "recombination",
        "other_source",
        "puff_source",
        "pump_source",
    ),
    "particles/components/neutral/boundary_inward": ("limited_diffusion",),
    "particles/exchange": ("charge_exchange",),
    "wall_closure/plasma_particles": (
        "diffusion_inward",
        "stabilization_inward",
        "residual",
    ),
    "wall_closure/neutral": (
        "source_inward",
        "physical_flux_inward",
        "stabilization_inward",
        "residual",
        "recycled_plasma_inward",
        "puff_source",
        "pump_sink",
    ),
    "wall_closure/neutral/recycled_plasma_components": (
        "parallel_source",
        "diffusion_source",
        "pinch_source",
    ),
    "wall_closure/neutral/physical_flux_components": (
        "limited_diffusion_inward",
    ),
}


def check_suite(summary_path: Path) -> dict[str, Any]:
    """Check every completed stage and retain its terminal history."""
    summary_path = require_file(summary_path, "balance diagnostics suite summary")
    summary = load_json(summary_path, "balance diagnostics suite summary")
    reports = []
    failures = []
    for result in summary.get("results", []):
        workflow = result.get("workflow_id")
        layout = result.get("layout_id")
        if result.get("run_status") == "completed":
            run_directory = require_directory(
                Path(result["run_directory"]),
                f"balance diagnostics run {workflow}/{layout}",
            )
            report = check_run(run_directory)
        else:
            failure = f"solver run did not complete: {workflow}/{layout}"
            report = {"status": "failed", "failures": [failure], "stages": []}
        reports.append({"workflow_id": workflow, "layout_id": layout, **report})
        failures.extend(
            f"{workflow}/{layout}: {failure}" for failure in report["failures"]
        )
    if not reports:
        failures.append("suite summary contains no run results")
    return {
        "suite_summary": str(summary_path),
        "identity_tolerance": IDENTITY_TOLERANCE,
        "terminal_tolerance": TERMINAL_TOLERANCE,
        "runs": reports,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def check_run(run_directory: Path) -> dict[str, Any]:
    """Check every completed stage output in one staged workflow."""
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    stages = metadata.get("stages")
    reports = []
    failures = []
    if not isinstance(stages, list) or not stages:
        failures.append("run metadata contains no staged outputs")
    for stage in stages or []:
        stage_id = stage.get("stage_id")
        if stage.get("status") != "completed" or not stage.get("selected_hdf5"):
            failure = f"stage is incomplete or has no selected HDF5: {stage_id}"
            report = {"stage_id": stage_id, "status": "failed", "failures": [failure]}
        else:
            stage_directory = require_directory(
                Path(stage["run_directory"]), f"diagnostic stage {stage_id}"
            )
            solution = Path(stage["selected_hdf5"])
            if not solution.is_absolute():
                solution = stage_directory / solution
            report = check_stage(
                stage_id,
                require_file(solution, f"diagnostic stage solution {stage_id}"),
                require_file(stage_directory / "stdout.log", f"stage log {stage_id}"),
            )
        reports.append(report)
        failures.extend(f"{stage_id}: {failure}" for failure in report["failures"])
    return {
        "run_directory": str(run_directory),
        "stages": reports,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def check_stage(stage_id: str, solution: Path, terminal: Path) -> dict[str, Any]:
    """Validate one detailed HDF5 file and its matching terminal block."""
    failures: list[str] = []
    values = _read_hdf5(solution, failures)
    if values:
        _check_particle_identities(values, failures)
        _check_wall_identities(values, failures)
    history = parse_terminal_history(terminal)
    if not history:
        failures.append("terminal log contains no detailed diagnostic blocks")
    elif values:
        _check_terminal_values(history[-1], values, failures)
    return {
        "stage_id": stage_id,
        "solution": str(solution),
        "terminal": str(terminal),
        "terminal_history": history,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def _read_hdf5(path: Path, failures: list[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    units: dict[str, str] = {}
    try:
        with h5py.File(path, "r") as handle:
            mode = _text(
                handle.get("simulation_parameters/switches/balance_diagnostics_mode")
            )
            if mode != "detailed":
                failures.append(
                    f"balance_diagnostics_mode is {mode!r}, expected 'detailed'"
                )
            root = handle.get("diagnostics")
            if not isinstance(root, h5py.Group):
                failures.append("required group is missing: /diagnostics")
                return {}

            def collect(name: str, item: h5py.Group | h5py.Dataset) -> None:
                if not isinstance(item, h5py.Dataset) or item.size != 1:
                    return
                value = _dataset_value(item)
                if name.endswith("/units") and isinstance(value, str):
                    units[name.removesuffix("/units")] = value
                elif isinstance(value, float):
                    values[name] = value

            root.visititems(collect)
    except OSError as exc:
        failures.append(f"cannot read HDF5 solution: {exc}")
        return {}

    expected_units = {
        **{
            f"particles/{group}": "particles" if group == "content" else "particles/s"
            for group in PARTICLE_GROUPS
        },
        **{group: "particles/s" for group in REQUIRED_COMPONENTS},
    }
    for group, expected in expected_units.items():
        if units.get(group) != expected:
            failures.append(
                f"/diagnostics/{group}/units is {units.get(group)!r}, "
                f"expected {expected!r}"
            )

    required = {
        *(f"particles/{group}/{species}" for group in PARTICLE_GROUPS for species in SPECIES),
        *(f"{group}/{name}" for group, names in REQUIRED_COMPONENTS.items() for name in names),
    }
    missing = sorted(required - values.keys())
    failures.extend(f"required scalar dataset is missing: /diagnostics/{key}" for key in missing)
    return values if not missing else {}


def _check_particle_identities(values: dict[str, float], failures: list[str]) -> None:
    for group in PARTICLE_GROUPS:
        species_terms = tuple(values[f"particles/{group}/{name}"] for name in SPECIES[:2])
        activity = _total_activity(group, values, species_terms)
        _identity(
            f"particles/{group}/total",
            values[f"particles/{group}/total"],
            sum(species_terms),
            activity,
            failures,
        )
    for species in SPECIES:
        temporal = values[f"particles/temporal/{species}"]
        volume = values[f"particles/volume/{species}"]
        boundary = values[f"particles/boundary_physical_inward/{species}"]
        physical = values[f"particles/physical_imbalance/{species}"]
        tau = values[f"particles/hdg_tau_inward/{species}"]
        _identity(
            f"particles/physical_imbalance/{species}",
            physical,
            temporal - volume - boundary,
            (temporal, volume, boundary),
            failures,
        )
        _identity(
            f"particles/conservation/{species}",
            values[f"particles/conservation/{species}"],
            physical - tau,
            (physical, tau),
            failures,
        )
    for aggregate, prefix in COMPONENT_GROUPS.items():
        terms = _terms(values, prefix)
        _identity(aggregate, values[aggregate], sum(terms), terms, failures)
    for reaction in ("ionization", "recombination"):
        terms = tuple(
            values[f"particles/components/{species}/volume/{reaction}"]
            for species in SPECIES[:2]
        )
        _identity(
            f"particles/components/{reaction}/plasma_plus_neutral",
            sum(terms),
            0.0,
            terms,
            failures,
        )


def _total_activity(
    group: str, values: dict[str, float], species_terms: Sequence[float]
) -> tuple[float, ...]:
    """Scale cancellation-prone totals by the terms that formed them."""
    dependencies = {
        "physical_imbalance": ("temporal", "volume", "boundary_physical_inward"),
        "conservation": ("physical_imbalance", "hdg_tau_inward"),
    }
    return tuple(
        values[f"particles/{term}/{species}"]
        for term in dependencies.get(group, ())
        for species in SPECIES[:2]
    ) or tuple(species_terms)


def _check_wall_identities(values: dict[str, float], failures: list[str]) -> None:
    plasma = "wall_closure/plasma_particles"
    plasma_terms = (
        values[f"{plasma}/diffusion_inward"],
        values[f"{plasma}/stabilization_inward"],
    )
    _identity(
        f"{plasma}/residual",
        values[f"{plasma}/residual"],
        sum(plasma_terms),
        plasma_terms,
        failures,
    )

    neutral = "wall_closure/neutral"
    source_terms = (
        values[f"{neutral}/recycled_plasma_inward"],
        values[f"{neutral}/puff_source"],
        -values[f"{neutral}/pump_sink"],
    )
    _identity(
        f"{neutral}/source_inward",
        values[f"{neutral}/source_inward"],
        sum(source_terms),
        source_terms,
        failures,
    )
    for aggregate, prefix in (
        ("recycled_plasma_inward", f"{neutral}/recycled_plasma_components"),
        ("physical_flux_inward", f"{neutral}/physical_flux_components"),
    ):
        terms = _terms(values, prefix)
        _identity(
            f"{neutral}/{aggregate}",
            values[f"{neutral}/{aggregate}"],
            sum(terms),
            terms,
            failures,
        )
    residual_terms = (
        values[f"{neutral}/source_inward"],
        -values[f"{neutral}/physical_flux_inward"],
        -values[f"{neutral}/stabilization_inward"],
    )
    _identity(
        f"{neutral}/residual",
        values[f"{neutral}/residual"],
        sum(residual_terms),
        residual_terms,
        failures,
    )


def parse_terminal_history(path: Path) -> list[dict[str, Any]]:
    """Return aggregate values from every detailed terminal block."""
    history = []
    current: dict[str, Any] | None = None
    section = ""
    species = ""
    wall_species = ""
    time_iteration = None
    newton_iteration = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := TIME_ITERATION.search(line):
            time_iteration = int(match.group(1))
        if match := NEWTON_ITERATION.search(line):
            newton_iteration = int(match.group(1))
        stripped = line.strip()
        if stripped == BLOCK_START:
            if current is not None:
                history.append(current)
            current = {
                "time_iteration": time_iteration,
                "newton_iteration": newton_iteration,
                "content": {},
                "conservation": {name: {} for name in SPECIES},
                "wall_residual": {},
            }
            section = ""
            species = ""
            wall_species = ""
            continue
        if current is None:
            continue
        if stripped == "Content [particles]":
            section = "content"
            continue
        if stripped.endswith(" conservation [particles/s]"):
            candidate = stripped.split()[0]
            if candidate in SPECIES:
                section, species = "conservation", candidate
            continue
        if stripped.startswith("Wall / HDG boundary-condition residuals"):
            section = "wall"
            continue
        if stripped.startswith("Plasma density BC:"):
            wall_species = "plasma"
            continue
        if stripped.startswith("Neutral density BC:"):
            wall_species = "neutral"
            continue
        match = VALUE_LINE.match(line)
        if match is None:
            continue
        label = match.group("label").strip()
        value = float(match.group("value").replace("D", "E").replace("d", "e"))
        if section == "content":
            labels = {"plasma n": "plasma", "neutral nn": "neutral", "total n+nn": "total"}
            if label in labels:
                current["content"][labels[label]] = value
        elif section == "conservation" and species:
            labels = {
                "temporal": "temporal",
                "volume": "volume",
                "boundary physical, inward": "boundary_physical_inward",
                "physical imbalance": "physical_imbalance",
                "HDG tau inward": "hdg_tau_inward",
                "discrete residual": "conservation",
            }
            if label in labels:
                current["conservation"][species][labels[label]] = value
        elif section == "wall" and label == "residual" and wall_species:
            current["wall_residual"][wall_species] = value
    if current is not None:
        history.append(current)
    return history


def _check_terminal_values(
    block: dict[str, Any], values: dict[str, float], failures: list[str]
) -> None:
    terminal_values = {
        **{
            f"particles/content/{species}": value
            for species, value in block["content"].items()
        },
        **{
            f"particles/{term}/{species}": value
            for species, terms in block["conservation"].items()
            for term, value in terms.items()
        },
        **{
            "wall_closure/"
            f"{'plasma_particles' if species == 'plasma' else 'neutral'}"
            f"/residual": value
            for species, value in block["wall_residual"].items()
        },
    }
    expected = {
        *(f"particles/{group}/{species}" for group in PARTICLE_GROUPS for species in SPECIES),
        "wall_closure/plasma_particles/residual",
        "wall_closure/neutral/residual",
    }
    failures.extend(
        f"terminal block is missing {key}" for key in sorted(expected - terminal_values.keys())
    )
    for key in sorted(expected & terminal_values.keys()):
        exact = values[key]
        printed = terminal_values[key]
        difference = abs(exact - printed) / max(abs(exact), abs(printed), 1.0)
        if not math.isfinite(printed) or difference > TERMINAL_TOLERANCE:
            failures.append(
                f"terminal/HDF5 mismatch for {key}: relative difference "
                f"{difference:.6e} exceeds {TERMINAL_TOLERANCE:.6e}"
            )


def _terms(values: dict[str, float], prefix: str) -> tuple[float, ...]:
    return tuple(value for key, value in values.items() if key.startswith(prefix + "/"))


def _identity(
    label: str,
    actual: float,
    expected: float,
    activity: Sequence[float],
    failures: list[str],
) -> None:
    if not all(math.isfinite(value) for value in (actual, expected, *activity)):
        return
    scale = max(sum(abs(value) for value in activity), abs(actual), 1.0)
    difference = abs(actual - expected) / scale
    if difference > IDENTITY_TOLERANCE:
        failures.append(
            f"identity failure for {label}: relative difference "
            f"{difference:.6e} exceeds {IDENTITY_TOLERANCE:.6e}"
        )


def _dataset_value(dataset: h5py.Dataset) -> float | str | None:
    value = dataset[()]
    if hasattr(value, "reshape"):
        value = value.reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode("utf-8").strip()
    if isinstance(value, str):
        return value.strip()
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result


def _text(dataset: h5py.Dataset | h5py.Group | None) -> str | None:
    if not isinstance(dataset, h5py.Dataset) or dataset.size != 1:
        return None
    value = _dataset_value(dataset)
    return value if isinstance(value, str) else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_summary", type=Path)
    args = parser.parse_args()
    try:
        report = check_suite(args.suite_summary)
        path = args.suite_summary.parent / "balance_diagnostics_check.json"
        write_json_atomic(path, report, "balance diagnostics report")
    except HarnessError as exc:
        parser.error(str(exc))
    print(f"balance diagnostics check: {report['status']}")
    print(f"report: {path}")
    for failure in report["failures"]:
        print(f"failure: {failure}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
