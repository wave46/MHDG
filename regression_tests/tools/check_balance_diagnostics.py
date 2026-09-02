#!/usr/bin/env python3
"""Validate detailed equation-oriented balance diagnostics."""

from __future__ import annotations

import argparse
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py

from comparison.shared.outputs import select_candidate
from support.documents import load_json, write_json_atomic
from support.errors import ComparisonError, HarnessError
from support.paths import require_directory, require_file


IDENTITY_TOLERANCE = 1.0e-12
TERMINAL_TOLERANCE = 5.1e-4
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
VALUE_LINE = re.compile(rf"^\s*(?P<label>.*?)\s+(?P<value>{NUMBER})\s*$")
TIME_ITERATION = re.compile(r"Time iteration\s*=\s*(\d+)")
NEWTON_ITERATION = re.compile(r"NR iteration:\s*(\d+)")
BLOCK_START = "Balance diagnostics (detailed)"
EQUATIONS = ("n", "n_n", "total_n")


@dataclass(frozen=True)
class DiagnosticContext:
    neutral_gamma: bool = False
    relocated_sources: bool = False
    configured_puff: float | None = None


def check_suite(summary_path: Path) -> dict[str, Any]:
    summary_path = require_file(summary_path, "balance diagnostics suite summary")
    summary = load_json(summary_path, "balance diagnostics suite summary")
    reports: list[dict[str, Any]] = []
    failures: list[str] = []
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
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    reports: list[dict[str, Any]] = []
    failures: list[str] = []
    stages = metadata.get("stages")
    if isinstance(stages, list) and stages:
        for stage in stages:
            report = _check_staged_output(stage)
            reports.append(report)
            failures.extend(
                f"{report['stage_id']}: {failure}" for failure in report["failures"]
            )
    else:
        try:
            solution = select_candidate(run_directory, metadata)
        except ComparisonError as exc:
            failures.append(f"warm: {exc}")
        else:
            report = check_stage(
                "warm",
                solution,
                require_file(run_directory / "stdout.log", "warm diagnostic log"),
            )
            reports.append(report)
            failures.extend(f"warm: {failure}" for failure in report["failures"])
    return {
        "run_directory": str(run_directory),
        "stages": reports,
        "status": "passed" if not failures else "failed",
        "failures": failures,
    }


def _check_staged_output(stage: dict[str, Any]) -> dict[str, Any]:
    stage_id = stage.get("stage_id")
    if stage.get("status") != "completed" or not stage.get("selected_hdf5"):
        failure = f"stage is incomplete or has no selected HDF5: {stage_id}"
        return {"stage_id": stage_id, "status": "failed", "failures": [failure]}

    stage_directory = require_directory(
        Path(stage["run_directory"]), f"diagnostic stage {stage_id}"
    )
    solution = Path(stage["selected_hdf5"])
    if not solution.is_absolute():
        solution = stage_directory / solution
    return check_stage(
        stage_id,
        require_file(solution, f"diagnostic stage solution {stage_id}"),
        require_file(stage_directory / "stdout.log", f"stage log {stage_id}"),
    )


def check_stage(stage_id: str, solution: Path, terminal: Path) -> dict[str, Any]:
    failures: list[str] = []
    values, context = _read_hdf5(solution, failures)
    if values:
        _check_physical(values, failures)
        _check_bc(values, failures)
        _check_relocated_sources(values, context, failures)
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


def _read_hdf5(
    path: Path, failures: list[str]
) -> tuple[dict[str, float], DiagnosticContext]:
    values: dict[str, float] = {}
    texts: dict[str, str] = {}
    context = DiagnosticContext()
    try:
        with h5py.File(path, "r") as handle:
            mode = _text(
                handle.get("simulation_parameters/switches/balance_diagnostics_mode")
            )
            if mode != "detailed":
                failures.append(
                    f"balance_diagnostics_mode is {mode!r}, expected 'detailed'"
                )
            names = _string_values(
                handle.get("simulation_parameters/physics/conservative_variable_names")
            )
            relocated = _optional_value(
                handle,
                "simulation_parameters/switches/neutral_wall_sources_in_elements",
            )
            configured_puff = _optional_value(
                handle, "simulation_parameters/physics/puff"
            )
            context = DiagnosticContext(
                neutral_gamma="Gamman" in names,
                relocated_sources=relocated == 1.0,
                configured_puff=(configured_puff if isinstance(configured_puff, float) else None),
            )
            root = handle.get("diagnostics")
            if not isinstance(root, h5py.Group):
                failures.append("required group is missing: /diagnostics")
                return {}, context

            def collect(name: str, item: h5py.Group | h5py.Dataset) -> None:
                if not isinstance(item, h5py.Dataset) or item.size != 1:
                    return
                item_value = _dataset_value(item)
                if isinstance(item_value, str):
                    texts[name] = item_value
                elif isinstance(item_value, float):
                    values[name] = item_value

            root.visititems(collect)
    except OSError as exc:
        failures.append(f"cannot read HDF5 solution: {exc}")
        return {}, context

    required = _required_paths(context)
    missing = sorted(required - values.keys())
    failures.extend(
        f"required scalar dataset is missing: /diagnostics/{name}" for name in missing
    )
    nonfinite = sorted(
        name for name in required & values.keys() if not math.isfinite(values[name])
    )
    failures.extend(
        f"required scalar dataset is non-finite: /diagnostics/{name}"
        for name in nonfinite
    )
    expected_units = {
        **{f"physical/{eq}/content_units": "particles" for eq in EQUATIONS},
        **{f"physical/{eq}/rate_units": "particles/s" for eq in EQUATIONS},
        **{f"discrete/{eq}/units": "particles/s" for eq in EQUATIONS},
        "bc/n/units": "particles/s",
        "bc/n_n/units": "particles/s",
    }
    for name, expected in expected_units.items():
        if texts.get(name) != expected:
            failures.append(
                f"/diagnostics/{name} is {texts.get(name)!r}, expected {expected!r}"
            )
    return (values, context) if not missing and not nonfinite else ({}, context)


def _required_paths(context: DiagnosticContext) -> set[str]:
    required = {
        *(f"physical/{eq}/{field}" for eq in EQUATIONS for field in (
            "content", "temporal", "volume", "boundary_physical_inward",
            "physical_imbalance",
        )),
        *(f"discrete/{eq}/{field}" for eq in EQUATIONS for field in (
            "hdg_tau_inward", "residual",
        )),
        "physical/n/volume_components/ionization",
        "physical/n/volume_components/recombination",
        "physical/n/volume_components/prescribed_source",
        "physical/n/boundary_components_inward/parallel",
        "physical/n/boundary_components_inward/diffusion",
        "physical/n/boundary_components_inward/pinch",
        "physical/n/exchange/charge_exchange_rate",
        "physical/n_n/volume_components/ionization",
        "physical/n_n/volume_components/recombination",
        "physical/n_n/volume_components/prescribed_source",
        "physical/n_n/volume_components/puff_source",
        "physical/n_n/volume_components/pump_source",
        "physical/n_n/boundary_components_inward/limited_diffusion",
        "physical/n_n/boundary_components_inward/limited_pressure",
        "bc/n/diffusion_inward",
        "bc/n/hdg_tau_inward",
        "bc/n/residual",
        "bc/n_n/imposed_source_inward",
        "bc/n_n/physical_flux_inward",
        "bc/n_n/hdg_tau_inward",
        "bc/n_n/residual",
        "bc/n_n/source_components/recycling_parallel_inward",
        "bc/n_n/source_components/recycling_diffusion_inward",
        "bc/n_n/source_components/recycling_pinch_inward",
        "bc/n_n/source_components/puff_source",
        "bc/n_n/source_components/pump_sink",
        "bc/n_n/physical_flux_components_inward/limited_diffusion",
        "bc/n_n/physical_flux_components_inward/limited_pressure",
    }
    if context.neutral_gamma:
        required.update(
            {
                "physical/n_n/boundary_components_inward/neutral_gamma_convection",
                "bc/n_n/physical_flux_components_inward/neutral_gamma_convection",
            }
        )
    return required


def _check_physical(values: dict[str, float], failures: list[str]) -> None:
    for equation in EQUATIONS:
        temporal = values[f"physical/{equation}/temporal"]
        volume = values[f"physical/{equation}/volume"]
        boundary = values[f"physical/{equation}/boundary_physical_inward"]
        imbalance = values[f"physical/{equation}/physical_imbalance"]
        tau = values[f"discrete/{equation}/hdg_tau_inward"]
        _identity(
            f"physical/{equation}/physical_imbalance",
            imbalance,
            temporal - volume - boundary,
            (temporal, volume, boundary),
            failures,
        )
        _identity(
            f"discrete/{equation}/residual",
            values[f"discrete/{equation}/residual"],
            imbalance - tau,
            (imbalance, tau),
            failures,
        )
    for field in (
        "content", "temporal", "volume", "boundary_physical_inward",
        "physical_imbalance",
    ):
        terms = tuple(values[f"physical/{eq}/{field}"] for eq in EQUATIONS[:2])
        _identity(
            f"physical/total_n/{field}",
            values[f"physical/total_n/{field}"],
            sum(terms),
            terms,
            failures,
        )
    for field in ("hdg_tau_inward", "residual"):
        terms = tuple(values[f"discrete/{eq}/{field}"] for eq in EQUATIONS[:2])
        _identity(
            f"discrete/total_n/{field}",
            values[f"discrete/total_n/{field}"],
            sum(terms),
            terms,
            failures,
        )
    _component_sum(values, "physical/n/volume", "physical/n/volume_components", failures)
    _component_sum(
        values,
        "physical/n/boundary_physical_inward",
        "physical/n/boundary_components_inward",
        failures,
    )
    _component_sum(values, "physical/n_n/volume", "physical/n_n/volume_components", failures)
    _component_sum(
        values,
        "physical/n_n/boundary_physical_inward",
        "physical/n_n/boundary_components_inward",
        failures,
    )
    for reaction in ("ionization", "recombination"):
        terms = tuple(
            values[f"physical/{eq}/volume_components/{reaction}"]
            for eq in EQUATIONS[:2]
        )
        _identity(f"physical/particle_exchange/{reaction}", sum(terms), 0.0, terms, failures)


def _check_bc(values: dict[str, float], failures: list[str]) -> None:
    density_terms = (values["bc/n/diffusion_inward"], values["bc/n/hdg_tau_inward"])
    _identity("bc/n/residual", values["bc/n/residual"], sum(density_terms), density_terms, failures)
    source_terms = tuple(_values_under(values, "bc/n_n/source_components"))
    source_expected = sum(
        -value if name.endswith("/pump_sink") else value
        for name, value in _items_under(values, "bc/n_n/source_components")
    )
    _identity(
        "bc/n_n/imposed_source_inward",
        values["bc/n_n/imposed_source_inward"],
        source_expected,
        source_terms,
        failures,
    )
    flux_terms = tuple(_values_under(values, "bc/n_n/physical_flux_components_inward"))
    _identity(
        "bc/n_n/physical_flux_inward",
        values["bc/n_n/physical_flux_inward"],
        sum(flux_terms),
        flux_terms,
        failures,
    )
    residual_terms = (
        values["bc/n_n/imposed_source_inward"],
        -values["bc/n_n/physical_flux_inward"],
        -values["bc/n_n/hdg_tau_inward"],
    )
    _identity("bc/n_n/residual", values["bc/n_n/residual"], sum(residual_terms), residual_terms, failures)


def _check_relocated_sources(
    values: dict[str, float], context: DiagnosticContext, failures: list[str]
) -> None:
    if not context.relocated_sources:
        return
    puff = values["physical/n_n/volume_components/puff_source"]
    pump = values["physical/n_n/volume_components/pump_source"]
    if context.configured_puff is None:
        failures.append("relocated-source check requires /simulation_parameters/physics/puff")
    else:
        _identity(
            "physical/n_n/volume_components/puff_source",
            puff,
            context.configured_puff,
            (puff, context.configured_puff),
            failures,
        )
    if pump > IDENTITY_TOLERANCE * max(abs(pump), 1.0):
        failures.append("relocated neutral pump must be a non-positive volume contribution")
    for name in ("puff_source", "pump_sink"):
        path = f"bc/n_n/source_components/{name}"
        _identity(path, values[path], 0.0, (values[path],), failures)


def parse_terminal_history(path: Path) -> list[dict[str, Any]]:
    history: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    equation = ""
    section = ""
    bc_equation = ""
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
                "values": {},
            }
            equation = section = bc_equation = ""
            continue
        if current is None:
            continue
        if stripped.startswith("Physical "):
            equation = stripped.split()[1]
            section = "physical"
            continue
        if stripped == "Discrete":
            section = "discrete"
            continue
        if stripped.startswith("Independent boundary-condition checks"):
            section = "bc"
            equation = ""
            continue
        if section == "bc" and stripped.startswith(("n:", "n_n:")):
            bc_equation = stripped.split(":", 1)[0]
            continue
        match = VALUE_LINE.match(line)
        if match is None:
            continue
        label = match.group("label").strip()
        number = float(match.group("value").replace("D", "E").replace("d", "e"))
        if section == "physical" and equation in EQUATIONS:
            field = _terminal_physical_field(label)
            if field:
                current["values"][f"physical/{equation}/{field}"] = number
        elif section == "discrete" and equation in EQUATIONS:
            field = {"HDG tau inward": "hdg_tau_inward", "residual": "residual"}.get(label)
            if field:
                current["values"][f"discrete/{equation}/{field}"] = number
        elif section == "bc" and bc_equation:
            path = _terminal_bc_path(bc_equation, label)
            if path:
                current["values"][path] = number
    if current is not None:
        history.append(current)
    return history


def _terminal_physical_field(label: str) -> str | None:
    if label.startswith("content ["):
        return "content"
    return {
        "temporal [particles/s]": "temporal",
        "volume": "volume",
        "boundary physical inward": "boundary_physical_inward",
        "physical imbalance": "physical_imbalance",
    }.get(label)


def _terminal_bc_path(equation: str, label: str) -> str | None:
    fields = {
        "n": {
            "diffusion inward": "diffusion_inward",
            "HDG tau inward": "hdg_tau_inward",
            "residual": "residual",
        },
        "n_n": {
            "recycling parallel inward": (
                "source_components/recycling_parallel_inward"
            ),
            "recycling diffusion inward": (
                "source_components/recycling_diffusion_inward"
            ),
            "recycling pinch inward": "source_components/recycling_pinch_inward",
            "puff source": "source_components/puff_source",
            "pump sink (subtracted)": "source_components/pump_sink",
            "imposed source inward": "imposed_source_inward",
            "limited diffusion": (
                "physical_flux_components_inward/limited_diffusion"
            ),
            "limited pressure": (
                "physical_flux_components_inward/limited_pressure"
            ),
            "neutral-gamma convection": (
                "physical_flux_components_inward/neutral_gamma_convection"
            ),
            "physical flux inward": "physical_flux_inward",
            "HDG tau inward": "hdg_tau_inward",
            "residual": "residual",
        },
    }
    field = fields.get(equation, {}).get(label)
    return f"bc/{equation}/{field}" if field else None


def _check_terminal_values(
    block: dict[str, Any], values: dict[str, float], failures: list[str]
) -> None:
    printed = block["values"]
    expected = {
        *(f"physical/{eq}/{field}" for eq in EQUATIONS for field in (
            "content", "temporal", "volume", "boundary_physical_inward",
            "physical_imbalance",
        )),
        *(f"discrete/{eq}/{field}" for eq in EQUATIONS for field in (
            "hdg_tau_inward", "residual",
        )),
        "bc/n/diffusion_inward",
        "bc/n/hdg_tau_inward",
        "bc/n/residual",
        "bc/n_n/imposed_source_inward",
        "bc/n_n/physical_flux_inward",
        "bc/n_n/hdg_tau_inward",
        "bc/n_n/residual",
        *(name for name in values if name.startswith("bc/n_n/source_components/")),
        *(
            name
            for name in values
            if name.startswith("bc/n_n/physical_flux_components_inward/")
        ),
    }
    failures.extend(f"terminal block is missing {name}" for name in sorted(expected - printed.keys()))
    for name in sorted(expected & printed.keys()):
        exact = values[name]
        shown = printed[name]
        difference = abs(exact - shown) / max(abs(exact), abs(shown), 1.0)
        if not math.isfinite(shown) or difference > TERMINAL_TOLERANCE:
            failures.append(
                f"terminal/HDF5 mismatch for {name}: relative difference "
                f"{difference:.6e} exceeds {TERMINAL_TOLERANCE:.6e}"
            )


def _component_sum(
    values: dict[str, float], aggregate: str, prefix: str, failures: list[str]
) -> None:
    terms = tuple(_values_under(values, prefix))
    _identity(aggregate, values[aggregate], sum(terms), terms, failures)


def _items_under(values: dict[str, float], prefix: str) -> list[tuple[str, float]]:
    return [(name, value) for name, value in values.items() if name.startswith(prefix + "/")]


def _values_under(values: dict[str, float], prefix: str) -> list[float]:
    return [value for _, value in _items_under(values, prefix)]


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
        return float(value)
    except (TypeError, ValueError):
        return None


def _optional_value(handle: h5py.File, path: str) -> float | str | None:
    dataset = handle.get(path)
    if not isinstance(dataset, h5py.Dataset) or dataset.size != 1:
        return None
    return _dataset_value(dataset)


def _string_values(dataset: h5py.Dataset | h5py.Group | None) -> tuple[str, ...]:
    if not isinstance(dataset, h5py.Dataset):
        return ()
    values = dataset[()]
    values = values.reshape(-1) if hasattr(values, "reshape") else (values,)
    return tuple(
        value.decode("utf-8").strip(" \x00") if isinstance(value, bytes) else value.strip(" \x00")
        for value in values
        if isinstance(value, (bytes, str))
    )


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
