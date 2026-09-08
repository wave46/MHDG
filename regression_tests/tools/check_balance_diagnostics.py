#!/usr/bin/env python3
"""Validate the detailed equation-oriented balance diagnostics contract."""

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
TERMINAL_TOLERANCE = 5.1e-3
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
CONTENT_PAIR = re.compile(
    rf"(?<!\S)(?P<equation>nEi\+nEe|n\+n_n|nEi|nEe|n_n|nu|n)\s+"
    rf"(?P<value>{NUMBER})(?=\s|$)"
)
AGGREGATE_HEADING = re.compile(
    r"^\s{4}(?P<equation>nEi\+nEe|n\+n_n|nEi|nEe|n_n|nu|n)\s+"
    r"\[[^]]+\](?:\s+\(derived\))?\s*$"
)
NUMERIC_ROW = re.compile(rf"^\s+(?P<values>{NUMBER}(?:\s+{NUMBER})+)\s*$")
TIME_ITERATION = re.compile(r"Time iteration\s*=\s*(\d+)")
NEWTON_ITERATION = re.compile(r"NR iteration:\s*(\d+)")
BLOCK_START = "Balance diagnostics (detailed)"
PRINTED_EQUATIONS = {"n+n_n": "total_n", "nEi+nEe": "total_E"}


@dataclass(frozen=True)
class EquationSpec:
    name: str
    content_units: str
    rate_units: str
    volume_components: tuple[str, ...]
    physical_boundary_components: tuple[str, ...]
    equation_boundary_components: tuple[str, ...]


PRIMARY_EQUATIONS = (
    EquationSpec(
        "n",
        "particles",
        "particles/s",
        ("ionization", "recombination", "prescribed_source"),
        ("parallel_convection", "pinch"),
        ("parallel_convection", "diffusion", "pinch"),
    ),
    EquationSpec(
        "nu",
        "kg m s^-1",
        "N",
        (
            "ionization",
            "recombination",
            "charge_exchange",
            "pressure_divergence",
            "prescribed_source",
        ),
        (),
        ("convection", "diffusion", "pinch"),
    ),
    EquationSpec(
        "nEi",
        "J",
        "W",
        (
            "ionization",
            "recombination",
            "charge_exchange",
            "parallel_electric_work",
            "temperature_exchange",
            "prescribed_source",
        ),
        ("sheath", "pinch"),
        ("convection", "diffusion", "parallel_conduction", "pinch"),
    ),
    EquationSpec(
        "nEe",
        "J",
        "W",
        (
            "ionization",
            "recombination",
            "radiation",
            "ohmic",
            "parallel_electric_work",
            "temperature_exchange",
            "prescribed_source",
        ),
        ("sheath", "pinch"),
        ("convection", "diffusion", "parallel_conduction", "pinch"),
    ),
    EquationSpec(
        "n_n",
        "particles",
        "particles/s",
        ("ionization", "recombination", "prescribed_source", "puff", "pump"),
        (
            "recycling_parallel",
            "recycling_diffusion",
            "recycling_pinch",
            "puff",
            "pump",
        ),
        ("limited_diffusion", "limited_pressure"),
    ),
)
DERIVED_EQUATIONS = {
    "total_n": ("n", "n_n", "particles", "particles/s"),
    "total_E": ("nEi", "nEe", "J", "W"),
}
PRIMARY_NAMES = tuple(spec.name for spec in PRIMARY_EQUATIONS)
EQUATIONS = (*PRIMARY_NAMES, *DERIVED_EQUATIONS)
MOMENTUM_BC_COMPONENTS = (
    "perpendicular_diffusion_inward",
    "split_diffusion_inward",
    "tau_stabilization_inward",
)
ENERGY_BC_COMPONENTS = (
    "perpendicular_diffusion_inward",
    "split_diffusion_inward",
    "parallel_conduction_inward",
    "sheath_minus_bulk_inward",
    "tau_stabilization_inward",
)


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
        if result.get("run_status") != "completed":
            report = {
                "status": "failed",
                "failures": [f"solver run did not complete: {workflow}/{layout}"],
                "stages": [],
            }
        else:
            run_directory = require_directory(
                Path(result["run_directory"]),
                f"balance diagnostics run {workflow}/{layout}",
            )
            report = check_run(run_directory)
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
        reports = [_check_staged_output(stage) for stage in stages]
    else:
        try:
            solution = select_candidate(run_directory, metadata)
        except ComparisonError as exc:
            failures.append(f"warm: {exc}")
        else:
            reports = [
                check_stage(
                    "warm",
                    solution,
                    require_file(run_directory / "stdout.log", "warm diagnostic log"),
                )
            ]
    for report in reports:
        failures.extend(
            f"{report['stage_id']}: {failure}" for failure in report["failures"]
        )
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
        str(stage_id),
        require_file(solution, f"diagnostic stage solution {stage_id}"),
        require_file(stage_directory / "stdout.log", f"stage log {stage_id}"),
    )


def check_stage(stage_id: str, solution: Path, terminal: Path) -> dict[str, Any]:
    failures: list[str] = []
    values, context = _read_hdf5(solution, failures)
    if values:
        _check_equation_identities(values, failures)
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
            puff = _optional_value(handle, "simulation_parameters/physics/puff")
            context = DiagnosticContext(
                neutral_gamma="Gamman" in names,
                relocated_sources=relocated == 1.0,
                configured_puff=puff if isinstance(puff, float) else None,
            )
            root = handle.get("diagnostics/equations")
            if not isinstance(root, h5py.Group):
                failures.append("required group is missing: /diagnostics/equations")
                return {}, context

            def collect(name: str, item: h5py.Group | h5py.Dataset) -> None:
                if not isinstance(item, h5py.Dataset) or item.size != 1:
                    return
                value = _dataset_value(item)
                if isinstance(value, str):
                    texts[name] = value
                elif isinstance(value, float):
                    values[name] = value

            root.visititems(collect)
    except OSError as exc:
        failures.append(f"cannot read HDF5 solution: {exc}")
        return {}, context

    required = _required_paths(context)
    missing = sorted(required - values.keys())
    failures.extend(
        f"required scalar dataset is missing: /diagnostics/equations/{name}"
        for name in missing
    )
    nonfinite = sorted(
        name for name in required & values.keys() if not math.isfinite(values[name])
    )
    failures.extend(
        f"required scalar dataset is non-finite: /diagnostics/equations/{name}"
        for name in nonfinite
    )
    _check_units(texts, failures)
    return (values, context) if not missing and not nonfinite else ({}, context)


def _required_paths(context: DiagnosticContext) -> set[str]:
    required = {
        *(
            f"{equation}/{field}"
            for equation in EQUATIONS
            for field in (
                "content",
                "physical/temporal",
                "physical/volume",
                "physical/boundary_inward",
                "physical/imbalance",
                "discrete/equation_boundary_inward",
                "discrete/tau_stabilization_inward",
                "discrete/numerical_boundary_inward",
                "discrete/residual",
            )
        ),
        *(
            f"{spec.name}/physical/volume_components/{component}"
            for spec in PRIMARY_EQUATIONS
            for component in spec.volume_components
        ),
        *(
            f"{spec.name}/physical/boundary_components_inward/{component}"
            for spec in PRIMARY_EQUATIONS
            for component in spec.physical_boundary_components
        ),
        *(
            f"{spec.name}/discrete/boundary_components_inward/{component}"
            for spec in PRIMARY_EQUATIONS
            for component in spec.equation_boundary_components
        ),
        "n/physical/exchange/charge_exchange_rate",
        "n/bc/diffusion_inward",
        "n/bc/tau_stabilization_inward",
        "n/bc/residual",
        *(f"nu/bc/{component}" for component in MOMENTUM_BC_COMPONENTS),
        "nu/bc/residual",
        *(
            f"{equation}/bc/{component}"
            for equation in ("nEi", "nEe", "total_E")
            for component in ENERGY_BC_COMPONENTS
        ),
        "nEi/bc/residual",
        "nEe/bc/residual",
        "total_E/bc/residual",
        "n_n/bc/imposed_source_inward",
        "n_n/bc/physical_flux_inward",
        "n_n/bc/tau_stabilization_inward",
        "n_n/bc/residual",
        "n_n/bc/source_components/recycling_parallel_inward",
        "n_n/bc/source_components/recycling_diffusion_inward",
        "n_n/bc/source_components/recycling_pinch_inward",
        "n_n/bc/source_components/puff",
        "n_n/bc/source_components/pump",
        "n_n/bc/physical_flux_components_inward/limited_diffusion",
        "n_n/bc/physical_flux_components_inward/limited_pressure",
    }
    if context.neutral_gamma:
        required.update(
            {
                "n_n/discrete/boundary_components_inward/neutral_gamma_convection",
                "n_n/bc/physical_flux_components_inward/neutral_gamma_convection",
            }
        )
    return required


def _check_units(texts: dict[str, str], failures: list[str]) -> None:
    specs = {
        **{spec.name: (spec.content_units, spec.rate_units) for spec in PRIMARY_EQUATIONS},
        **{name: units[2:] for name, units in DERIVED_EQUATIONS.items()},
    }
    for equation, (content_units, rate_units) in specs.items():
        expected = {
            f"{equation}/content_units": content_units,
            f"{equation}/rate_units": rate_units,
            f"{equation}/physical/units": rate_units,
            f"{equation}/discrete/units": rate_units,
        }
        if equation in ("n", "nu", "nEi", "nEe", "n_n", "total_E"):
            expected[f"{equation}/bc/units"] = rate_units
        for name, value in expected.items():
            if texts.get(name) != value:
                failures.append(
                    f"/diagnostics/equations/{name} is {texts.get(name)!r}, "
                    f"expected {value!r}"
                )


def _check_equation_identities(
    values: dict[str, float], failures: list[str]
) -> None:
    for equation in EQUATIONS:
        prefix = equation
        temporal = values[f"{prefix}/physical/temporal"]
        volume = values[f"{prefix}/physical/volume"]
        boundary = values[f"{prefix}/physical/boundary_inward"]
        equation_boundary = values[f"{prefix}/discrete/equation_boundary_inward"]
        tau = values[f"{prefix}/discrete/tau_stabilization_inward"]
        numerical = values[f"{prefix}/discrete/numerical_boundary_inward"]
        _identity(
            f"{prefix}/physical/imbalance",
            values[f"{prefix}/physical/imbalance"],
            volume + boundary - temporal,
            (volume, boundary, temporal),
            failures,
        )
        _identity(
            f"{prefix}/discrete/numerical_boundary_inward",
            numerical,
            equation_boundary + tau,
            (equation_boundary, tau),
            failures,
        )
        _identity(
            f"{prefix}/discrete/residual",
            values[f"{prefix}/discrete/residual"],
            volume + numerical - temporal,
            (volume, numerical, temporal),
            failures,
        )

    additive_fields = (
        "content",
        "physical/temporal",
        "physical/volume",
        "physical/boundary_inward",
        "physical/imbalance",
        "discrete/equation_boundary_inward",
        "discrete/tau_stabilization_inward",
        "discrete/numerical_boundary_inward",
        "discrete/residual",
    )
    for total, (first, second, _, _) in DERIVED_EQUATIONS.items():
        for field in additive_fields:
            terms = (values[f"{first}/{field}"], values[f"{second}/{field}"])
            _identity(f"{total}/{field}", values[f"{total}/{field}"], sum(terms), terms, failures)

    for spec in PRIMARY_EQUATIONS:
        _component_sum(
            values,
            f"{spec.name}/physical/volume",
            f"{spec.name}/physical/volume_components",
            failures,
        )
        if spec.physical_boundary_components:
            _component_sum(
                values,
                f"{spec.name}/physical/boundary_inward",
                f"{spec.name}/physical/boundary_components_inward",
                failures,
            )
        _component_sum(
            values,
            f"{spec.name}/discrete/equation_boundary_inward",
            f"{spec.name}/discrete/boundary_components_inward",
            failures,
        )
    _identity(
        "nu physical/equation boundary agreement",
        values["nu/physical/boundary_inward"],
        values["nu/discrete/equation_boundary_inward"],
        (values["nu/discrete/equation_boundary_inward"],),
        failures,
    )
    for reaction in ("ionization", "recombination"):
        terms = tuple(
            values[f"{equation}/physical/volume_components/{reaction}"]
            for equation in ("n", "n_n")
        )
        _identity(f"particle exchange/{reaction}", sum(terms), 0.0, terms, failures)
    for exchange in ("parallel_electric_work", "temperature_exchange"):
        terms = tuple(
            values[f"{equation}/physical/volume_components/{exchange}"]
            for equation in ("nEi", "nEe")
        )
        _identity(f"energy exchange/{exchange}", sum(terms), 0.0, terms, failures)


def _check_bc(values: dict[str, float], failures: list[str]) -> None:
    density_terms = (
        values["n/bc/diffusion_inward"],
        values["n/bc/tau_stabilization_inward"],
    )
    _identity("n/bc/residual", values["n/bc/residual"], sum(density_terms), density_terms, failures)
    _check_bc_component_sum(values, "nu", MOMENTUM_BC_COMPONENTS, failures)
    for equation in ("nEi", "nEe"):
        _check_bc_component_sum(values, equation, ENERGY_BC_COMPONENTS, failures)
    for component in ENERGY_BC_COMPONENTS:
        terms = tuple(values[f"{equation}/bc/{component}"] for equation in ("nEi", "nEe"))
        _identity(
            f"total_E/bc/{component}",
            values[f"total_E/bc/{component}"],
            sum(terms),
            terms,
            failures,
        )
    _check_bc_component_sum(values, "total_E", ENERGY_BC_COMPONENTS, failures)
    source_prefix = "n_n/bc/source_components"
    source_terms = tuple(_values_under(values, source_prefix))
    source_expected = sum(
        -value if name.endswith("/pump") else value
        for name, value in _items_under(values, source_prefix)
    )
    _identity(
        "n_n/bc/imposed_source_inward",
        values["n_n/bc/imposed_source_inward"],
        source_expected,
        source_terms,
        failures,
    )
    flux_prefix = "n_n/bc/physical_flux_components_inward"
    flux_terms = tuple(_values_under(values, flux_prefix))
    _identity(
        "n_n/bc/physical_flux_inward",
        values["n_n/bc/physical_flux_inward"],
        sum(flux_terms),
        flux_terms,
        failures,
    )
    residual_terms = (
        values["n_n/bc/imposed_source_inward"],
        -values["n_n/bc/physical_flux_inward"],
        -values["n_n/bc/tau_stabilization_inward"],
    )
    _identity("n_n/bc/residual", values["n_n/bc/residual"], sum(residual_terms), residual_terms, failures)


def _check_bc_component_sum(
    values: dict[str, float],
    equation: str,
    components: Sequence[str],
    failures: list[str],
) -> None:
    terms = tuple(values[f"{equation}/bc/{component}"] for component in components)
    _identity(
        f"{equation}/bc/residual",
        values[f"{equation}/bc/residual"],
        sum(terms),
        terms,
        failures,
    )


def _check_relocated_sources(
    values: dict[str, float], context: DiagnosticContext, failures: list[str]
) -> None:
    if not context.relocated_sources:
        return
    puff = values["n_n/physical/volume_components/puff"]
    pump = values["n_n/physical/volume_components/pump"]
    if context.configured_puff is None:
        failures.append("relocated-source check requires configured puff")
    else:
        _identity("relocated puff", puff, context.configured_puff, (puff,), failures)
    if pump > IDENTITY_TOLERANCE * max(abs(pump), 1.0):
        failures.append("relocated neutral pump must be non-positive")
    for name in ("puff", "pump"):
        path = f"n_n/bc/source_components/{name}"
        _identity(path, values[path], 0.0, (values[path],), failures)


def parse_terminal_history(path: Path) -> list[dict[str, Any]]:
    history: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    time_iteration = None
    newton_iteration = None
    section = None
    equation = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if match := TIME_ITERATION.search(line):
            time_iteration = int(match.group(1))
        if match := NEWTON_ITERATION.search(line):
            newton_iteration = int(match.group(1))
        if line.strip() == BLOCK_START:
            if current is not None:
                history.append(current)
            current = {
                "time_iteration": time_iteration,
                "newton_iteration": newton_iteration,
                "values": {},
            }
            section = None
            equation = None
            continue
        if current is None:
            continue
        if line.strip() == "Content":
            section = "content"
            equation = None
            continue
        if line.strip() == "Physical balances":
            section = "physical"
            equation = None
            continue
        if line.strip() == "Discrete equations":
            section = "discrete"
            equation = None
            continue
        if line.strip() == "Independent boundary-condition checks":
            section = "bc"
            equation = None
            continue
        if section == "content":
            for match in CONTENT_PAIR.finditer(line):
                name = PRINTED_EQUATIONS.get(
                    match.group("equation"), match.group("equation")
                )
                current["values"][f"{name}/content"] = _number(
                    match.group("value")
                )
            continue
        if section not in ("physical", "discrete"):
            continue
        if match := AGGREGATE_HEADING.match(line):
            equation = PRINTED_EQUATIONS.get(
                match.group("equation"), match.group("equation")
            )
            continue
        if equation is None or (match := NUMERIC_ROW.match(line)) is None:
            continue
        shown = tuple(_number(value) for value in match.group("values").split())
        if section == "physical" and len(shown) == 4:
            current["values"][f"{equation}/physical/imbalance"] = shown[-1]
            equation = None
        elif section == "discrete" and len(shown) == 6:
            current["values"][f"{equation}/discrete/residual"] = shown[-1]
            equation = None
    if current is not None:
        history.append(current)
    return history


def _check_terminal_values(
    block: dict[str, Any], values: dict[str, float], failures: list[str]
) -> None:
    printed = block["values"]
    expected = {
        *(f"{equation}/{field}" for equation in EQUATIONS for field in (
            "content", "physical/imbalance", "discrete/residual"
        )),
    }
    failures.extend(
        f"terminal block is missing {name}" for name in sorted(expected - printed.keys())
    )
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
    return [
        (name, value) for name, value in values.items() if name.startswith(prefix + "/")
    ]


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


def _number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _dataset_value(dataset: h5py.Dataset) -> float | str | None:
    value = dataset[()]
    if hasattr(value, "reshape"):
        value = value.reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode("utf-8").strip(" \x00")
    if isinstance(value, str):
        return value.strip(" \x00")
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
        value.decode("utf-8").strip(" \x00")
        if isinstance(value, bytes)
        else value.strip(" \x00")
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
