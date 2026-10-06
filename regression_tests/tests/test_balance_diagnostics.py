"""Diagnostic output behavior; deliberately no physical-identity fixture."""

import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from regression_tests.diagnostics import check_output, compare_outputs
from regression_tests.tests.fixtures.harness import run_command


@pytest.fixture
def output(tmp_path):
    def write(mode="equations"):
        solution, terminal = tmp_path / f"{mode}.h5", tmp_path / "stdout.log"
        with h5py.File(solution, "w") as handle:
            handle["simulation_parameters/switches/balance_diagnostics_mode"] = mode
            contents = {"n": 12., "n_n": 3., "nu": 7., "nEi": 9., "nEe": 8., "total_n": 15., "total_E": 17.}
            lines = [f"Balance diagnostics ({mode})"] if mode != "off" else []
            if mode in ("summary", "detailed"):
                lines.append("  Content")
            for name, value in contents.items():
                family = "particles" if name in ("n", "n_n", "total_n") else "momentum" if name == "nu" else "plasma_energy"
                content_unit, rate_unit = {"particles": ("particles", "particles/s"), "momentum": ("kg m s^-1", "N"), "plasma_energy": ("J", "W")}[family]
                if mode == "off":
                    continue
                if mode == "summary":
                    handle[f"diagnostics/summary/content/{family}/{name}"] = value
                    units = f"diagnostics/summary/content/{family}/units"
                    if units not in handle:
                        handle[units] = content_unit
                    if name in ("total_n", "total_E"):
                        handle[f"diagnostics/summary/balances/{name}/physical_imbalance"] = 4.
                        handle[f"diagnostics/summary/balances/{name}/units"] = rate_unit
                else:
                    root = f"diagnostics/equations/{name}"
                    for field, number in {"content": value, "physical/imbalance": 4., "discrete/residual": 6.}.items():
                        handle[f"{root}/{field}"] = number
                    for field, unit in {"content_units": content_unit, "rate_units": rate_unit, "physical/units": rate_unit, "discrete/units": rate_unit}.items():
                        handle[f"{root}/{field}"] = unit
                    if name != "total_n":
                        handle[f"{root}/bc/residual"] = 2.
                        handle[f"{root}/bc/units"] = rate_unit
                    if mode == "detailed":
                        for field in ("physical/temporal", "physical/volume", "physical/boundary_inward", "discrete/equation_boundary_inward", "discrete/tau_stabilization_inward", "discrete/numerical_boundary_inward"):
                            handle[f"{root}/{field}"] = 5.
                printed = {"total_n": "n+n_n", "total_E": "nEi+nEe"}.get(name, name) if mode != "equations" else name
                lines.append(f"    {printed} {value:.2E} 4.00E+00 6.00E+00")
        terminal.write_text("\n".join(lines) + "\n")
        return solution, terminal
    return write


@pytest.mark.parametrize("mode", ["off", "summary", "equations", "detailed"])
def test_formats_and_requested_mode(output, mode):
    solution, terminal = output(mode)
    assert check_output(solution, terminal, mode)["status"] == "passed"
    other = "off" if mode != "off" else "detailed"
    assert check_output(solution, terminal, other)["status"] == "failed"


def test_missing_nonfinite_units_and_terminal_failures(output):
    solution, terminal = output()
    with h5py.File(solution, "r+") as handle:
        del handle["diagnostics/equations/n/content"]
        handle["diagnostics/equations/n_n/content"][()] = np.nan
        handle["diagnostics/equations/nu/rate_units"][()] = "wrong"
    failures = check_output(solution, terminal)["failures"]
    assert any("missing scalar diagnostic: n/content" in item for item in failures)
    assert any("non-finite" in item for item in failures)
    assert any("units: nu/rate_units" in item for item in failures)
    solution, terminal = output()
    terminal.write_text(terminal.read_text().replace("1.20E+01", "1.30E+01"))
    assert "terminal/HDF5 content mismatch: n" in check_output(solution, terminal)["failures"]


def test_parallel_scaling_and_changed_scalar(output, tmp_path):
    off, _ = output("off")
    candidate = tmp_path / "candidate.h5"
    candidate.write_bytes(off.read_bytes())
    assert compare_outputs(off, candidate)["status"] == "skipped"
    solution, _ = output()
    assert compare_outputs(off, solution)["status"] == "failed"
    candidate.write_bytes(solution.read_bytes())
    for path in (solution, candidate):
        with h5py.File(path, "r+") as handle:
            # Large particle content must not hide a rate discrepancy.
            handle["diagnostics/equations/n/content"][()] = 1e20
            handle["diagnostics/equations/n/bc/residual"][()] = 0.
    with h5py.File(candidate, "r+") as handle:
        handle["diagnostics/equations/n/bc/residual"][()] += 1e-10
    assert compare_outputs(solution, candidate)["status"] == "passed"
    with h5py.File(candidate, "r+") as handle:
        handle["diagnostics/equations/n/physical/imbalance"][()] += .01
    assert "diagnostic differs across layouts: n/physical/imbalance" in compare_outputs(solution, candidate)["failures"]
    for path in (solution, candidate):
        with h5py.File(path, "r+") as handle:
            del handle["diagnostics/equations"]
            handle.create_group("diagnostics/equations")
    assert compare_outputs(solution, candidate)["status"] == "failed"


def test_known_puff_and_cli_attach_to_existing_output(output, tmp_path):
    solution, terminal = output("detailed")
    with h5py.File(solution, "r+") as handle:
        handle["simulation_parameters/switches/neutral_wall_sources_in_elements"] = 1
        handle["simulation_parameters/physics/puff"] = 100.
        handle["diagnostics/equations/n_n/physical/volume_components/puff"] = 100.
        handle["diagnostics/equations/n_n/bc/source_components/puff"] = 0.
        handle["diagnostics/equations/n_n/bc/source_components/pump"] = 0.
    assert check_output(solution, terminal)["status"] == "passed"
    (tmp_path / "run_metadata.json").write_text(json.dumps({"hdf5_outputs": [solution.name]}))
    summary = tmp_path / "suite_summary.json"
    summary.write_text(json.dumps({"results": [{"workflow_id": "example", "layout_id": "serial", "run_status": "completed", "run_directory": str(tmp_path)}]}))
    result = run_command("compare", "--suite", "--diagnostics", str(summary))
    assert result.returncode == 0, result.stderr
    assert "balance diagnostics passed:" in result.stdout
    with h5py.File(solution, "r+") as handle:
        handle["diagnostics/equations/n_n/physical/volume_components/puff"][()] = 90.
    assert run_command("compare", "--suite", "--diagnostics", str(summary)).returncode == 1


@pytest.mark.parametrize("mode", ["summary", "detailed"])
def test_neutral_wall_absorption_contract(output, mode):
    solution, terminal = output(mode)
    root = "diagnostics/summary/" if mode == "summary" else "diagnostics/equations/"
    source = root + ("external_sources/particles/" if mode == "summary" else "n_n/bc/source_components/")
    physical = root + "n_n/physical/boundary_components_inward/"
    with h5py.File(solution, "r+") as handle:
        for name in ("recycling_neutral", "recycling_neutral_pump"):
            handle["simulation_parameters/physics/" + name] = .99
        handle[source + "neutral_wall_absorption"] = 6.
        handle[source + "units"] = "particles/s"
        if mode == "detailed":
            handle[physical + "neutral_wall_absorption"] = -6.
            handle[physical + "units"] = "particles/s"
    terminal.write_text(terminal.read_text() + ("neutral wall absorption -6.00E+00\n" if mode == "detailed" else "")
                        + "neutral wall absorption 6.00E+00\n")
    assert check_output(solution, terminal)["status"] == "passed"
    with h5py.File(solution, "r+") as handle:
        if mode == "detailed":
            handle[physical + "neutral_wall_absorption"][()] = 6.
        else:
            del handle[source + "neutral_wall_absorption"]
    assert check_output(solution, terminal)["status"] == "failed"


def test_unit_recycling_requires_zero_wall_absorption(output):
    solution, terminal = output("summary")
    with h5py.File(solution, "r+") as handle:
        for name in ("recycling_neutral", "recycling_neutral_pump"):
            handle["simulation_parameters/physics/" + name] = 1.
        handle["diagnostics/summary/external_sources/particles/neutral_wall_absorption"] = 1.
        handle["diagnostics/summary/external_sources/particles/units"] = "particles/s"
    terminal.write_text(terminal.read_text() + "neutral wall absorption 1.00E+00\n")
    assert "neutral wall absorption is nonzero with unit recycling" in check_output(solution, terminal)["failures"]
