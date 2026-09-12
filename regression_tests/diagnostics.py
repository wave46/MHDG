"""Small output checks and parallel comparisons for solver diagnostics."""

from pathlib import Path
import re

import h5py
import numpy as np

from .compare_common import select_candidate
from .documents import load_json
from .support import HarnessError

UNITS = {
    "n": ("particles", "particles/s"), "n_n": ("particles", "particles/s"),
    "nu": ("kg m s^-1", "N"), "nEi": ("J", "W"), "nEe": ("J", "W"),
    "total_n": ("particles", "particles/s"), "total_E": ("J", "W"),
}
CONTENT = {name: f"content/{family}/{name}" for family, names in (
    ("particles", ("n", "n_n", "total_n")), ("momentum", ("nu",)),
    ("plasma_energy", ("nEi", "nEe", "total_E")),
) for name in names}
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
TERMINAL_RTOL = 5.1e-3
PARALLEL_RTOL = 1e-9
PARALLEL_ATOL = 1e-12  # Floor in each dataset's reported physical units.


def _text(item):
    if not isinstance(item, h5py.Dataset) or item.size != 1:
        return None
    value = np.asarray(item[()]).reshape(-1)[0]
    return value.decode().strip(" \x00") if isinstance(value, bytes) else str(value).strip(" \x00")


def _read(handle, failures):
    mode = _text(handle.get("simulation_parameters/switches/balance_diagnostics_mode"))
    values, texts = {}, {}
    if mode in ("summary", "equations", "detailed"):
        root = handle.get("diagnostics/summary" if mode == "summary" else "diagnostics/equations")
        if not isinstance(root, h5py.Group):
            failures.append(f"missing diagnostic output for {mode}")
        else:
            def collect(name, item):
                if isinstance(item, h5py.Dataset):
                    if item.dtype.kind in "iuf":
                        if not np.isfinite(item[()]).all():
                            failures.append(f"non-finite diagnostic: {name}")
                        if item.size == 1:
                            values[name] = float(np.asarray(item[()]).reshape(-1)[0])
                    elif item.size == 1:
                        texts[name] = _text(item)
            root.visititems(collect)
    return mode, values, texts


def check_output(solution, terminal, expected_mode=None):
    """Check output usability, representative printed content and known puff input."""
    failures = []
    with h5py.File(solution) as handle:
        mode, values, texts = _read(handle, failures)
        if mode is None and expected_mode is None:
            return None  # No diagnostic mode recorded or requested.
        if expected_mode is not None and mode != expected_mode:
            failures.append(f"diagnostics mode {mode!r}, expected {expected_mode!r}")
        if mode not in ("off", "summary", "equations", "detailed"):
            failures.append(f"unsupported diagnostics mode: {mode!r}")
        log = terminal.read_text(errors="replace")
        if mode == "off":
            if any(path in handle for path in ("diagnostics/summary", "diagnostics/equations")) or "Balance diagnostics (" in log:
                failures.append("diagnostics output is present with mode off")
        elif mode in ("summary", "equations", "detailed"):
            unexpected = "diagnostics/equations" if mode == "summary" else "diagnostics/summary"
            if unexpected in handle or (mode == "equations" and "n/physical/volume" in values):
                failures.append(f"unexpected diagnostic detail for mode {mode}")
            required, units = set(), {}
            for name, (content_unit, rate_unit) in UNITS.items():
                if mode == "summary":
                    path = CONTENT[name]
                    required.add(path)
                    units[path.rsplit("/", 1)[0] + "/units"] = content_unit
                else:
                    required.update(f"{name}/{field}" for field in ("content", "physical/imbalance", "discrete/residual"))
                    units.update({
                        f"{name}/content_units": content_unit, f"{name}/rate_units": rate_unit,
                        f"{name}/physical/units": rate_unit, f"{name}/discrete/units": rate_unit,
                    })
                    if name != "total_n":
                        required.add(f"{name}/bc/residual")
                        units[f"{name}/bc/units"] = rate_unit
                    if mode == "detailed":
                        required.update(f"{name}/physical/{field}" for field in ("temporal", "volume", "boundary_inward"))
                        required.update(f"{name}/discrete/{field}" for field in (
                            "equation_boundary_inward", "tau_stabilization_inward", "numerical_boundary_inward",
                        ))
            if mode == "summary":
                for name in ("total_n", "total_E"):
                    required.add(f"balances/{name}/physical_imbalance")
                    units[f"balances/{name}/units"] = UNITS[name][1]
            if mode == "detailed":
                names = handle.get("simulation_parameters/physics/conservative_variable_names")
                names = [] if names is None else np.asarray(names[()]).reshape(-1)
                if "Gamman" in [v.decode().strip() if isinstance(v, bytes) else str(v).strip() for v in names]:
                    required.update(("n_n/discrete/boundary_components_inward/neutral_gamma_convection",
                                     "n_n/bc/physical_flux_components_inward/neutral_gamma_convection"))
                relocated = handle.get("simulation_parameters/switches/neutral_wall_sources_in_elements")
                if relocated is not None and np.asarray(relocated[()]).item() == 1:
                    puff = handle.get("simulation_parameters/physics/puff")
                    actual = values.get("n_n/physical/volume_components/puff", np.nan)
                    expected = float(np.asarray(puff[()]).item()) if puff is not None else np.nan
                    if not np.isfinite(expected) or not np.isfinite(actual) or abs(actual - expected) > 1e-12 * max(abs(expected), 1.):
                        failures.append("integrated relocated puff differs from configured puff")
                    for field in ("puff", "pump"):
                        value = values.get(f"n_n/bc/source_components/{field}", np.nan)
                        if not np.isfinite(value) or abs(value) > 1e-12:
                            failures.append(f"relocated wall {field} is missing or nonzero")
            failures.extend(f"missing scalar diagnostic: {name}" for name in sorted(required - values.keys()))
            failures.extend(f"wrong or missing units: {name}" for name, unit in units.items() if texts.get(name) != unit)
            _check_terminal(log, mode, values, failures)
    return {"solution": str(solution), "mode": mode, "status": "failed" if failures else "passed", "failures": failures}


def _check_terminal(log, mode, values, failures):
    block = log.rsplit("Balance diagnostics (", 1)[-1]
    if not block.startswith(mode + ")"):
        failures.append(f"missing final {mode} terminal block")
        return
    if mode == "equations":
        pairs = re.findall(rf"^\s*(n|n_n|nu|nEi|nEe|total_n|total_E)\s+({NUMBER})", block, re.M)
    else:
        content = block.split("  Content\n", 1)[-1]
        content = re.split(r"\n  [A-Z]", content, maxsplit=1)[0]
        pairs = re.findall(rf"(?<!\S)(nEi\+nEe|n\+n_n|nEi|nEe|n_n|nu|n)\s+({NUMBER})(?=\s|$)", content)
    aliases = {"n+n_n": "total_n", "nEi+nEe": "total_E"}
    printed = {
        aliases.get(name, name): float(value.replace("D", "E").replace("d", "e"))
        for name, value in pairs
    }
    for name in UNITS:
        actual = printed.get(name, np.nan)
        expected = values.get(CONTENT[name] if mode == "summary" else f"{name}/content", np.nan)
        if not np.isfinite(actual) or not np.isfinite(expected) or abs(actual - expected) > TERMINAL_RTOL * max(abs(actual), abs(expected), 1.):
            failures.append(f"terminal/HDF5 content mismatch: {name}")


def check_suite(summary, *, required=True):
    """Attach checks to existing outputs; no solver launch or duplicate history."""
    if isinstance(summary, Path):
        summary = load_json(summary, "suite summary")
    reports, failures = [], []
    for result in summary.get("results", []):
        label = f"{result.get('workflow_id')}/{result.get('layout_id')}"
        if result.get("run_status") != "completed":
            failures.append(f"{label}: solver run did not complete")
            continue
        try:
            run = Path(result["run_directory"])
            metadata = load_json(run / "run_metadata.json", "run metadata")
            for stage in metadata.get("stages") or [{"run_directory": str(run), "status": "completed"}]:
                directory = Path(stage["run_directory"])
                if stage.get("status") != "completed":
                    failures.append(f"{label}: incomplete diagnostic stage")
                    continue
                plan_path = directory / "run_plan.json"
                plan = load_json(plan_path, "run plan") if plan_path.exists() else {}
                expected = plan.get("parameter_overrides", {}).get("balance_diagnostics_mode")
                parameters = directory / "param.txt"
                if parameters.exists():
                    match = re.search(r"(?im)^\s*balance_diagnostics_mode\s*=\s*['\"]([^'\"]+)['\"]", parameters.read_text())
                    if match:
                        expected = match.group(1)
                if stage.get("selected_hdf5"):
                    solution = Path(stage["selected_hdf5"])
                    if not solution.is_absolute():
                        solution = directory / solution
                else:
                    solution = select_candidate(directory, metadata)
                report = check_output(solution, directory / "stdout.log", expected)
                if report:
                    reports.append({"run": label, **report})
                    failures.extend(f"{label}: {failure}" for failure in report["failures"])
        except (HarnessError, OSError, ValueError) as exc:
            failures.append(f"{label}: {exc}")
    if required and not reports:
        failures.append("no diagnostic outputs or explicit off-mode evidence found")
    return {"status": "failed" if failures else "passed", "outputs": reports, "failures": failures}


def compare_outputs(reference, candidate):
    """Compare scalar diagnostics across matching-state parallel runs."""
    failures = []
    with h5py.File(reference) as first, h5py.File(candidate) as second:
        left_mode, left, _ = _read(first, failures)
        right_mode, right, _ = _read(second, failures)
    if left_mode != right_mode or left.keys() != right.keys():
        failures.append("diagnostic mode or scalar dataset selection differs across layouts")
    scales = {}
    def group(name):
        if left_mode == "summary":
            return name  # Summary exposes individual aggregates, not equation rate components.
        return (name.split("/")[0], "content" if name.endswith("/content") else "rate")
    for name in sorted(left.keys() & right.keys()):
        key = group(name)
        scales[key] = max(scales.get(key, 0.), abs(left[name]), abs(right[name]))
    for name in sorted(left.keys() & right.keys()):
        if abs(left[name] - right[name]) > max(PARALLEL_ATOL, PARALLEL_RTOL * scales[group(name)]):
            failures.append(f"diagnostic differs across layouts: {name}")
    return {"status": "failed" if failures else "passed", "failures": failures,
            "relative_tolerance": PARALLEL_RTOL, "absolute_floor": PARALLEL_ATOL}
