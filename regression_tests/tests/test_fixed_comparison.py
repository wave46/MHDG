"""Direct field/mesh contracts and the public comparison command."""

import json

import h5py
import numpy as np
import pytest

from regression_tests.compare_fixed import compare_hdf5_files, validate_solution_file
from regression_tests.compare_common import read_newton_convergence
from regression_tests.compare import compare_completed_run
from regression_tests.tests.fixtures.harness import REGRESSION_ROOT, run_command
from regression_tests.tests.fixtures.solutions import write_solution

TOLERANCES = {
    "newton_error_max": 2e-4, "mesh_coordinate_atol": 1e-12,
    "relative_l2_max": 1e-10, "normalized_linf_max": 1e-9,
}


@pytest.fixture
def files(tmp_path):
    reference, candidate = tmp_path / "reference.h5", tmp_path / "candidate.h5"
    write_solution(reference)
    write_solution(candidate)
    return reference, candidate


def test_missing_candidate_equation_names_use_reference_names(files):
    reference, candidate = files
    with h5py.File(candidate, "r+") as handle:
        del handle["simulation_parameters/physics/conservative_variable_names"]
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == "passed"
    assert report["solution"]["equation_names"] == ["rho", "Gamma"]
    assert not validate_solution_file(reference)
    assert not validate_solution_file(candidate)


@pytest.mark.parametrize("dataset,index,value,section", [
    ("solution/u", 1, 3.0, "solution"),
    ("solution/q", 0, np.nan, "solution"),
    ("mesh/T", (0, 0), 99, "mesh"),
    ("transport_1d/coefficients/d_fs", 0, 2.0, "transport_1d"),
])
def test_field_nonfinite_mesh_and_transport_failures(files, dataset, index, value, section):
    reference, candidate = files
    with h5py.File(candidate, "r+") as handle:
        handle[dataset][index] = value
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == "failed"
    assert not report[section]["passed"]
    # A finite changed field is valid data even though it fails regression;
    # bad connectivity and nonfinite fields are invalid on their own.
    assert bool(validate_solution_file(candidate)) == (dataset in {"solution/q", "mesh/T"})
    if section == "transport_1d":
        with h5py.File(candidate, "r+") as handle:
            handle[dataset][index] = np.nan
        assert any("transport_1d" in failure for failure in validate_solution_file(candidate))
    if dataset == "solution/u":
        equations = report["solution"]["datasets"]["u"]["equations"]
        assert equations["rho"]["passed"] and not equations["Gamma"]["passed"]
    elif dataset == "solution/q":
        assert not report["solution"]["datasets"]["q"]["equations"]["rho"]["finite"]
    elif dataset == "mesh/T":
        assert report["solution"]["reason"] == "mesh comparison failed"


@pytest.mark.parametrize("difference", [None, "region", "missing"])
def test_magnetic_contract(files, difference):
    reference, candidate = files
    for path in files:
        if path == candidate and difference == "missing":
            continue
        with h5py.File(path, "r+") as handle:
            magnetic = handle.create_group("magnetic")
            for name, values in {
                "topology": [b"limited"], "topology_id": [1],
                "topology_region": [1, 1, 2, 2], "rho_pol_norm": [0., .5, 1., 1.2],
                "topology_normal": [[1., 1., 1., 1.], [0., 0., 0., 0.]],
            }.items():
                magnetic[name] = np.asarray(values)
            if path == candidate and difference == "region":
                magnetic["topology_region"][2] = 3
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == ("failed" if difference else "passed")
    assert not validate_solution_file(candidate)
    if difference == "missing":
        assert "candidate is missing magnetic data" in report["failures"]
    elif difference == "region":
        assert "magnetic/topology_region differs" in report["failures"]
    else:
        with h5py.File(candidate, "r+") as handle:
            del handle["transport_1d"]  # Optional data need not exist for file validity.
        assert not validate_solution_file(candidate)
        with h5py.File(candidate, "r+") as handle:
            handle["magnetic/rho_pol_norm"][0] = np.inf
        assert any("magnetic/rho_pol_norm" in failure for failure in validate_solution_file(candidate))


@pytest.mark.parametrize("last_value,expected", [
    (None, None), ("", None), ("NaN", None), ("Infinity", None),
    ("********", None), ("1.0E309", None),
    ("1.0D-2", 1e-2), ("1.0d-5", 1e-5), ("1.0E-2", 1e-2),
])
def test_newton_acceptance_uses_last_record(tmp_path, last_value, expected):
    log = tmp_path / "stdout.log"
    log.write_text("no Newton iteration\n" if last_value is None else (
        f" Error: 1.0E-5\n\tError: {last_value}\nOutput written to file result.h5\n"
    ))
    for maximum in (2e-4, None):
        convergence = read_newton_convergence(log, maximum)
        assert convergence.final_error == expected
        assert convergence.passed == (
            expected is not None and (maximum is None or expected <= maximum)
        )
        # Failed parsing/nonfinite values must still produce a valid JSON report.
        json.dumps(convergence.as_report(), allow_nan=False)


@pytest.mark.parametrize("mode", ["same_layout", "cross_layout", "protected_output"])
def test_final_output_selection_tolerance_profile_and_cli(tmp_path, mode):
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    (run / "outputs").mkdir()
    reference, checkpoint, final = (run / name for name in (
        "inputs/reference.h5", "outputs/result_0000.h5", "outputs/result.h5",
    ))
    write_solution(reference)
    write_solution(checkpoint, solution_offset=1.0)
    write_solution(final)
    (run / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": "baseline_warm",
        "layout_id": "serial_omp16" if mode == "cross_layout" else "mpi4_omp4",
    }))
    (run / "run_metadata.json").write_text(json.dumps({
        "status": "completed", "hdf5_outputs": ["outputs/result.h5", "outputs/result_0000.h5"],
    }))
    (run / "stdout.log").write_text(
        f"Error: 1.0E-5\nOutput written to file {checkpoint}\nOutput written to file {final}\n"
    )
    report_path = final if mode == "protected_output" else tmp_path / "comparison.json"
    original = final.read_bytes()
    if mode == "cross_layout":
        compare_completed_run(run, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json",
                              report_override=report_path)
    else:
        completed = run_command("compare", str(run), "--report", str(report_path))
        assert completed.returncode == int(mode == "protected_output"), completed.stderr
    assert final.read_bytes() == original
    if mode == "protected_output":
        assert "cannot replace a comparison input" in completed.stderr
        return
    report = json.loads(report_path.read_text())
    assert report["candidate"] == str(final)
    assert report["status"] == "passed"
    assert report["tolerance_profile"]["id"] == (
        "fixed_cross_layout" if mode == "cross_layout" else "fixed_same_layout"
    )
