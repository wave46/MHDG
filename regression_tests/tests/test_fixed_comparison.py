"""Direct field/mesh contracts and the public comparison command."""

import json

import h5py
import numpy as np
import pytest

from regression_tests.compare_fixed import compare_hdf5_files
from regression_tests.tests.fixtures.harness import run_command
from regression_tests.tests.fixtures.solutions import write_solution

TOLERANCES = {
    "newton_error_max": 2e-4, "mesh_coordinate_atol": 1e-12,
    "relative_l2_max": 1e-10, "normalized_linf_max": 1e-9,
}


@pytest.fixture
def files(tmp_path):
    reference, candidate = tmp_path / "reference.h5", tmp_path / "candidate.h5"
    write_solution(reference, grouped=True)
    write_solution(candidate, grouped=False)
    return reference, candidate


def test_grouped_and_flat_fields_match_with_reference_equation_names(files):
    reference, candidate = files
    with h5py.File(candidate, "r+") as handle:
        del handle["conservative_variable_names"]
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == "passed"
    assert report["formats"]["candidate"] == "flat_solution+flat_mesh"
    assert report["solution"]["equation_names"] == ["rho", "Gamma"]
    assert report["solution"]["datasets"]["u"]["equations"]["Gamma"]["relative_l2"] == 0


@pytest.mark.parametrize("dataset,index,value,section", [
    ("u", 1, 3.0, "solution"),
    ("q", 0, np.nan, "solution"),
    ("T", (0, 0), 99, "mesh"),
    ("transport_1d/coefficients/d_fs", 0, 2.0, "transport_1d"),
])
def test_field_nonfinite_mesh_and_transport_failures(files, dataset, index, value, section):
    reference, candidate = files
    with h5py.File(candidate, "r+") as handle:
        handle[dataset][index] = value
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == "failed"
    assert not report[section]["passed"]
    if dataset == "u":
        equations = report["solution"]["datasets"]["u"]["equations"]
        assert equations["rho"]["passed"] and not equations["Gamma"]["passed"]
    elif dataset == "q":
        assert not report["solution"]["datasets"]["q"]["equations"]["rho"]["finite"]
    elif dataset == "T":
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
    if difference == "missing":
        assert "candidate is missing magnetic data" in report["failures"]
    elif difference == "region":
        assert "magnetic/topology_region differs" in report["failures"]
    else:
        assert report["magnetic"]["datasets"]["topology"]["mode"] == "exact"
        assert report["magnetic"]["datasets"]["rho_pol_norm"]["relative_l2"] == 0


@pytest.mark.parametrize("mode", ["same_layout", "cross_layout", "nonconverged", "protected_output"])
def test_compare_command_selects_final_output_and_checks_convergence(tmp_path, mode):
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    (run / "outputs").mkdir()
    reference, checkpoint, final = (run / name for name in (
        "inputs/reference.h5", "outputs/result_0000.h5", "outputs/result.h5",
    ))
    write_solution(reference)
    write_solution(checkpoint, grouped=False, solution_offset=1.0)
    write_solution(final, grouped=False)
    (run / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": "warm",
        "layout_id": "serial_omp16" if mode == "cross_layout" else "mpi4_omp4",
    }))
    (run / "run_metadata.json").write_text(json.dumps({
        "status": "completed", "hdf5_outputs": ["outputs/result.h5", "outputs/result_0000.h5"],
    }))
    error = "3.0E-4" if mode == "nonconverged" else "1.0E-5"
    (run / "stdout.log").write_text(
        f"Error: 8.0E-4\nError: {error}\nOutput written to file {checkpoint}\nOutput written to file {final}\n"
    )
    report_path = final if mode == "protected_output" else tmp_path / "comparison.json"
    original = final.read_bytes()
    completed = run_command("compare", str(run), "--report", str(report_path))
    assert completed.returncode == int(mode in {"nonconverged", "protected_output"}), completed.stderr
    assert final.read_bytes() == original
    if mode == "protected_output":
        assert "cannot replace a comparison input" in completed.stderr
        return
    report = json.loads(report_path.read_text())
    assert report["candidate"] == str(final)
    assert report["convergence"]["final_newton_error"] == float(error)
    assert report["tolerance_profile"]["id"] == (
        "fixed_cross_layout" if mode == "cross_layout" else "fixed_same_layout"
    )
    assert len(report["files"]["reference"]["sha256"]) == 64
    assert report["convergence"]["passed"] == (mode != "nonconverged")
    for text in ("Mesh: PASS", "solution/u: PASS", "transport_1d: PASS"):
        assert text in completed.stdout
    assert ("Newton error: FAIL" if mode == "nonconverged" else "Newton error: PASS") in completed.stdout
