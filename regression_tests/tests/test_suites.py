"""Suite execution, trustworthy resume and saved scientific comparisons."""

import json
import shutil
from pathlib import Path
from unittest.mock import Mock

import pytest
import h5py

from regression_tests import suites
from regression_tests.compare import compare_generated_meshes, producer_converged
from support.errors import BundleError
from bundle.settings import read_settings
from tests.fixtures.harness import create_harness, run_command, REGRESSION_ROOT as ROOT
from tests.fixtures.solutions import write_solution

SOLVER = """#!/usr/bin/env bash
set -euo pipefail
cp inputs/reference.h5 outputs/result.h5
printf 'Error: 1.0E-5\\nOutput written to file outputs/result.h5\\n'
"""


def write_off_solution(path):
    write_solution(path)
    with h5py.File(path, "r+") as handle:
        handle["simulation_parameters/switches/balance_diagnostics_mode"] = b"off"


@pytest.fixture
def harness(tmp_path):
    return create_harness(tmp_path, solver=SOLVER, reference_writer=write_off_solution)


def run(harness, name="warm", run_id="check", **kwargs):
    kwargs.setdefault("parameter_overrides", {"balance_diagnostics_mode": "off"})
    return suites.run_suite(harness.settings, name, ROOT / "cases", ROOT / "layouts.json",
                            ROOT / "suites.json", ROOT / "tolerances.json", run_id, case_id="legacy_case", **kwargs)


def test_resume_reuses_good_outputs_retries_failed_and_changed_runs(harness):
    # Same executable throughout; one serial launch fails once.
    harness.install_solver(SOLVER.replace('cp inputs', '''if [[ "$0" == */serial && "$OMP_NUM_THREADS" == 1 && ! -e "$(dirname "$0")/failed-once" ]]; then
  touch "$(dirname "$0")/failed-once"
  exit 7
fi
cp inputs'''))
    _, first = run(harness, "warm_parallelism", compare=False)
    assert first["status"] == "failed"
    failed = next(item for item in first["results"] if item["run_status"] == "solver_failed")
    good = {item["layout_id"]: item["run_directory"] for item in first["results"] if item["run_status"] == "completed"}
    assert good
    _, resumed = run(harness, "warm_parallelism", compare=False, resume=True)
    assert resumed["status"] == "passed"
    assert len(resumed["results"]) == len(first["results"])
    for item in resumed["results"]:
        if item["layout_id"] in good:
            assert item["run_directory"] == good[item["layout_id"]]
        else:
            assert Path(item["run_directory"]).name == "check-resume-1"
    assert (Path(failed["run_directory"]) / "run_metadata.json").is_file()
    damaged = next(iter(good))
    output = Path(good[damaged]) / "outputs/result.h5"
    output.write_bytes(b"changed output")
    _, repaired = run(harness, "warm_parallelism", compare=False, resume=True)
    changed = next(item for item in repaired["results"] if item["layout_id"] == damaged)
    assert repaired["status"] == "passed"
    assert changed["run_directory"] != good[damaged]
    assert output.read_bytes() == b"changed output"  # Old evidence is not overwritten.


def test_interruption_preserves_partial_directory(harness, monkeypatch):
    partial = harness.run_directory("warm", "mpi4_omp4", "check")
    def interrupt(*args):
        partial.mkdir(parents=True)
        (partial / "evidence.txt").write_text("keep")
        raise KeyboardInterrupt
    with monkeypatch.context() as patch:
        patch.setattr(suites, "run_cell", interrupt)
        with pytest.raises(KeyboardInterrupt):
            run(harness, compare=False)
    _, summary = run(harness, compare=False, resume=True)
    assert summary["status"] == "passed"
    assert (partial / "evidence.txt").read_text() == "keep"
    assert Path(summary["results"][0]["run_directory"]).name == "check-resume-1"


def test_resume_rechecks_comparisons_and_rejects_changed_inputs(harness, monkeypatch, tmp_path):
    catalog = tmp_path / "catalog"
    shutil.copytree(ROOT / "cases", catalog / "cases")
    shutil.copytree(ROOT / "schemas", catalog / "schemas")
    shutil.copy2(ROOT / "workflows.json", catalog / "workflows.json")
    selected = read_settings(harness.settings)
    def check(resume=False):
        return suites.run_suite(selected, "warm", catalog / "cases", ROOT / "layouts.json",
                                 ROOT / "suites.json", ROOT / "tolerances.json", "checked", resume=resume, case_id="legacy_case")
    _, first = check()
    assert first["status"] == "passed"
    # An unused executable and build-only preferences do not affect these runs.
    selected["MHDG_SERIAL_EXECUTABLE"] = "/unused/serial"
    selected["MHDG_REGRESSION_BUILD_JOBS"] = "3"
    monkeypatch.setattr(suites, "run_cell", lambda *args: pytest.fail("valid output should be reused"))
    comparison = Mock(return_value=("fixed_hdf5", tmp_path / "comparison.json", {
        "status": "failed", "failures": ["changed acceptance result"],
    }))
    monkeypatch.setattr(suites, "compare_completed_run", comparison)
    _, resumed = check(resume=True)
    comparison.assert_called_once()
    assert resumed["status"] == "failed"
    assert resumed["results"][0]["run_directory"] == first["results"][0]["run_directory"]
    shared = catalog / "workflows.json"
    original = shared.read_text()
    shared.write_text(original + "\n")
    with pytest.raises(BundleError, match="workflow_catalog"):
        check(resume=True)
    shared.write_text(original)
    runtime = harness.runtime_file.read_bytes()
    harness.runtime_file.write_bytes(b"changed runtime input")
    with pytest.raises(BundleError, match="runtime_files"):
        check(resume=True)
    harness.runtime_file.write_bytes(runtime)
    harness.install_solver("#!/usr/bin/env bash\nexit 7\n", "parallel")
    with pytest.raises(BundleError, match="parallel_executable"):
        check(resume=True)


def test_cli_requires_golden_and_records_diagnostic_overrides(harness):
    result = run_command("check", "warm", "--case", "legacy_case", "--settings", str(harness.settings))
    assert result.returncode == 1 and "requires bundle_class=golden" in result.stderr
    harness.set_bundle_class("golden")
    args = ("check", "warm", "--case", "legacy_case", "--settings", str(harness.settings), "--run-id", "mode", "--run-only")
    result = run_command(*args, "--diagnostics", "detailed")
    assert result.returncode == 0, result.stderr
    directory = harness.run_directory("warm", "mpi4_omp4", "mode")
    assert "balance_diagnostics_mode = 'detailed'" in (directory / "param.txt").read_text()
    result = run_command(*args, "--diagnostics", "off", "--resume")
    assert result.returncode == 1 and "parameter_overrides" in result.stderr
    result = run_command(*args, "--build", "--resume")
    assert result.returncode == 1 and "--build cannot be used with --resume" in result.stderr


def test_parallel_comparison_and_saved_recheck(harness):
    write_off_solution(harness.serial_executable.parent / "race_result.h5")
    harness.install_solver(SOLVER.replace('cp inputs/reference.h5', 'cp "$(dirname "$0")/race_result.h5"').replace('1.0E-5', '1.0E5'))
    path, summary = run(harness, "parallel")
    assert summary["status"] == "passed"
    assert summary["comparisons"]
    report = json.loads(Path(summary["comparisons"][0]["comparison_report"]).read_text())
    assert report["convergence"] == {"passed": True, "final_newton_error": 1e5, "maximum": None}
    result = run_command("compare", "--suite", str(path))
    assert result.returncode == 0, result.stderr


def test_generated_mesh_comparison_is_byte_exact(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    for root in (left, right):
        (root / "res").mkdir(parents=True)
        (root / "res/temp.msh").write_text("mesh\n")
    assert compare_generated_meshes(left, right)["passed"]
    (right / "res/temp.msh").write_text("mesh\n\n")
    assert not compare_generated_meshes(left, right)["passed"]


def test_offline_checks_references_pairs_and_diagnostics(tmp_path, monkeypatch):
    passed = {"status": "passed", "failures": [], "convergence": {"passed": True}}
    compare = Mock(return_value=("fixed_hdf5", tmp_path / "compare.json", passed))
    pairs = Mock(return_value=[{"status": "passed"}])
    diagnostics = Mock(return_value={"status": "passed", "outputs": [], "failures": []})
    monkeypatch.setattr(suites, "compare_completed_run", compare)
    monkeypatch.setattr(suites, "compare_layout_pairs", pairs)
    monkeypatch.setattr("regression_tests.diagnostics.check_suite", diagnostics)
    results = [dict(workflow_id="cold_fixed", layout_id="serial_omp1", run_status="completed", run_directory=str(tmp_path))]
    source = tmp_path / "suite_summary.json"
    data = {"schema_version": 2, "suite_id": "example", "run_id": "example", "case_id": "legacy_case",
            "workflow_ids": ["cold_fixed"], "reference_comparisons": True, "results": results,
            "layout_comparisons": [{"baseline": "serial_omp1", "candidate": "serial_omp16"}], "tolerance_profile": "cold_cross_layout"}
    def check():
        source.write_text(json.dumps(data))
        return suites.verify_suite(source, ROOT / "cases", ROOT / "tolerances.json")[1]
    assert check()["status"] == "passed"
    compare.assert_called_once()
    pairs.assert_called_once()
    diagnostics.return_value = {"status": "failed", "outputs": [], "failures": ["missing diagnostics"]}
    assert check()["status"] == "failed"
    results.append({**results[0], "run_status": "solver_failed"})
    compare.reset_mock()
    report = check()
    compare.assert_called_once()  # Incomplete runs do not reach the field comparator.
    assert report["results"][-1]["convergence_status"] is None
    assert report["status"] == "failed"
    pairs.reset_mock()
    suites.verify_suite(source, ROOT / "cases", ROOT / "tolerances.json", include_layout_pairs=False)
    pairs.assert_not_called()


def test_staged_convergence_is_independent_of_old_reference_failure(tmp_path):
    stages = [{"stage_id": "initial", "newton_check": "finite_only"}, {"stage_id": "settle", "newton_check": "bounded"}]
    records = []
    for stage, error in zip(stages, (1., 1e-5)):
        directory = tmp_path / stage["stage_id"]
        directory.mkdir()
        (directory / "stdout.log").write_text(f"Error: {error}\n")
        records.append({**stage, "status": "completed", "run_directory": str(directory)})
    (tmp_path / "run_metadata.json").write_text(json.dumps({"stages": records}))
    source = {"run_directory": str(tmp_path), "layout_id": "serial_omp1"}
    workflow = {"stages": stages, "comparison_policy": "fixed_hdf5", "stage_tolerance_profile": "fixed_stage_reference"}
    report = {"status": "failed", "failures": ["old fields differ"]}
    def converged():
        return producer_converged(source, workflow, "reference_matrix", report, ROOT / "tolerances.json")
    assert converged()
    (tmp_path / "settle/stdout.log").write_text("Error: 1e-3\n")
    assert not converged()


def test_profile_resumes_across_cases_without_repeating_success(harness, monkeypatch):
    from bundle.creation import create_bundle

    diverted = harness.root / "diverted"
    create_bundle("diverted_case", harness.source, diverted, ROOT / "cases")
    values = read_settings(harness.settings)
    settings = {"legacy_case": values,
                "diverted_case": {**values, "MHDG_REGRESSION_DATA_ROOT": str(diverted)}}
    checks = [{"suite_id": "warm", "case_id": case} for case in settings]
    harness.install_solver(SOLVER.replace('cp inputs', '''if [[ "$PWD" == */diverted_case/* && ! -e "$(dirname "$0")/failed-once" ]]; then
  touch "$(dirname "$0")/failed-once"
  exit 7
fi
cp inputs'''))
    path, first = suites.run_profile("example", checks, settings, ROOT, "profile")
    assert first["status"] == "failed"
    assert [item["status"] for item in first["results"]] == ["passed", "failed"]
    calls = []
    run_cell = suites.run_cell
    def record(inputs, *args):
        calls.append(inputs.case_id)
        return run_cell(inputs, *args)
    monkeypatch.setattr(suites, "run_cell", record)
    resumed_path, resumed = suites.run_profile("example", checks, settings, ROOT, "profile", resume=True)
    assert resumed_path == path
    assert resumed["status"] == "passed"
    assert calls == ["diverted_case"]
    assert len({item["summary"] for item in resumed["results"]}) == 2
