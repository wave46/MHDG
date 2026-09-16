"""Suite execution, trustworthy resume and saved scientific comparisons."""

import json
from pathlib import Path
from unittest.mock import Mock

import pytest
import h5py

from regression_tests import reporting, suites
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.harness import create_harness, run_command as package_command, REGRESSION_ROOT as ROOT

SOLVER = """#!/usr/bin/env bash
set -euo pipefail
cp inputs/reference.h5 outputs/result.h5
printf 'Error: 1.0E-5\\nOutput written to file outputs/result.h5\\n'
"""


def run_command(*args):
    return package_command(*args, catalog=ROOT)


@pytest.fixture
def harness(tmp_path, monkeypatch):
    result = create_harness(tmp_path, solver=SOLVER)
    monkeypatch.setitem(globals(), "ROOT", result.catalog)
    return result


def run(harness, name="warm", run_id="check", **kwargs):
    kwargs.setdefault("parameter_overrides", {"balance_diagnostics_mode": "off"})
    return suites.run_suite(harness.values, name, ROOT / "cases", ROOT / "layouts.json",
                            ROOT / "suites.json", ROOT / "tolerances.json", run_id, case_id="legacy_case", **kwargs)


def test_resume_reuses_good_outputs_preserves_interruption_and_retries_damage(harness, monkeypatch):
    # Same executable throughout; one serial launch fails once.
    harness.install_solver(SOLVER.replace('cp inputs', '''if [[ "$0" == */serial && "$OMP_NUM_THREADS" == 1 && ! -e "$(dirname "$0")/failed-once" ]]; then
  touch "$(dirname "$0")/failed-once"
  exit 7
fi
cp inputs'''))
    partial = harness.run_directory("warm", "mpi4_omp4", "check")
    execute = suites.run_cell
    def interrupt(inputs, workflow, layout):
        if layout == "mpi4_omp4":
            partial.mkdir(parents=True)
            (partial / "evidence.txt").write_text("keep")
            raise KeyboardInterrupt
        return execute(inputs, workflow, layout)
    with monkeypatch.context() as patch:
        patch.setattr(suites, "run_cell", interrupt)
        with pytest.raises(KeyboardInterrupt):
            run(harness, "warm_parallelism")
    first = json.loads((harness.run_root / "suites/warm_parallelism/legacy_case/check/suite_summary.json").read_text())
    failed = next(item for item in first["results"] if item["run_status"] == "solver_failed")
    good = {item["layout_id"]: item["run_directory"] for item in first["results"] if item["run_status"] == "completed"}
    assert good
    _, resumed = run(harness, "warm_parallelism", resume=True)
    assert resumed["status"] == "passed"
    assert len(resumed["results"]) == len(first["results"]) + 1
    assert (partial / "evidence.txt").read_text() == "keep"
    for item in resumed["results"]:
        if item["layout_id"] in good:
            assert item["run_directory"] == good[item["layout_id"]]
        else:
            assert Path(item["run_directory"]).name == "check-resume-1"
    assert (Path(failed["run_directory"]) / "run_metadata.json").is_file()
    damaged = next(iter(good))
    output = Path(good[damaged]) / "outputs/result.h5"
    output.write_bytes(b"changed output")
    _, repaired = run(harness, "warm_parallelism", resume=True)
    changed = next(item for item in repaired["results"] if item["layout_id"] == damaged)
    assert repaired["status"] == "passed"
    assert changed["run_directory"] != good[damaged]
    assert output.read_bytes() == b"changed output"  # Old evidence is not overwritten.

    # Reassessment must use current checks without executing valid outputs again.
    monkeypatch.setattr(suites, "run_cell", lambda *args: pytest.fail("valid output should be reused"))
    comparison = Mock(return_value=("fixed_hdf5", harness.root / "comparison.json", {
        "status": "passed", "failures": [], "convergence": {"passed": False}}))
    monkeypatch.setattr(suites, "compare_completed_run", comparison)
    path, rejected = run(harness, "warm_parallelism", resume=True)
    assert rejected["status"] == "failed" and comparison.call_count == 4
    assert [row["run_directory"] for row in rejected["results"]] == [row["run_directory"] for row in repaired["results"]]
    verified = suites.verify_suite(path, ROOT / "cases", ROOT / "tolerances.json")[1]
    assert verified["results"] == rejected["results"]
    shared = ROOT / "workflows.json"
    original = shared.read_bytes()
    shared.write_bytes(original + b"\n")
    with pytest.raises(BundleError, match="workflow_catalog"):
        run(harness, "warm_parallelism", resume=True)
    shared.write_bytes(original)
    runtime = harness.runtime_file.read_bytes()
    harness.runtime_file.write_bytes(b"changed runtime input")
    with pytest.raises(BundleError, match="checksum/size changed"):
        run(harness, "warm_parallelism", resume=True)
    harness.runtime_file.write_bytes(runtime)
    harness.install_solver("#!/bin/sh\nexit 7\n", "parallel")
    with pytest.raises(BundleError, match="executables"):
        run(harness, "warm_parallelism", resume=True)


def test_cli_requires_golden_and_records_diagnostic_overrides(harness):
    result = run_command("check", "warm", "--case", "legacy_case", "--settings", str(harness.settings))
    assert result.returncode == 1 and "requires bundle_class=golden" in result.stderr
    harness.set_bundle_class("golden")
    args = ("check", "warm", "--case", "legacy_case", "--settings", str(harness.settings),
            "--build-manifest", str(harness.build_manifest), "--run-id", "mode", "--run-only")
    result = run_command(*args, "--diagnostics", "detailed")
    assert result.returncode == 0, result.stderr
    assert "suite deferred:" in result.stdout
    assert "balance diagnostics deferred:" in result.stdout
    directory = harness.run_directory("warm", "mpi4_omp4", "mode")
    assert "balance_diagnostics_mode = 'detailed'" in (directory / "param.txt").read_text()
    # Explicit selection survives unrelated/default settings changes on resume.
    machine = json.loads(harness.settings.read_text())
    machine["defaults"]["build"] = "unavailable-build.json"
    machine["defaults"]["bundles"]["unselected_case"] = "unavailable-bundle"
    machine["build_jobs"] = 3
    harness.settings.write_text(json.dumps(machine))
    result = run_command(*args, "--diagnostics", "detailed", "--resume")
    assert result.returncode == 0 and "reusing completed" in result.stdout, result.stderr
    result = run_command(*args, "--diagnostics", "off", "--resume")
    assert result.returncode == 1 and "parameter_overrides" in result.stderr
    result = run_command(*args, "--build", "--resume")
    assert result.returncode == 1 and "--build cannot be used with --resume" in result.stderr


def test_parallel_comparison_and_saved_recheck(harness):
    solver = SOLVER.replace('cp inputs/reference.h5', 'cp "$(dirname "$0")/seed.h5"').replace('1.0E-5', '1.0E5')
    harness.install_solver(solver.replace("cp ", "printf 'mesh\\n' > res/temp.msh\ncp ", 1))
    manifest_path = harness.bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["roles"]["warm_reference"]  # Parallel checks have their own run as reference.
    manifest_path.write_text(json.dumps(manifest))
    path, summary = run(harness, "parallel")
    assert summary["status"] == "passed"
    assert len(summary["comparisons"]) == 1
    report = json.loads(Path(summary["comparisons"][0]["comparison_report"]).read_text())
    assert report["convergence"] == {"passed": True, "final_newton_error": 1e5, "maximum": None}
    result = run_command("compare", "--suite", str(path))
    assert result.returncode == 0, result.stderr
    assert "diagnostic comparison skipped: diagnostics off" in result.stdout
    assert "mode off: output absence checked" in result.stdout
    mesh = next(Path(summary["comparisons"][0]["candidate_run_directory"]).glob("**/res/temp.msh"))
    mesh.write_text("mesh\n\n")
    _, verification = suites.verify_suite(path, ROOT / "cases", ROOT / "tolerances.json")
    assert verification["status"] == "failed"
    assert any("generated mesh differs" in failure
               for pair in verification["comparisons"] for failure in pair["failures"])


def test_offline_checks_references_pairs_and_diagnostics(tmp_path, monkeypatch, capsys):
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
    # Individual failures must stay visible when the report also contains pairs.
    passed.update(status="failed", failures=["individual field failure"])
    pairs.return_value = [{"workflow_id": "cold_fixed", "baseline_layout_id": "serial_omp1",
                           "candidate_layout_id": "serial_omp16", "comparison_policy": "fixed_hdf5",
                           "status": "failed", "failures": ["parallel field failure"]}]
    report = check()
    reporting.print_verification_summary(report, source)
    output = capsys.readouterr().out
    assert all(message in output for message in ("individual field failure", "parallel field failure", "missing diagnostics"))
    pairs.reset_mock()
    suites.verify_suite(source, ROOT / "cases", ROOT / "tolerances.json", include_layout_pairs=False)
    pairs.assert_not_called()


@pytest.mark.parametrize("suite,damage", [
    ("initialization", "nonfinite"),
    ("bootstrap", "nonconvergence"),
])
def test_reference_free_runs_validate_outputs_and_stage_convergence(harness, suite, damage, monkeypatch):
    seed = harness.serial_executable.parent / "seed.h5"
    solver = r'''#!/usr/bin/env bash
set -euo pipefail
cp "$(dirname "$0")/seed.h5" outputs/result.h5
error=1.0E-5
[[ "$PWD" != *01_initial ]] || error=1.0
printf 'Error: %s\nOutput written to file outputs/result.h5\n' "$error"
'''
    harness.install_solver(solver)
    path, summary = run(harness, suite)
    assert summary["status"] == "passed"
    first = summary["results"][0]
    assert first["comparison_status"] == "not_run"
    stages = first["validation"]["stages"]
    assert stages[0]["convergence"] == {"passed": True, "final_newton_error": 1., "maximum": None}
    assert suites.verify_suite(path, ROOT / "cases", ROOT / "tolerances.json")[1]["status"] == "passed"
    if damage == "nonconvergence":
        # An intermediate cold stage must converge even if the last stage does.
        log = Path(stages[1]["output"]).parent.parent / "stdout.log"
        log.write_text("Error: 1e-3\n")
        expected = "continued: final Newton error exceeds tolerance"
        monkeypatch.setattr(suites, "run_cell", lambda *args: pytest.fail("reuse completed outputs"))
        assert run(harness, suite, resume=True)[1]["status"] == "failed"
    else:
        def corrupt(output):
            with h5py.File(output, "r+") as handle:
                handle["solution/u"][0] = float("nan")
        corrupt(Path(stages[0]["output"]))
        corrupt(seed)
        expected = "solution/u"
        assert run(harness, suite, run_id="bad-output")[1]["status"] == "failed"
    _, verification = suites.verify_suite(path, ROOT / "cases", ROOT / "tolerances.json")
    assert verification["status"] == "failed"
    assert any(expected in failure for result in verification["results"] for failure in result["failures"])
    if damage == "nonconvergence":
        completed = run_command("compare", "--suite", str(path))
        assert completed.returncode == 1 and expected in completed.stdout


def test_profile_preflights_cases_and_forwards_suite_resume(harness, monkeypatch):
    from regression_tests.bundles import create_bundle

    diverted = harness.root / "diverted"
    create_bundle("diverted_case", harness.source, diverted, ROOT / "cases")
    values = harness.values
    settings = {"legacy_case": values,
                "diverted_case": {**values, "MHDG_REGRESSION_DATA_ROOT": str(diverted)}}
    checks = [{"suite_id": "warm", **suites.load_suite_definition(
        "warm", ROOT / "suites.json", ROOT / "layouts.json", ROOT / "cases", case_id=case,
    )} for case in settings]
    calls = []
    def selected_suite(values, suite_id, *args, case_id, **kwargs):
        resumed = args[7]
        calls.append((case_id, resumed))
        path = harness.run_root / f"suites/{suite_id}/{case_id}/profile/suite_summary.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        report = {"status": "failed" if case_id == "diverted_case" and not resumed else "passed",
                  "results": [], "duration_seconds": 0.}
        path.write_text(json.dumps(report))
        return path, report
    monkeypatch.setattr(suites, "run_suite", selected_suite)

    manifest = diverted / "manifest.json"
    original = manifest.read_text()
    missing = json.loads(original)
    del missing["roles"]["warm_reference"]
    manifest.write_text(json.dumps(missing))
    with monkeypatch.context() as patch:
        patch.setattr(suites, "run_suite", lambda *args, **kwargs: pytest.fail("must preflight later suites first"))
        with pytest.raises(BundleError, match="missing required artifact roles.*warm_reference"):
            suites.run_profile("example", checks, settings, ROOT, "profile")
    assert not harness.run_root.exists()
    manifest.write_text(original)
    path, first = suites.run_profile("example", checks, settings, ROOT, "profile")
    assert first["status"] == "failed"
    assert [item["status"] for item in first["results"]] == ["passed", "failed"]
    assert calls == [("legacy_case", False), ("diverted_case", False)]
    calls.clear()
    resumed_path, resumed = suites.run_profile("example", checks, settings, ROOT, "profile", resume=True)
    assert resumed_path == path
    assert resumed["status"] == "passed"
    assert calls == [("legacy_case", True), ("diverted_case", True)]
    assert len({item["summary"] for item in resumed["results"]}) == 2
