"""Settings ownership, build selection and read-only setup checks."""

import json

import pytest


from regression_tests.config import settings, build_settings
from regression_tests.support import BundleError, DocumentError
from regression_tests.tests.fixtures.harness import create_harness, run_command


@pytest.fixture
def setup(tmp_path):
    harness = create_harness(tmp_path, solver=(
        '#!/bin/sh\ncp inputs/reference.h5 outputs/result.h5\n'
        'echo "Error: 1.0E-5"\necho "Output written to file outputs/result.h5"\n'
    ))

    manifest, path = harness.build_manifest, harness.settings
    document = json.loads(path.read_text())
    return harness, path, manifest, document


def test_explicit_selections_override_case_defaults_and_builds(setup, monkeypatch):
    harness, path, manifest, document = setup
    monkeypatch.setenv("MHDG_REGRESSION_SETTINGS", str(path))
    values = settings(case="legacy_case")
    assert values["MHDG_REGRESSION_DATA_ROOT"] == str(harness.bundle)
    assert values["MHDG_EXECUTABLES"]["NGammaTiTeNeutral/serial"] == str(harness.serial_executable)
    assert values["MHDG_REGRESSION_RUN_ROOT"] == str(harness.run_root)
    document["defaults"]["build"] = "unavailable-build.json"
    document["defaults"]["bundles"]["another_case"] = "separate-data"
    path.write_text(json.dumps(document))
    values = settings(path, case="another_case", build_manifest=manifest)
    assert values["MHDG_REGRESSION_DATA_ROOT"] == str(path.parent / "separate-data")
    values = settings(path, case="another_case", bundle=harness.bundle, build_manifest=manifest)
    assert values["MHDG_REGRESSION_DATA_ROOT"] == str(harness.bundle)
    assert "MHDG_EXECUTABLES" not in settings(path, use_build=False)


def test_configuration_rejects_user_owned_generated_fields(setup):
    _, path, _, document = setup
    document["solver_revision"] = "user-entered"
    path.write_text(json.dumps(document))
    with pytest.raises(BundleError, match="unknown machine settings fields: solver_revision"):
        settings(path)
    path.write_text("MHDG_REGRESSION_DATA_ROOT=/legacy/data\n")
    with pytest.raises(DocumentError, match="invalid JSON in machine settings"):
        settings(path)


def test_build_selection_rejects_failed_or_changed_artifacts(setup):
    harness, _, manifest, _ = setup
    document = json.loads(manifest.read_text())
    manifest.write_text(json.dumps({**document, "status": "failed"}))
    with pytest.raises(BundleError, match="completed"):
        build_settings(manifest)
    manifest.write_text(json.dumps(document))
    harness.runtime_file.write_text("changed runtime input")
    with pytest.raises(BundleError, match="checksum/size changed"):
        build_settings(manifest)


def test_json_check_resumes_with_same_effective_selection(setup):
    harness, path, manifest, document = setup
    arguments = ("check", "warm", "--case", "legacy_case", "--allow-candidate", "--run-only", "--run-id", "json-settings",
                 "--settings", str(path), "--build-manifest", str(manifest))
    result = run_command(*arguments)
    assert result.returncode == 0, result.stderr
    document["defaults"]["bundles"]["unselected_case"] = "unavailable-bundle"
    document["defaults"]["build"] = "unavailable-build.json"
    document["build_jobs"] = 3
    path.write_text(json.dumps(document))
    resumed = run_command(*arguments, "--resume")
    assert resumed.returncode == 0, resumed.stderr
    assert "reusing completed" in resumed.stdout
    summary = json.loads((harness.run_root / "suites/warm/legacy_case/json-settings/suite_summary.json").read_text())
    assert summary["execution_inputs"]["build_manifest"]["path"] == str(manifest)


def test_doctor_checks_artifacts_without_creating_scratch(setup):
    harness, path, _, _ = setup
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path))
    assert result.returncode == 0, result.stdout + result.stderr
    assert "doctor: passed" in result.stdout
    assert not harness.run_root.exists()
    manifest = json.loads((harness.bundle / "manifest.json").read_text())
    restart = harness.bundle / manifest["artifacts"][manifest["roles"]["warm_restart"]]["path"]
    original = restart.read_bytes()
    restart.unlink()  # Optional for the base bundle, required by the selected warm check.
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path))
    assert result.returncode == 1 and "missing required artifact roles: warm_restart" in result.stdout
    assert not harness.run_root.exists()
    restart.write_bytes(original)
    artifact = next(iter(manifest["artifacts"].values()))
    (harness.bundle / artifact["path"]).write_text("changed bundle input")
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path))
    assert result.returncode == 1
    assert "FAIL legacy_case bundle:" in result.stdout
    assert not harness.run_root.exists()
    path.write_text("{}")
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path))
    assert result.returncode == 1
    assert "no bundle selected" in result.stdout
    assert "no build selected" in result.stdout
