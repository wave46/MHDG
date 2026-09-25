"""Settings ownership, build selection and read-only setup checks."""

import json

import pytest


from regression_tests.config import settings, build_settings
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.harness import create_harness, run_command


@pytest.fixture
def setup(tmp_path):
    harness = create_harness(tmp_path)
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


def test_doctor_checks_artifacts_without_creating_scratch(setup):
    harness, path, _, _ = setup
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path), catalog=harness.catalog)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "doctor: passed" in result.stdout
    assert not harness.run_root.exists()
    manifest = json.loads((harness.bundle / "manifest.json").read_text())
    restart = harness.bundle / manifest["artifacts"][manifest["roles"]["warm_restart"]]["path"]
    restart.unlink()  # Optional for the base bundle, required by the selected warm check.
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path), catalog=harness.catalog)
    assert result.returncode == 1 and "missing required artifact roles: warm_restart" in result.stdout
    assert not harness.run_root.exists()
    path.write_text("{}")
    result = run_command("doctor", "warm", "--case", "legacy_case", "--settings", str(path), catalog=harness.catalog)
    assert result.returncode == 1
    assert "no bundle selected" in result.stdout
    assert "no build selected" in result.stdout
