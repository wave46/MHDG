"""Settings ownership, build selection and read-only setup checks."""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / "tools"))

from regression_tests.config import settings, build_settings  # noqa: E402
from support.errors import BundleError  # noqa: E402
from support.files import file_identity  # noqa: E402
from tests.fixtures.harness import create_harness, run_command  # noqa: E402


@pytest.fixture
def setup(tmp_path):
    harness = create_harness(tmp_path, solver=(
        '#!/bin/sh\ncp inputs/reference.h5 outputs/result.h5\n'
        'echo "Error: 1.0E-5"\necho "Output written to file outputs/result.h5"\n'
    ))

    def record(path):
        return {"path": str(path.relative_to(tmp_path)), **file_identity(path)}

    manifest = tmp_path / "build_metadata.json"
    manifest.write_text(json.dumps({
        "schema_version": 2, "status": "completed", "build_id": "test-build",
        "repository": {"revision": "test-revision"},
        "profile": {"model": "NGammaTiTeNeutral", "dimension": "2D"},
        "artifacts": {
            "serial": record(harness.serial_executable),
            "parallel": record(harness.parallel_executable),
        },
        "runtime_files": {harness.runtime_file.name: record(harness.runtime_file)},
    }))
    document = {
        "run_root": "runs", "mpi_launcher": str(harness.mpi_launcher),
        "defaults": {"build": manifest.name, "bundles": {"legacy_case": "bundle"}},
    }
    path = tmp_path / "machine.json"
    path.write_text(json.dumps(document))
    return harness, path, manifest, document


def test_explicit_selections_override_case_defaults_and_builds(setup, monkeypatch):
    harness, path, manifest, document = setup
    monkeypatch.setenv("MHDG_REGRESSION_SETTINGS", str(path))
    values = settings(case="legacy_case")
    assert values["MHDG_REGRESSION_DATA_ROOT"] == str(harness.bundle)
    assert values["MHDG_SERIAL_EXECUTABLE"] == str(harness.serial_executable)
    assert values["MHDG_REGRESSION_RUN_ROOT"] == str(harness.run_root)
    document["defaults"]["build"] = "unavailable-build.json"
    document["defaults"]["bundles"]["another_case"] = "separate-data"
    path.write_text(json.dumps(document))
    values = settings(path, case="another_case", build_manifest=manifest)
    assert values["MHDG_REGRESSION_DATA_ROOT"] == str(path.parent / "separate-data")
    values = settings(path, case="another_case", bundle=harness.bundle, build_manifest=manifest)
    assert values["MHDG_REGRESSION_DATA_ROOT"] == str(harness.bundle)
    assert "MHDG_SERIAL_EXECUTABLE" not in settings(path, use_build=False)


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


def test_json_check_resumes_with_same_effective_selection(setup):
    harness, path, manifest, document = setup
    arguments = ("check", "--allow-candidate", "--run-only", "--run-id", "json-settings",
                 "--settings", str(path), "--build-manifest", str(manifest))
    result = run_command(*arguments)
    assert result.returncode == 0, result.stderr
    document["defaults"]["bundles"]["unselected_case"] = "unavailable-bundle"
    document["defaults"]["build"] = "unavailable-build.json"
    document["build_jobs"] = 3
    path.write_text(json.dumps(document))
    resumed = run_command(*arguments, "--resume")
    assert resumed.returncode == 0, resumed.stderr
    assert "skipping recorded" in resumed.stdout
    summary = json.loads((harness.run_root / "suites/warm/json-settings/suite_summary.json").read_text())
    assert summary["execution_inputs"]["settings"]["MHDG_BUILD_MANIFEST"] == str(manifest)


def test_doctor_checks_artifacts_without_creating_scratch(setup):
    harness, path, _, _ = setup
    result = run_command("doctor", "--settings", str(path))
    assert result.returncode == 0, result.stdout + result.stderr
    assert "doctor: passed" in result.stdout
    assert not harness.run_root.exists()
    manifest = json.loads((harness.bundle / "manifest.json").read_text())
    artifact = next(iter(manifest["artifacts"].values()))
    (harness.bundle / artifact["path"]).write_text("changed bundle input")
    result = run_command("doctor", "--settings", str(path))
    assert result.returncode == 1
    assert "FAIL bundle:" in result.stdout
    assert not harness.run_root.exists()
    path.write_text("{}")
    result = run_command("doctor", "--settings", str(path))
    assert result.returncode == 1
    assert "no bundle selected" in result.stdout
    assert "no build selected" in result.stdout
