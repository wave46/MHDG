"""Bundle readiness, standalone copies and validation of selected inputs."""

import json
import pytest

from regression_tests.bundles import bundle_readiness, create_bundle, validate_bundle_root
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.case_data import write_case_source
from regression_tests.tests.fixtures.harness import run_command


@pytest.fixture
def data(tmp_path, request):
    return getattr(request, "param", "legacy_case"), write_case_source(tmp_path / "source"), tmp_path / "bundle"


@pytest.mark.parametrize("data", ["legacy_case", "diverted_case"], indirect=True)
def test_readiness_and_selected_creation_require_actual_files(data):
    case, source, output = data
    (source / "restart.h5").unlink()
    base = bundle_readiness(case, source)
    selected = bundle_readiness(case, source, workflows=["warm"])
    rows = {row["role"]: row for row in selected["artifacts"]}
    assert base["status"] == "ready"
    assert selected["missing"] == ["warm_restart"]
    assert rows["warm_restart"]["origin"] == "workflow-producible"
    assert rows["warm_restart"]["presence"] == "missing"
    assert rows["warm_restart"]["producers"]
    assert rows["geometry"]["origin"] == "user-supplied"
    command = run_command("bundle", "readiness", case, "--source", str(source), "--workflow", "warm")
    assert command.returncode == 1 and "workflow-producible" in command.stdout
    with pytest.raises(BundleError, match="required files missing.*warm_restart"):
        create_bundle(case, source, output, workflows=["warm"])
    assert not output.exists()
    create_bundle(case, source, output)
    validate_bundle_root(output)
    with pytest.raises(BundleError, match="missing required artifact roles.*warm_restart"):
        validate_bundle_root(output, workflows=["warm"])
    with pytest.raises(BundleError, match="no workflow"):
        bundle_readiness(case, source, workflows=["unknown"])


@pytest.mark.parametrize("data", ["legacy_case", "diverted_case"], indirect=True)
def test_creation_copies_links_and_validation_uses_manifest_paths(data, tmp_path):
    case, source, output = data
    shared = tmp_path / "equilibrium.h5"
    (source / "equilibrium.h5").rename(shared)
    (source / "equilibrium.h5").symlink_to(shared)
    result = run_command("bundle", "create", "--case", case, "--source", str(source),
                         "--output", str(output), "--workflow", "cold_step_adaptive")
    assert result.returncode == 0, result.stderr
    copied = output / "inputs/equilibrium.h5"
    assert copied.read_bytes() == shared.read_bytes() and not copied.is_symlink()
    shared.unlink()
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    artifact = manifest["artifacts"][manifest["roles"]["equilibrium_magnetic_field"]]
    copied.rename(output / "different.h5")
    artifact["path"] = "different.h5"
    manifest_path.write_text(json.dumps(manifest))
    assert validate_bundle_root(output, workflows=["cold_step_adaptive"]).case_id == case


@pytest.mark.parametrize("damage, message", [("checksum", "sha256"), ("escape", "resolves outside")])
def test_validation_rejects_changed_or_external_artifacts(data, tmp_path, damage, message):
    case, source, output = data
    create_bundle(case, source, output)
    mesh = output / "inputs/mesh.msh"
    if damage == "checksum":
        mesh.write_bytes(b"X" * mesh.stat().st_size)
    else:
        outside = tmp_path / "outside.msh"
        mesh.rename(outside)
        mesh.symlink_to(outside)
    with pytest.raises(BundleError, match=message):
        validate_bundle_root(output)


def test_optional_file_becomes_required_for_selected_workflow(data):
    case, source, output = data
    create_bundle(case, source, output)
    (output / "inputs/param_cold_fixed_time_init.txt").unlink()
    assert validate_bundle_root(output).warnings
    with pytest.raises(BundleError, match="cold_fixed_time_init_parameters"):
        validate_bundle_root(output, workflows=["cold_step_adaptive"])


def test_failed_creation_leaves_existing_data_and_no_partial_bundle(data):
    case, source, output = data
    (source / "geometry.geo").unlink()
    with pytest.raises(BundleError, match="required files missing.*geometry"):
        create_bundle(case, source, output)
    assert not output.exists()
    output.symlink_to(output.parent / "absent")
    with pytest.raises(BundleError, match="output already exists"):
        create_bundle(case, source, output)
    assert output.is_symlink()
    output.unlink()
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("keep")
    with pytest.raises(BundleError, match="output already exists"):
        create_bundle(case, source, output)
    assert marker.read_text() == "keep"
