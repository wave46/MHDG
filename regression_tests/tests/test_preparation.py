"""Preparation behavior, without duplicating every production workflow declaration."""

import json

import pytest

from regression_tests.prepare import prepare_run
from regression_tests.parameters import render_parameter_file, read_selected_input_values
from regression_tests.catalog import load_case_definition
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.harness import create_harness

@pytest.fixture
def harness(tmp_path):
    return create_harness(tmp_path)


def prepare(harness, workflow="warm", run_id="probe", **options):
    return prepare_run(harness.values, "legacy_case", workflow, "mpi4_omp4",
                       harness.catalog / "cases", harness.catalog / "layouts.json", run_id, **options)


def test_warm_links_immutable_inputs_and_renders_isolated_paths(harness):
    source = harness.bundle / "inputs/param.txt"
    original = source.read_bytes()
    run = prepare(harness)
    assert run.command[-2] == str(run.path / "inputs/mesh")  # Required CLI argument, unused file.
    assert not (run.path / "inputs/mesh.msh").exists()
    assert run.command[-1] == str(run.path / "inputs/restart")
    assert (run.path / "inputs/equilibrium.h5").is_symlink()
    assert (run.path / "positionFeketeNodesTri2D.h5").resolve() == harness.runtime_file
    rendered = (run.path / "param.txt").read_text()
    assert str(run.path / "inputs/equilibrium.h5") in rendered
    assert str(run.path / "outputs") in rendered
    assert source.read_bytes() == original
    assert not (run.path / "param.txt").is_symlink()
    with pytest.raises(BundleError, match="already exists"):
        prepare(harness)
    assert (run.path / "param.txt").read_text() == rendered


def test_restart_and_reference_roles_are_selected_independently(harness, monkeypatch):
    case = load_case_definition("legacy_case", harness.catalog / "cases")
    workflow = case["workflows"]["warm"]
    workflow.update(restart="warm_reference", reference="warm_restart")
    monkeypatch.setattr("regression_tests.prepare.load_case_definition", lambda *_, **__: case)
    run = prepare(harness)
    assert (run.path / "inputs/restart.h5").resolve() == harness.bundle / "inputs/reference_mpi4_omp4.h5"
    assert (run.path / "inputs/reference.h5").resolve() == harness.bundle / "inputs/restart.h5"

    manifest_path = harness.bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["roles"]["warm_restart"]  # Reference-only role in this swapped workflow.
    manifest_path.write_text(json.dumps(manifest))
    producer = prepare(harness, run_id="producer", require_reference=False)
    assert (producer.path / "inputs/restart.h5").resolve() == harness.bundle / "inputs/reference_mpi4_omp4.h5"
    assert not (producer.path / "inputs/reference.h5").exists()
    with pytest.raises(BundleError, match="does not provide role warm_restart"):
        prepare(harness, run_id="needs-reference")


@pytest.mark.parametrize("adaptive", [False, True])
def test_cold_preparation_preserves_mesh_and_restart_policy(harness, adaptive):
    workflow = "cold_adaptive" if adaptive else "cold_fixed"
    run = prepare(harness, workflow)
    assert [stage.stage_id for stage in run.stages] == ["initial", "continued", "final"]
    assert run.stages[0].restart_from == "analytical"
    assert all(stage.restart_from == "previous_stage" for stage in run.stages[1:])
    mesh = "mesh_adaptive_initial.msh" if adaptive else "mesh.msh"
    assert (run.stages[0].run.path / "inputs/mesh.msh").resolve() == harness.bundle / "inputs" / mesh
    assert not (run.stages[0].run.path / "inputs/transport_model.nml").exists()
    for stage in run.stages[1:]:
        assert not (stage.run.path / "inputs/mesh.msh").exists()
        assert (stage.run.path / "inputs/transport_model.nml").is_symlink()
        assert not (stage.run.path / "inputs/restart.h5").exists()
    first = (run.stages[0].run.path / "param.txt").read_text()
    assert f"rest_adapt = .{'true' if adaptive else 'false'}." in first


def test_disabled_features_prepare_without_optional_files(harness):
    path = harness.catalog / "workflows.json"
    catalog = json.loads(path.read_text())
    catalog["workflows"]["minimal"] = {
        "extends": "warm", "transport": None, "impurity_configuration": None,
        "parameters": "initial_parameters",
        "parameter_overrides": {"transport_1d": False, "impurity_radiation": False},
    }
    path.write_text(json.dumps(catalog))
    path = harness.catalog / "cases/legacy_case.json"
    case = json.loads(path.read_text())
    case["workflows"]["minimal"] = {}
    path.write_text(json.dumps(case))
    for name in ("mesh.msh", "transport_model.nml", "impurity_model_w.nml"):
        (harness.bundle / "inputs" / name).unlink()
    template = harness.root / "minimal.txt"
    template.write_text("".join(line for line in (harness.bundle / "inputs/param_cold_fixed_time_init.txt").read_text().splitlines(True)
                                if "transport_model_path" not in line and "impurity_model_path" not in line))
    run = prepare(harness, "minimal", artifact_overrides={"initial_parameters": template})
    assert not any((run.path / "inputs" / name).exists() for name in
                   ("mesh.msh", "transport_model.nml", "impurity_model.nml"))
    assert (run.path / "inputs/restart.h5").is_symlink()
    assert json.loads((run.path / "run_plan.json").read_text())["artifacts"]["initial_parameters"] == str(template)
    with pytest.raises(BundleError, match="no declared input"):
        prepare(harness, "minimal", run_id="enabled", requested_overrides={"transport_1d": True})
    with pytest.raises(BundleError, match="declare a mesh"):
        prepare(harness, "minimal", run_id="mesh", requested_overrides={"readMeshFromSol": False})


def test_analytical_start_cannot_read_a_nonexistent_restart(harness):
    with pytest.raises(BundleError, match="requires a restart"):
        prepare(harness, "cold_fixed", requested_overrides={"readMeshFromSol": True})


def test_failed_render_does_not_publish_a_partial_run(harness):
    source = harness.bundle / "inputs/param.txt"
    source.write_text(source.read_text() + "field_path = 'duplicate'\n")
    with pytest.raises(BundleError, match="appear once"):
        prepare(harness, validate_bundle=False)
    parent = harness.run_directory("warm", "mpi4_omp4", "probe").parent
    assert not list(parent.iterdir())


def test_parameter_renderer_uses_declared_namelists_for_new_scalars(tmp_path):
    source, target = tmp_path / "input.nml", tmp_path / "output.nml"
    source.write_text("&existing\n label = 'old!text' ! retain comment\n untouched = 3\n/\n&new_section\n/\n")
    render_parameter_file(source, target, {},
                          {"label": "new!text", "new_knob": 0.25, "enabled": True},
                          {"new_knob": "new_section", "enabled": "new_section"})
    text = target.read_text()
    assert "label = 'new!text' ! retain comment" in text
    assert " untouched = 3\n" in text
    assert "new_knob = 0.25" in text.split("&new_section")[1]
    assert "enabled = .true." in text.split("&new_section")[1]
    selected = read_selected_input_values(target, {"New_Knob", "ENABLED"})
    assert selected == {"new_knob": 0.25, "enabled": True}


@pytest.mark.parametrize("text,namelists", [
    ("&n\n x = 1\n x = 2\n/\n", {}),
    ("&n\n/\n", {}),
    ("&n\n/\n", {"x": "absent_group"}),
    ("&n\n x = 1, y = 2\n/\n", {}),
])
def test_ambiguous_or_undeclared_parameter_edits_fail(tmp_path, text, namelists):
    source, target = tmp_path / "input.nml", tmp_path / "output.nml"
    source.write_text(text)
    with pytest.raises(BundleError):
        render_parameter_file(source, target, {}, {"x": 3}, namelists)
    assert not target.exists()
