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


@pytest.fixture
def me_harness(harness):
    import h5py
    import numpy as np
    from regression_tests.bundles import create_bundle

    source = harness.source
    (source / "magnetic").mkdir()
    frames = []
    for n in range(26):
        frame = {}
        for kind in ("magnetic", "current"):
            path = source / f"magnetic/{kind}_{n + 1:04d}.h5"
            with h5py.File(path, "w") as handle:
                handle["time"], handle["dt"] = float(n * .02), .02
                handle["field"] = np.arange(4).reshape(2, 2) + n * 10
            frame[kind] = path.relative_to(source).as_posix()
        frames.append(frame)
    (source / "series.json").write_text(json.dumps({"schema_version": 1, "dt_s": .02, "frames": frames}))
    for name in ("target_density", "impurity_concentration"):
        with h5py.File(source / f"{name}.h5", "w") as handle:
            handle["time"], handle[name] = [0., .5, .52], [2., 10., 10.]
    with h5py.File(source / "restart.h5", "r+") as handle:
        handle["simulation_parameters/adimensionalization/time_scale"] = 2e-7
    (source / "impurity_model_w.nml").write_text("&IMPURITY_RADIATION_LST\n"
        " impurity_concentration_path = '/old/concentration.h5'\n impurity_concentration_dimension = 3\n/\n")
    render_parameter_file(source / "param.txt", source / "param_me.txt", {"target_density_path": "/old/density.h5"},
                          {"target_density_dimension": 3},
                          {"target_density_path": "input_lst", "target_density_dimension": "input_lst"})
    path = harness.catalog / "cases/legacy_case.json"
    case = json.loads(path.read_text())
    case["files"]["optional"].update(equilibrium_series="series.json", me_parameters="param_me.txt",
                                    density="target_density.h5", concentration="impurity_concentration.h5")
    case["workflows"]["warm"].update(parameters="me_parameters", moving_equilibrium={
        "series": "equilibrium_series", "time_step_s": .08,
        "controls": {"target_density": "density", "impurity_concentration": "concentration"}})
    path.write_text(json.dumps(case))
    create_bundle("legacy_case", source, harness.root / "me_bundle", harness.catalog / "cases")
    settings = json.loads(harness.settings.read_text())
    settings["defaults"]["bundles"]["legacy_case"] = "me_bundle"
    harness.settings.write_text(json.dumps(settings))
    return harness


@pytest.mark.parametrize("dt,steps,mapping", [
    (.08, None, [1, 5, 9, 13, 17, 21, 25, 26]),
    (.1, None, [1, 6, 11, 16, 21, 26]), (.02, 3, [1, 2, 3, 4]),
])
def test_me_preparation_retimes_pairs_and_reuses_canonical_controls(me_harness, dt, steps, mapping):
    import h5py
    import numpy as np
    from regression_tests.files import file_identity

    path = me_harness.catalog / "cases/legacy_case.json"
    case = json.loads(path.read_text())
    case["workflows"]["warm"]["moving_equilibrium"]["time_step_s"] = dt
    path.write_text(json.dumps(case))
    bundle = me_harness.root / "me_bundle"
    original = {str(p): file_identity(p) for p in bundle.rglob("*") if p.is_file()}
    run = prepare(me_harness, requested_overrides={"nts": steps} if steps else {})
    plan = json.loads((run.path / "run_plan.json").read_text())
    assert plan["moving_equilibrium_inputs"]["canonical_frames"] == mapping
    assert plan["moving_equilibrium_inputs"]["source_manifest"] == file_identity(bundle / "manifest.json")
    values = read_selected_input_values(run.path / "param.txt", {"nts", "dt0", "tfi", "field_path", "jtor_path", "target_density_dimension"})
    assert values["nts"] == len(mapping) - 1
    assert values["dt0"] * 2e-7 == pytest.approx(dt)
    assert values["tfi"] * 2e-7 == pytest.approx((len(mapping) - 1) * dt)
    assert values["field_path"] == str(run.path / "inputs/magnetic/equilibrium")
    assert values["jtor_path"] == str(run.path / "inputs/magnetic/current_density")
    for name, identity in plan["moving_equilibrium_inputs"]["generated_files"].items():
        assert file_identity(run.path / name) == identity
        number = int(name[-7:-3])
        with h5py.File(run.path / name) as handle:
            assert handle["time"][()] == pytest.approx((number - 1) * dt)
            assert handle["dt"][()] == dt
            np.testing.assert_array_equal(handle["field"][()], np.arange(4).reshape(2, 2) + (mapping[number - 1] - 1) * 10)
    assert values["target_density_dimension"] == 3
    assert (run.path / "inputs/target_density.h5").resolve() == bundle / "inputs/target_density.h5"
    assert (run.path / "inputs/restart.h5").is_symlink()
    impurity = read_selected_input_values(run.path / "inputs/impurity_model.nml", {"impurity_concentration_path", "impurity_concentration_dimension"})
    assert impurity == {"impurity_concentration_path": str(run.path / "inputs/impurity_concentration.h5"), "impurity_concentration_dimension": 3}
    assert original == {str(p): file_identity(p) for p in bundle.rglob("*") if p.is_file()}


@pytest.mark.parametrize("failure", ["timestep", "metadata", "normalized_override", "control_range", "control_dimension"])
def test_me_preparation_failure_is_atomic(me_harness, failure):
    import h5py

    options = {"validate_bundle": False}
    if failure == "timestep":
        path = me_harness.catalog / "cases/legacy_case.json"
        case = json.loads(path.read_text())
        case["workflows"]["warm"]["moving_equilibrium"]["time_step_s"] = .07
        path.write_text(json.dumps(case))
    elif failure == "metadata":
        with h5py.File(me_harness.root / "me_bundle/inputs/equilibrium_series/magnetic/current_0001.h5", "r+") as handle:
            handle["time"][...] = .3
    elif failure == "control_range":
        with h5py.File(me_harness.root / "me_bundle/inputs/target_density.h5", "r+") as handle:
            handle["time"][...] = [0., .1, .2]
    elif failure == "control_dimension":
        options["requested_overrides"] = {"target_density_dimension": 27}
    else:
        options["requested_overrides"] = {"dt0": 123.}
    with pytest.raises(BundleError):
        prepare(me_harness, **options)
    parent = me_harness.run_directory("warm", "mpi4_omp4", "probe").parent
    assert not parent.exists() or not list(parent.iterdir())
