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
    assert read_selected_input_values(run.path / "param.txt", {"jtor_path"})["jtor_path"] == ""
    assert not (run.path / "inputs/current_density.h5").exists()
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
    with h5py.File(source / "trajectory.h5", "w") as handle:
        handle["time"] = np.arange(29)*.02
        handle["grid/r"], handle["grid/z"] = [2., 3.], [-1., 1.]
        handle["equilibrium/psiSep"] = np.arange(29)+40.
        handle["equilibrium/cocos"] = [12]
        for name in ("psi", "Bphi", "Br", "Bz", "Jtor"):
            handle[f"equilibrium/{name}"] = np.array([
                np.arange(4).reshape(2, 2)+min(n, 25)*10 for n in range(29)])
        handle["equilibrium/Br"].attrs["units"] = "T"
    for name in ("target_density", "impurity_concentration"):
        with h5py.File(source / f"{name}.h5", "w") as handle:
            handle["time"], handle[name] = [0., .5, .56], [2., 10., 10.]
    with h5py.File(source / "restart.h5", "r+") as handle:
        handle["simulation_parameters/adimensionalization/time_scale"] = 2e-7
    (source / "impurity_model_w.nml").write_text("&IMPURITY_RADIATION_LST\n"
        " impurity_concentration_path = '/old/concentration.h5'\n impurity_concentration_dimension = 3\n/\n")
    render_parameter_file(source / "param.txt", source / "param_me.txt", {"target_density_path": "/old/density.h5"},
                          {"target_density_dimension": 3},
                          {"target_density_path": "input_lst", "target_density_dimension": "input_lst"})
    path = harness.catalog / "cases/legacy_case.json"
    case = json.loads(path.read_text())
    case["files"]["optional"].update(equilibrium_series="trajectory.h5", me_parameters="param_me.txt",
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
    (.08, None, [1, 5, 9, 13, 17, 21, 25, 29]),
    (.1, None, [1, 6, 11, 16, 21, 26]), (.02, 3, [1, 2, 3, 4]),
])
def test_me_preparation_selects_trajectory_endpoints_and_reuses_canonical_controls(me_harness, dt, steps, mapping):
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
    assert values["field_path"] == str(run.path / "inputs/equilibrium.h5")
    assert values["jtor_path"] == ""
    for name, identity in plan["moving_equilibrium_inputs"]["generated_files"].items():
        assert file_identity(run.path / name) == identity
        with h5py.File(run.path / name) as handle:
            np.testing.assert_allclose(handle["time"][:], np.arange(len(mapping))*dt)
            assert handle["equilibrium/psi"].shape == (len(mapping), 2, 2)
            assert handle["equilibrium/Br"].attrs["units"] == "T"
            np.testing.assert_array_equal(handle["equilibrium/cocos"][:], [12])
            for number, frame in enumerate(mapping):
                for component in ("psi", "Jtor"):
                    np.testing.assert_array_equal(handle[f"equilibrium/{component}"][number],
                        np.arange(4).reshape(2, 2)+min(frame-1, 25)*10)
    assert values["target_density_dimension"] == 3
    assert (run.path / "inputs/target_density.h5").resolve() == bundle / "inputs/target_density.h5"
    assert (run.path / "inputs/restart.h5").is_symlink()
    impurity = read_selected_input_values(run.path / "inputs/impurity_model.nml", {"impurity_concentration_path", "impurity_concentration_dimension"})
    assert impurity == {"impurity_concentration_path": str(run.path / "inputs/impurity_concentration.h5"), "impurity_concentration_dimension": 3}
    assert original == {str(p): file_identity(p) for p in bundle.rglob("*") if p.is_file()}


@pytest.mark.parametrize("bad_time", [False, True])
def test_me_preparation_restores_checkpoint_clock(me_harness, bad_time):
    import h5py

    with h5py.File(me_harness.root / "me_bundle/inputs/restart.h5", "r+") as handle:
        handle['simulation_parameters/switches/ME'] = 1
        handle['simulation_parameters/time/Current_time_step_number'] = 1
        handle['simulation_parameters/time/Current_time'] = (.1 if bad_time else .08) / 2e-7
    options = {"validate_bundle": False, "requested_overrides": {"nts": 3}}
    if bad_time:
        with pytest.raises(BundleError, match="restart time differs"):
            prepare(me_harness, **options)
    else:
        run = prepare(me_harness, **options)
        plan = json.loads((run.path / 'run_plan.json').read_text())
        assert plan['moving_equilibrium_inputs']['initial_step'] == 1
        assert plan['parameter_overrides']['nts'] == 3


@pytest.mark.parametrize("failure", ["timestep", "metadata", "initial_time", "normalized_override", "control_range"])
def test_me_preparation_failure_is_atomic(me_harness, failure):
    import h5py

    options = {"validate_bundle": False}
    if failure == "timestep":
        path = me_harness.catalog / "cases/legacy_case.json"
        case = json.loads(path.read_text())
        case["workflows"]["warm"]["moving_equilibrium"]["time_step_s"] = .07
        path.write_text(json.dumps(case))
    elif failure in {"metadata", "initial_time"}:
        with h5py.File(me_harness.root / "me_bundle/inputs/trajectory.h5", "r+") as handle:
            if failure == "metadata":
                handle["time"][2] = .3
            else:
                handle["time"][...] += .02
    elif failure == "control_range":
        with h5py.File(me_harness.root / "me_bundle/inputs/target_density.h5", "r+") as handle:
            handle["time"][...] = [0., .1, .2]
    else:
        options["requested_overrides"] = {"dt0": 123.}
    with pytest.raises(BundleError):
        prepare(me_harness, **options)
    parent = me_harness.run_directory("warm", "mpi4_omp4", "probe").parent
    assert not parent.exists() or not list(parent.iterdir())


@pytest.mark.parametrize("fixed", [False, True])
def test_physical_preparation_preserves_intermediate_clock_and_supports_constant_controls(me_harness, fixed):
    import h5py
    import numpy as np

    path = me_harness.catalog / "cases/legacy_case.json"
    case = json.loads(path.read_text())
    moving = case["workflows"]["warm"]["moving_equilibrium"]
    moving["time_step_s"] = .08 if fixed else .1
    initial, final = (1.36, 2.3) if fixed else (.225, .3)
    if fixed:
        moving["equilibrium_time_s"] = .56
        moving["controls"] = {"target_density": "density"}
        with h5py.File(me_harness.root / "me_bundle/inputs/target_density.h5", "r+") as handle:
            handle["time"][...] = [1.36, 1.8, 2.4]
    path.write_text(json.dumps(case))
    with h5py.File(me_harness.root / "me_bundle/inputs/restart.h5", "r+") as handle:
        handle["simulation_parameters/switches/ME"] = 1
        handle["simulation_parameters/time/Current_time_step_number"] = 9
        handle["simulation_parameters/time/Current_time"] = initial / 2e-7
        handle["simulation_parameters/time/Current_time_seconds"] = initial
    run = prepare(me_harness, validate_bundle=False, requested_overrides={
        "time_mode": "adaptive" if fixed else "prescribed", "t_final": final})
    record = json.loads((run.path / "run_plan.json").read_text())["moving_equilibrium_inputs"]
    assert record["initial_step"] == 9 and record["initial_time_s"] == initial
    assert record["final_time_s"] == final
    values = read_selected_input_values(run.path / "param.txt", {"time_units", "dt0", "tfi", "t_final", "me", "tis"})
    assert values == {"time_units": "seconds", "dt0": moving["time_step_s"], "tfi": final,
                      "t_final": final, "me": not fixed, "tis": 1}
    with h5py.File(run.path / "inputs/equilibrium.h5") as handle:
        np.testing.assert_allclose(handle["time"][:], [.56] if fixed else [0., .1, .2, .3])
        assert handle["equilibrium/psi"].shape == (1 if fixed else 4, 2, 2)
    if fixed:
        assert read_selected_input_values(run.path / "inputs/impurity_model.nml", {"impurity_concentration_path"}) == {
            "impurity_concentration_path": ""}
        assert read_selected_input_values(run.path / "param.txt", {"puff_path", "zeff_path"}) == {
            "puff_path": "", "zeff_path": ""}


def test_physical_preparation_starts_steady_seed_at_configured_start(me_harness):
    import h5py

    with h5py.File(me_harness.root / "me_bundle/inputs/restart.h5", "r+") as handle:
        handle["simulation_parameters/time/Physical_time_evolution"] = 0
        handle["simulation_parameters/time/Current_time_step_number"] = 70
        handle["simulation_parameters/time/Current_time"] = 1e10
    run = prepare(me_harness, validate_bundle=False, requested_overrides={
        "time_mode": "adaptive", "t_start": .12, "t_final": .3})
    record = json.loads((run.path / "run_plan.json").read_text())["moving_equilibrium_inputs"]
    assert record["initial_step"] == 0 and record["initial_time_s"] == .12
