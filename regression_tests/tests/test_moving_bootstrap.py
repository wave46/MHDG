"""Moving runs share a pinned state and generate their own anchored controls."""

import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest

from regression_tests import config, golden
from regression_tests.bundles import create_bundle, preflight_bundles
from regression_tests.catalog import load_selection
from regression_tests.files import file_identity
from regression_tests.moving_bootstrap import bootstrap_pin, resolve_bootstrap, measure_initial_line_density
from regression_tests.parameters import read_selected_input_values
from regression_tests.prepare import prepare_run
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.case_data import PARAMETERS
from regression_tests.tests.fixtures.harness import create_harness
from regression_tests.tests.fixtures.solutions import write_solution


ROOT = Path(__file__).resolve().parents[1]


def control(path, name, times, values):
    with h5py.File(path, 'w') as handle:
        handle['time'], handle[name] = np.asarray(times, dtype=float), np.asarray(values, dtype=float)
        if name == 'target_density':
            for key, value in zip(('x_lower', 'x_upper', 'y_lower', 'y_upper'), (-.2, 1.2, .25, .25)):
                handle[key] = [value]


@pytest.fixture
def data(tmp_path):
    harness = create_harness(tmp_path)
    limited = tmp_path / 'limited_source'
    limited.mkdir()
    for name in ('geometry.geo', 'equilibrium.h5'):
        (limited / name).write_text('supplied input')
    seed = limited / 'reference_bootstrap.h5'
    write_solution(seed, model=True)
    with h5py.File(seed, 'r+') as handle:
        values = handle['solution/u'][()]
        values[0::5] = 2.
        handle['solution/u'][...] = values
        handle['simulation_parameters/adimensionalization/length_scale'] = 1.901e-3
        handle['simulation_parameters/adimensionalization/density_scale'] = 1e19
        handle['simulation_parameters/adimensionalization/time_scale'] = 2e-7
        handle['simulation_parameters/physics/puff'] = 9e20
    limited_bundle = tmp_path / 'limited'
    create_bundle('legacy_case', limited, limited_bundle)

    source = tmp_path / 'moving_source'
    source.mkdir()
    (source / 'geometry.geo').write_text('geometry')
    parameters = PARAMETERS.replace('readMeshFromSol = .false.', 'readMeshFromSol = .true.').replace(
        'transport_1d = .true.', 'transport_1d = .false.').replace(
        '&SWITCH_LST\n', '&SWITCH_LST\n target_variable = 1\n').replace(
        '&PHYS_LST\n', '&PHYS_LST\n puff = 5e20\n feedback_propotional_gain = 1\n'
        ' feedback_integral_gain = 0\n feedback_derivative_gain = 0\n').replace(
        '&ADAPT_LST\n', '&ADAPT_LST\n check_neutral_positivity = .true.\n')
    (source / 'param.txt').write_text(parameters)
    (source / 'impurity_model.nml').write_text('&IMPURITY_RADIATION_LST\n/\n')
    (source / 'bootstrap.json').write_text(json.dumps(bootstrap_pin(limited_bundle, ROOT / 'cases')))
    with h5py.File(source / 'equilibrium.h5', 'w') as handle:
        handle['time'] = np.arange(69)*.02
        handle['grid/r'], handle['grid/z'] = [0., 1.], [0., 1.]
        handle['equilibrium/psiSep'] = np.ones(69)
        for name in ('psi', 'Bphi', 'Br', 'Bz', 'Jtor'):
            handle[f'equilibrium/{name}'] = np.ones((69, 2, 2))
    control(source / 'puff.h5', 'puff', [0., .5, 1.36], [5e20, 5e21, 5e21])
    control(source / 'target_density.h5', 'target_density', [0., .5, 1.36], [2e19, 8e19, 8e19])
    control(source / 'short_puff.h5', 'puff', [0., .05, .13, .21], [5e20, 9.5e20, 1.67e21, 2.39e21])
    control(source / 'short_target_density.h5', 'target_density', [0., .05, .13, .21], [2e19, 2.6e19, 3.56e19, 4.52e19])
    control(source / 'ramp_target_density.h5', 'target_density', [1.36, 2.3], [8e19, 9e19])
    control(source / 'zeff.h5', 'zeff', [0., 1.36], [1., 1.])
    control(source / 'impurity_concentration.h5', 'concentration', [0., 1.36], [1e-4, 1e-4])
    shutil.copy2(seed, source / 'restart_mode1_step6.h5')
    shutil.copy2(seed, source / 'reference_mode1.h5')
    with h5py.File(source / 'restart_mode1_step6.h5', 'r+') as handle:
        handle['solution/u'][...] *= 100
        handle['simulation_parameters/time/Current_time_step_number'] = 6
        handle['simulation_parameters/time/Current_time'] = .48 / 2e-7
        handle['simulation_parameters/switches/ME'] = 1
    moving = tmp_path / 'moving'
    create_bundle('moving_equilibrium', source, moving)
    values = {**harness.values, 'MHDG_REGRESSION_DATA_ROOT': str(moving),
              'MHDG_REGRESSION_BOOTSTRAP_BUNDLE': str(limited_bundle)}
    return harness, limited_bundle, moving, values


def prepare(data, workflow, name):
    return prepare_run(data[3], 'moving_equilibrium', workflow, 'mpi4_omp4',
                       ROOT / 'cases', ROOT / 'layouts.json', name, require_reference=False)


def test_main_short_and_restart_controls_share_the_original_anchor(data):
    _, limited, bundle, _ = data
    original = {path.name: file_identity(path) for path in (bundle / 'inputs').iterdir()}
    main = prepare(data, 'mode1', 'main')
    short = prepare(data, 'mode2_short', 'short')
    restart = prepare(data, 'mode1_restart', 'restart')
    assert (main.path / 'inputs/restart.h5').resolve() == limited / 'inputs/reference_bootstrap.h5'
    assert (restart.path / 'inputs/restart.h5').resolve() == bundle / 'inputs/restart_mode1_step6.h5'
    assert not (bundle / 'inputs/restart.h5').exists()
    measured = measure_initial_line_density(limited / 'inputs/reference_bootstrap.h5', bundle / 'inputs/target_density.h5')
    with h5py.File(main.path / 'inputs/target_density.h5') as full:
        assert full['target_density'][0] == 2*measured
        assert full['target_density'][-1] == 8e19
        np.testing.assert_array_equal(full['time'][()], [0., .5, 1.36])
        with h5py.File(short.path / 'inputs/target_density.h5') as small:
            np.testing.assert_array_equal(small['time'][()], [0., .05, .13, .21])
            np.testing.assert_allclose(small['target_density'][()], np.interp(
                small['time'][()], full['time'][()], full['target_density'][()]))
    for name in ('target_density.h5', 'puff.h5'):
        assert file_identity(main.path / 'inputs' / name) == file_identity(restart.path / 'inputs' / name)
    parameters = read_selected_input_values(main.path / 'param.txt', ['puff'])
    assert parameters['puff'] == 9e20
    assert original == {path.name: file_identity(path) for path in (bundle / 'inputs').iterdir()}
    plan = json.loads((main.path / 'run_plan.json').read_text())
    for name, identity in plan['moving_equilibrium_inputs']['generated_files'].items():
        assert file_identity(main.path / name) == identity


def test_wrong_bootstrap_is_rejected_before_a_run_is_created(data):
    harness, limited, bundle, values = data
    descriptor = bundle / 'inputs/bootstrap.json'
    state = limited / 'inputs/reference_bootstrap.h5'
    with h5py.File(state, 'r+') as handle:
        handle['solution/u'][...] *= 2
    manifest = json.loads((limited / 'manifest.json').read_text())
    manifest['artifacts'][manifest['roles']['bootstrap_reference']].update(file_identity(state))
    (limited / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(BundleError, match='differs from the pinned'):
        prepare(data, 'mode1', 'bad')
    assert not harness.run_root.exists()


def test_publication_or_relocation_can_preserve_the_pinned_state(data, tmp_path):
    _, limited, moving, values = data
    relocated = tmp_path / 'published_limited'
    shutil.copytree(limited, relocated)
    manifest = json.loads((relocated / 'manifest.json').read_text())
    manifest.update(bundle_class='golden', bundle_version='reviewed')
    (relocated / 'manifest.json').write_text(json.dumps(manifest))
    selected = {**values, 'MHDG_REGRESSION_BOOTSTRAP_BUNDLE': str(relocated)}
    state, _ = resolve_bootstrap(moving / 'inputs/bootstrap.json', selected, ROOT / 'cases')
    assert state == relocated / 'inputs/reference_bootstrap.h5'


def test_preflight_and_settings_select_the_dependency_explicitly(data, tmp_path):
    _, limited, moving, values = data
    _, checks = load_selection('moving_equilibrium_parallel', ROOT / 'suites.json', ROOT / 'layouts.json', ROOT / 'cases')
    preflight_bundles(checks, {'moving_equilibrium': values}, ROOT / 'cases')
    machine = tmp_path / 'machine.json'
    machine.write_text(json.dumps({'defaults': {'bundles': {'legacy_case': str(limited)}}}))
    resolved = config.settings(machine, case='moving_equilibrium', bundle=moving, use_build=False)
    assert resolved['MHDG_REGRESSION_BOOTSTRAP_BUNDLE'] == str(limited)
    with pytest.raises(BundleError, match='select the limited bootstrap'):
        preflight_bundles(checks, {'moving_equilibrium': {'MHDG_REGRESSION_DATA_ROOT': str(moving)}}, ROOT / 'cases')


def test_refresh_selects_new_bootstrap_and_collects_only_its_descriptor(data, tmp_path, monkeypatch):
    from types import SimpleNamespace
    from regression_tests.bundles import load_manifest
    from regression_tests.catalog import load_case_definition

    harness, limited, moving, values = data
    new_limited = tmp_path / 'new_limited'
    shutil.copytree(limited, new_limited)
    seed = new_limited / 'inputs/reference_bootstrap.h5'
    with h5py.File(seed, 'r+') as handle:
        handle['solution/u'][...] *= 2
    manifest = json.loads((new_limited / 'manifest.json').read_text())
    manifest['artifacts'][manifest['roles']['bootstrap_reference']].update(file_identity(seed))
    (new_limited / 'manifest.json').write_text(json.dumps(manifest))
    selected = {**values, 'MHDG_REGRESSION_BOOTSTRAP_BUNDLE': str(new_limited)}
    original = file_identity(moving / 'inputs/bootstrap.json')
    monkeypatch.setattr(golden, 'build_solver', lambda *a, **k: SimpleNamespace(metadata_path=harness.build_manifest))
    captured = []
    def stop_before_execution(prepared, settings):
        captured.append(prepared)
        raise KeyboardInterrupt
    monkeypatch.setattr(golden, 'execute_prepared', stop_before_execution)
    workspace = tmp_path / 'refresh'
    with pytest.raises(KeyboardInterrupt):
        golden.refresh('moving_equilibrium', selected, workspace)
    report = json.loads((workspace / 'refresh.json').read_text())
    assert report['status'] == 'failed'
    descriptor = workspace / report['bootstrap_dependency']['path']
    assert json.loads(descriptor.read_text())['solution'] == file_identity(seed)
    assert (captured[0].path / 'inputs/restart.h5').resolve() == seed
    assert original == file_identity(moving / 'inputs/bootstrap.json')
    candidate = workspace / 'candidate'
    case = load_case_definition('moving_equilibrium', ROOT / 'cases')
    golden._collect(moving, candidate, load_manifest(moving), {'bootstrap': descriptor}, report, case)
    assert file_identity(candidate / 'inputs/bootstrap.json') == file_identity(descriptor)
    assert not (candidate / 'inputs/restart.h5').exists()
