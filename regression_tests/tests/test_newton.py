"""Reject misleading Newton comparisons even when final fields agree."""
import json
from pathlib import Path

import h5py
import pytest

from regression_tests.newton import compare_newton_runs
from regression_tests.suites import run_suite
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.solutions import write_solution

CATALOG = Path(__file__).resolve().parents[1]
LOG = ('HDG residual before update: norm= 1.0; target= 1e-8; within target=no\n'
       'Newton linear check: initial=1.0; trial=1e-10; relative=1e-10\n'
       'Error: 1e-2\nError: 1e-5\n')


def test_regular_build_is_rejected_before_any_newton_runs(tmp_path):
    manifest = tmp_path/'build.json'
    manifest.write_text('{}')
    settings = {'MHDG_BUILD_MANIFEST': str(manifest)}
    with pytest.raises(BundleError, match='requires a linear-observation test build'):
        run_suite(settings, 'newton_formulations', CATALOG/'cases', CATALOG/'layouts.json',
                  CATALOG/'suites.json', CATALOG/'tolerances.json')


@pytest.fixture
def pair(tmp_path):
    runs = [tmp_path/'absolute', tmp_path/'incremental']
    for run in runs:
        (run/'outputs').mkdir(parents=True)
        (run/'inputs').mkdir()
        for name in ('restart.h5', 'geometry.geo', 'equilibrium.h5'):
            (run/'inputs'/name).write_text('identical input')
        plan = {'layout': {'execution': 'mpi', 'mpi_ranks': 4, 'omp_threads': 4},
                'parameter_overrides': {'newton_formulation': run.name, 'stab': 2}}
        (run/'run_plan.json').write_text(json.dumps(plan))
        (run/'stdout.log').write_text(LOG)
        write_solution(run/'outputs/state_0000.h5', model=True)
        write_solution(run/'outputs/state_NR0001.h5', model=True, solution_offset=.1)
        with h5py.File(run/'outputs/state_NR0001.h5', 'r+') as handle:
            handle['simulation_parameters/numerics/newton_formulation'] = run.name
    return runs


def compare(pair):
    return compare_newton_runs(*pair, CATALOG/'tolerances.json', 'newton_formulations')


def test_complete_update_and_linear_observations(pair):
    report = compare(pair)
    assert report['status'] == 'passed', report['failures']
    assert all(report['updates'][name]['passed'] for name in ('u', 'q', 'u_tilde'))
    assert report['histories']['incremental']['linear'][0]['trial'] == 1e-10


@pytest.mark.parametrize('update_error, offset', [('1e-5', 1e-7), ('0.0', 0.0)])
def test_one_step_warm_reconvergence_is_accepted(pair, update_error, offset):
    for run in pair:
        (run/'stdout.log').write_text(LOG.replace('Error: 1e-2\nError: 1e-5\n',
                                                f'Error: {update_error}\n'))
        write_solution(run/'outputs/state_NR0001.h5', model=True, solution_offset=offset)
        with h5py.File(run/'outputs/state_NR0001.h5', 'r+') as handle:
            handle['simulation_parameters/numerics/newton_formulation'] = run.name
    report = compare(pair)
    assert report['status'] == 'passed', report['failures']
    assert len(report['histories']['absolute']['update_errors']) == 1
    assert len(report['histories']['incremental']['update_errors']) == 1


@pytest.mark.parametrize('field', ['q', 'u_tilde'])
def test_corrupted_gradient_or_trace_update_is_rejected(pair, field):
    with h5py.File(pair[1]/'outputs/state_NR0001.h5', 'r+') as handle:
        handle['solution/'+field][0] += .5
    report = compare(pair)
    assert report['status'] == 'failed'
    assert not report['updates'][field]['passed']


def test_equal_updates_from_different_initial_states_are_rejected(pair):
    for filename in ('state_0000.h5', 'state_NR0001.h5'):
        with h5py.File(pair[1]/'outputs'/filename, 'r+') as handle:
            handle['solution/u'][0] += .01
    assert any('different assembled initial u' in failure for failure in compare(pair)['failures'])


def test_large_energy_scale_cannot_hide_a_neutral_update_error(pair):
    for run in pair:
        for filename in ('state_0000.h5', 'state_NR0001.h5'):
            with h5py.File(run/'outputs'/filename, 'r+') as handle:
                values = handle['solution/u'][()]
                values[2::5] *= 1e12
                handle['solution/u'][...] = values
    with h5py.File(pair[1]/'outputs/state_NR0001.h5', 'r+') as handle:
        handle['solution/u'][4] += .5
    report = compare(pair)
    assert report['updates']['u']['relative_l2'] < 1e-8
    assert not report['updates']['u']['passed']


@pytest.mark.parametrize('replacement', ['', LOG.replace('relative=1e-10', 'relative=0.5'),
                                       LOG.replace('trial=1e-10', 'trial=nan'),
                                       LOG.replace('Error: 1e-2\nError: 1e-5\n', '')])
def test_missing_or_invalid_observations_are_rejected(pair, replacement):
    (pair[1]/'stdout.log').write_text(replacement)
    assert compare(pair)['status'] == 'failed'


def test_different_solver_controls_are_rejected(pair):
    path = pair[1]/'run_plan.json'
    plan = json.loads(path.read_text())
    plan['parameter_overrides']['stab'] = 1
    path.write_text(json.dumps(plan))
    assert any('different numerical/physics controls' in failure for failure in compare(pair)['failures'])
