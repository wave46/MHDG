"""Generate both warm references and validate them with a separate observer."""

import json
from pathlib import Path
import shutil
import sys
from unittest.mock import Mock

import pytest

from regression_tests import golden, golden_continue, suites
from regression_tests.catalog import load_selection
from regression_tests.files import file_identity
from regression_tests.support import BundleError
from regression_tests.tests.test_golden import setup


@pytest.fixture
def campaign(setup):
    path = setup.root / 'workflows.json'
    document = json.loads(path.read_text())
    document['parameter_namelists'].update(newton_formulation='numer_lst', residual_diagnostics_mode='utils_lst')
    document['workflows']['warm']['parameter_overrides']['newton_formulation'] = 'absolute'
    document['workflows']['incremental'] = {
        'extends': 'warm', 'reference': 'incremental_reference', 'outputs': ['incremental_reference'],
        'parameter_overrides': {'newton_formulation': 'incremental'},
        'comparison': {'method': 'fixed_hdf5', 'profile': 'feature_converged'}}
    document['workflows']['ordinary'] = {'extends': 'warm', 'outputs': []}
    path.write_text(json.dumps(document))
    path = setup.root / 'cases/legacy_case.json'
    document = json.loads(path.read_text())
    document['files']['optional']['incremental_reference'] = 'reference_incremental.h5'
    document['workflows']['incremental'] = {}
    document['workflows']['ordinary'] = {}
    path.write_text(json.dumps(document))
    path = setup.root / 'suites.json'
    document = json.loads(path.read_text())
    document['suites']['newton_formulations'] = {
        'description': 'Both warm formulations', 'workflows': ['warm', 'incremental'],
        'workflow_comparison': ['warm', 'incremental'], 'newton_formulation_comparison': True,
        'tolerance_profile': 'newton_formulations', 'layouts': ['serial_omp1']}
    document['suites']['ordinary'] = {'description': 'Ordinary production check', 'workflows': ['ordinary']}
    document['profiles'] = {'combined': {'description': 'Observer then ordinary check',
                                        'checks': [{'suite': 'ordinary'}, {'suite': 'newton_formulations'}]}}
    path.write_text(json.dumps(document))
    (setup.root / 'golden.json').write_text(json.dumps({'cases': {'legacy_case': {
        'producers': ['warm', 'incremental'], 'checks': ['newton_formulations', 'ordinary']}}}))

    program = '#!' + sys.executable + '\n' + SOLVER
    setup.harness.install_solver(program, 'serial')
    normal = json.loads(setup.harness.build_manifest.read_text())
    observed_directory = setup.harness.serial_executable.parent / 'observer'
    observed_directory.mkdir()
    for name in ('serial', 'parallel', 'seed.h5', setup.harness.runtime_file.name):
        shutil.copy2(setup.harness.serial_executable.parent / name, observed_directory / name)
    for name in ('serial', 'parallel'):
        executable = observed_directory / name
        executable.write_text(program.replace('OBSERVE = False', 'OBSERVE = True'))
    observed = {**normal, 'artifacts': {
        key: {'path': item['path'], **file_identity(observed_directory / item['path'])}
        for key, item in normal['artifacts'].items()},
        'newton_probe': {'kind': 'mixed_linear_residual_observation'}}
    observed_manifest = observed_directory / 'build_metadata.json'
    observed_manifest.write_text(json.dumps(observed))
    normal['newton_probe'] = {'manifest': {'path': str(observed_manifest), **file_identity(observed_manifest)}}
    setup.harness.build_manifest.write_text(json.dumps(normal))
    setup.observed = observed_directory
    return setup


def test_refresh_saves_both_references_and_limits_observation_to_pair(campaign):
    setup = campaign
    _, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    assert report['status'] == 'ready'
    assert setup.build.call_args.kwargs['newton_probe'] == {('NGammaTiTeNeutral', 'serial')}
    assert [item['workflow_id'] for item in report['producers']] == ['warm', 'incremental']
    for producer in report['producers']:
        log = Path(producer['run_directory']) / 'stdout.log'
        assert 'Newton linear check' not in log.read_text()
    manifest = json.loads((setup.workspace / 'candidate/manifest.json').read_text())
    assert {'warm_reference', 'incremental_reference'} <= manifest['roles'].keys()
    pair = json.loads((setup.workspace / report['checks'][0]['summary']).read_text())
    assert pair['parameter_overrides']['saveNR'] is True
    assert pair['parameter_overrides']['residual_diagnostics_mode'] == 'equations'
    assert pair['reference_comparisons'] is True
    for result in pair['results']:
        assert result['comparison_status'] == 'passed'
        assert 'Newton linear check' in (Path(result['run_directory']) / 'stdout.log').read_text()
        comparison = json.loads(Path(result['comparison_report']).read_text())
        assert comparison['status'] == 'passed'
        role, profile = ('warm_reference', 'fixed_same_layout') if result['workflow_id'] == 'warm' else (
            'incremental_reference', 'feature_converged')
        expected = setup.workspace / 'candidate' / manifest['artifacts'][manifest['roles'][role]]['path']
        assert Path(comparison['reference']) == expected
        assert comparison['tolerance_profile']['id'] == profile
    cross = json.loads(Path(pair['comparisons'][0]['comparison_report']).read_text())
    assert cross['status'] == 'passed'
    assert cross['tolerance_profile']['id'] == 'newton_formulations'
    assert pair['comparisons'][0]['newton']['histories']['absolute']['update_errors'] == [1e-5]
    ordinary = json.loads((setup.workspace / report['checks'][1]['summary']).read_text())
    assert ordinary['execution_inputs']['executables'] != pair['execution_inputs']['executables']
    assert 'Newton linear check' not in (Path(ordinary['results'][0]['run_directory']) / 'stdout.log').read_text()


def test_continuation_observes_incremental_rerun_and_retains_absolute(campaign, monkeypatch):
    setup = campaign
    failure = setup.observed / 'bad_linear'
    failure.touch()
    with pytest.raises(BundleError, match='golden validation failed: newton_formulations'):
        golden.refresh('legacy_case', setup.values, setup.workspace)
    first = json.loads((setup.workspace / 'refresh.json').read_text())
    old_pair = json.loads((setup.workspace / first['checks'][0]['summary']).read_text())
    failure.unlink()
    monkeypatch.setattr(golden_continue, 'build_solver', setup.build)
    _, report = golden_continue.continue_refresh(setup.workspace, setup.values, ['incremental'], catalog_root=setup.root)
    assert report['status'] == 'ready'
    assert setup.build.call_args.kwargs['newton_probe'] == {('NGammaTiTeNeutral', 'serial')}
    pair = json.loads((setup.workspace / report['checks'][0]['summary']).read_text())
    assert pair['results'][0]['reused_execution'] is True
    assert pair['results'][0]['run_directory'] == old_pair['results'][0]['run_directory']
    assert pair['results'][1]['run_directory'] != old_pair['results'][1]['run_directory']


def test_profile_rejects_missing_observer_before_running_any_case(campaign, monkeypatch):
    setup = campaign
    build = json.loads(setup.harness.build_manifest.read_text())
    build.pop('newton_probe')
    setup.harness.build_manifest.write_text(json.dumps(build))
    _, checks = load_selection('combined', setup.root / 'suites.json', setup.root / 'layouts.json',
                               setup.root / 'cases')
    # This preflight only needs the restart inputs; references are produced later.
    for check in checks:
        check['reference_comparisons'] = False
    execute = Mock(side_effect=AssertionError('must not launch a workflow'))
    monkeypatch.setattr(suites, 'run_suite', execute)
    with pytest.raises(BundleError, match='requires a linear-observation test build'):
        suites.run_profile('combined', checks, {'legacy_case': setup.values}, setup.root, 'preflight')
    execute.assert_not_called()


SOLVER = '''import pathlib, re, shutil
import h5py
OBSERVE = False
parameters = pathlib.Path('param.txt').read_text()
formulation = re.search(r"newton_formulation\\s*=\\s*'([^']+)'", parameters, re.I).group(1)
shutil.copy2('inputs/restart.h5', 'outputs/result.h5')
with h5py.File('outputs/result.h5', 'r+') as handle:
    handle['simulation_parameters/numerics/newton_formulation'] = formulation
if OBSERVE:
    for name in ('state_0000.h5', 'state_NR0001.h5'):
        shutil.copy2('outputs/result.h5', 'outputs/' + name)
    bad = formulation == 'incremental' and (pathlib.Path(__file__).parent/'bad_linear').exists()
    relative = 0.5 if bad else 1e-10
    print('HDG residual before update: norm= 1.0; target= 1e-8; within target=no')
    print(f'Newton linear check: initial=1.0; trial={relative}; relative={relative}')
print('Error: 1.0E-5')
print('Output written to file outputs/result.h5')
'''
