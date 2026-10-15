"""Exercise refresh/publication with real files, preparation and comparisons."""

import json
import shutil
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import h5py
import pytest

from regression_tests import cli, golden
from regression_tests.bundles import create_bundle, validate_bundle_root, series_artifact_ids
from regression_tests.catalog import load_case_definition
from regression_tests.files import file_identity
from regression_tests.support import BundleError
from regression_tests.tests.fixtures.harness import create_harness

SOLVER = '''#!/usr/bin/env bash
set -euo pipefail
if [[ -e inputs/restart.h5 ]]; then
  cp inputs/restart.h5 outputs/result.h5
else
  cp "$(dirname "$0")/seed.h5" outputs/result.h5
fi
printf 'Error: 1.0E-8\\nOutput written to file outputs/result.h5\\n'
'''


@pytest.fixture
def setup(tmp_path, monkeypatch):
    harness = create_harness(tmp_path, solver=SOLVER)
    root = harness.catalog
    workflows = json.loads((root / 'workflows.json').read_text())
    for name in ('warm', 'cold_fixed'):
        workflows['workflows'][name]['layout'] = 'serial_omp1'
    workflows['workflows']['cold_fixed'].pop('reference')
    (root / 'workflows.json').write_text(json.dumps(workflows))
    (root / 'suites.json').write_text(json.dumps({
        'schema_version': 2, 'defaults': {'case': 'legacy_case', 'layout': 'serial_omp1'},
        'suites': {'warm': {'description': 'Recheck collected reference', 'workflows': ['warm']}}}))
    recipe = {'cases': {'legacy_case': {'producers': ['warm'], 'checks': ['warm']}}}
    (root / 'golden.json').write_text(json.dumps(recipe))
    monkeypatch.setattr(golden, 'ROOT', root)
    monkeypatch.setattr(cli, 'ROOT', root)
    with h5py.File(harness.serial_executable.parent / 'seed.h5', 'r+') as handle:
        handle['solution/u'][...] += 1
    build = Mock(return_value=SimpleNamespace(metadata_path=harness.build_manifest))
    monkeypatch.setattr(golden, 'build_solver', build)
    return SimpleNamespace(harness=harness, root=root, build=build,
                           values=harness.values, workspace=tmp_path / 'refresh', output=tmp_path / 'golden')


def publish(setup):
    return golden.publish(setup.workspace, setup.output, 'reviewed-1', 'Accepted solver change', 'Review record 42')


def test_refresh_collects_a_checkpoint_separately_from_the_final_reference(setup):
    import shlex
    import sys

    path = setup.root / 'cases/legacy_case.json'
    case = json.loads(path.read_text())
    case['files']['optional']['step6_checkpoint'] = 'restart_step6.h5'
    case['workflows']['warm']['checkpoint_outputs'] = {'step6_checkpoint': 6}
    path.write_text(json.dumps(case))
    checkpoint = f'''{shlex.quote(sys.executable)} - <<'PY'
import h5py, shutil
shutil.copy2('outputs/result.h5', 'outputs/checkpoint.h5')
with h5py.File('outputs/checkpoint.h5', 'r+') as handle:
    handle['simulation_parameters/time/Current_time_step_number'] = 6
    handle['solution/u'][...] += 6
PY
printf 'Error: 1.0E-8\\nOutput written to file outputs/checkpoint.h5\\n'
'''
    solver = SOLVER.replace("printf 'Error:", checkpoint + "printf 'Error:")
    setup.harness.install_solver(solver, 'serial')
    _, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    candidate = setup.workspace / 'candidate'
    manifest = json.loads((candidate / 'manifest.json').read_text())
    record = manifest['artifacts'][manifest['roles']['step6_checkpoint']]
    produced = Path(report['producers'][0]['run_directory']) / 'outputs/checkpoint.h5'
    assert file_identity(candidate / record['path']) == file_identity(produced)
    with h5py.File(candidate / record['path']) as checkpoint, h5py.File(candidate / 'inputs/reference_mpi4_omp4.h5') as final:
        assert checkpoint['simulation_parameters/time/Current_time_step_number'][()] == 6
        assert (checkpoint['solution/u'][()] == final['solution/u'][()] + 6).all()
    publish(setup)
    assert file_identity(setup.output / record['path']) == file_identity(produced)


def mesh_producers(setup):
    path = setup.root / 'workflows.json'
    workflows = json.loads(path.read_text())
    workflows['workflows']['cold_adaptive']['outputs'] = ['warm_restart', 'mesh']
    workflows['workflows']['cold_fixed']['outputs'] = []
    path.write_text(json.dumps(workflows))
    (setup.root / 'golden.json').write_text(json.dumps({
        'cases': {'legacy_case': {'producers': ['cold_adaptive', 'cold_fixed', 'warm'],
                                  'checks': ['warm']}}}))
    # Known files from a tiny executable; this fixture does not simulate meshing.
    return SOLVER + '''
if [[ "$PWD" == */cold_adaptive/* && "$PWD" != */03_final ]]; then
  printf '%s\\n' "${PWD##*/}" > res/temp.msh
fi
'''


def test_refresh_handoff_review_and_explicit_self_contained_publication(setup):
    setup.harness.install_solver(mesh_producers(setup), 'serial')
    recipe_path = setup.root / 'golden.json'
    recipe = json.loads(recipe_path.read_text())
    producers = recipe['cases']['legacy_case']['producers']
    producers[0], producers[1] = producers[1], producers[0]
    recipe_path.write_text(json.dumps(recipe))
    with pytest.raises(BundleError, match='producer cold_fixed needs: mesh'):
        golden.refresh('legacy_case', setup.values, setup.workspace)
    assert not setup.workspace.exists()
    producers[0], producers[1] = producers[1], producers[0]
    recipe_path.write_text(json.dumps(recipe))
    manifest_path = setup.harness.bundle / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    artifact = manifest['artifacts'].pop(manifest['roles'].pop('warm_restart'))
    (setup.harness.bundle / artifact['path']).unlink()
    # Historical artifacts must not become requirements of the new candidate.
    manifest['roles']['retired_restart'] = 'retired_restart'
    manifest['artifacts']['retired_restart'] = dict(manifest['artifacts'][manifest['roles']['warm_reference']])
    manifest['artifacts']['retired_restart']['path'] = 'inputs/retired.h5'
    shutil.copy2(setup.harness.bundle / 'inputs/reference_mpi4_omp4.h5', setup.harness.bundle / 'inputs/retired.h5')
    manifest_path.write_text(json.dumps(manifest))
    original = manifest_path.read_bytes()
    machine = setup.workspace.parent / 'settings.json'
    machine.write_text(json.dumps({'defaults': {'bundles': {'legacy_case': str(setup.harness.bundle)}}}))
    assert cli.main(['golden', 'refresh', 'legacy_case', '--settings', str(machine),
                     '--workspace', str(setup.workspace), '--jobs', '2']) == 0
    report = json.loads((setup.workspace / 'refresh.json').read_text())
    assert report['status'] == 'ready' and not setup.output.exists()
    setup.build.assert_called_once()
    assert setup.build.call_args.kwargs['requirements'] == {('NGammaTiTeNeutral', 'serial')}
    assert report['producers'][0]['old_reference']['status'] == 'unavailable'
    assert report['producers'][-1]['old_reference']['status'] == 'failed'
    warm = Path(report['producers'][-1]['run_directory'])
    restart = warm / 'inputs/restart.h5'
    cold = Path(report['producers'][0]['run_directory'])
    assert restart.is_symlink() and restart.resolve().is_relative_to(cold)
    assert not (warm / 'inputs/reference.h5').exists()
    fixed = Path(report['producers'][-2]['run_directory'])
    mesh = fixed / 'stages/01_initial/inputs/mesh.msh'
    assert mesh.is_symlink() and mesh.resolve() == cold / 'stages/02_continued/res/temp.msh'
    assert not (fixed / 'stages/01_initial/inputs/restart.h5').exists()
    assert report['producers'][-2]['mesh_handoff']['status'] == 'passed'
    # The same bytes are protected until publication, not merely linked once.
    original_mesh = mesh.read_bytes()
    mesh.write_bytes(b'changed mesh')
    with pytest.raises(BundleError, match='evidence changed'):
        publish(setup)
    mesh.write_bytes(original_mesh)
    result = cli.main(['golden', 'publish', str(setup.workspace), '--output', str(setup.output),
                      '--bundle-version', 'reviewed-1', '--reason', 'Accepted solver change', '--provenance', 'Review record 42'])
    assert result == 0
    with pytest.raises(BundleError, match='output already exists'):
        publish(setup)
    assert manifest_path.read_bytes() == original
    assert not any(item.is_symlink() for item in setup.output.rglob('*'))
    manifest = json.loads((setup.output / 'manifest.json').read_text())
    assert manifest['bundle_class'] == 'golden'
    assert not (setup.output / 'provenance/refresh').exists()
    assert 'retired_restart' not in manifest['roles']
    assert 'retired_restart' not in manifest['artifacts']
    assert not (setup.output / 'inputs/retired.h5').exists()
    reference = setup.output / manifest['artifacts'][manifest['roles']['warm_reference']]['path']
    assert reference.read_bytes() == (warm / "outputs/result.h5").read_bytes()
    published_mesh = setup.output / manifest['artifacts'][manifest['roles']['mesh']]['path']
    assert published_mesh.read_bytes() == original_mesh
    receipt = json.loads((setup.output / 'provenance/publication.json').read_text())
    assert receipt['reason'] == 'Accepted solver change' and receipt['provenance'] == 'Review record 42'
    assert receipt['build']['repository']['revision'] == 'fixture-revision'
    assert receipt['refresh_report'] == file_identity(setup.workspace / 'refresh.json')
    shutil.rmtree(setup.workspace)
    shutil.rmtree(setup.harness.bundle)
    assert validate_bundle_root(setup.output, setup.root / 'cases').case_id == 'legacy_case'
    # An ordinary fixed run uses only the published bundle, after refresh is gone.
    from regression_tests.prepare import prepare_run
    values = {**setup.values, 'MHDG_REGRESSION_DATA_ROOT': str(setup.output)}
    prepared = prepare_run(values, 'legacy_case', 'cold_fixed', 'serial_omp1',
                           setup.root / 'cases', setup.root / 'layouts.json', require_reference=False)
    assert (prepared.stages[0].run.path / 'inputs/mesh.msh').resolve() == published_mesh
    assert not (prepared.stages[0].run.path / 'inputs/restart.h5').exists()


def test_missing_generated_mesh_is_rejected(tmp_path):
    (tmp_path / 'run_metadata.json').write_text(json.dumps({
        'stages': [{'run_directory': str(tmp_path / 'stage')}],
    }))
    with pytest.raises(BundleError, match='retained no res/temp.msh'):
        golden._generated_mesh(tmp_path)


@pytest.mark.parametrize('failure', ['exit', 'nonconvergence', 'validation', 'endpoint'])
def test_failed_producer_or_validation_blocks_publication(setup, failure):
    solver = SOLVER
    if failure == 'exit':
        solver = '#!/usr/bin/env bash\nexit 7\n'
    elif failure == 'nonconvergence':
        solver = SOLVER.replace('1.0E-8', '1.0E5')
    elif failure == 'validation':
        solver = SOLVER.replace('set -euo pipefail', 'set -euo pipefail\n[[ "$PWD" != *-verify ]] || exit 7')
    else:
        mesh_solver = mesh_producers(setup)
        seed = setup.harness.serial_executable.parent / 'mismatched_result.h5'
        shutil.copy2(seed.with_name('seed.h5'), seed)
        with h5py.File(seed, 'r+') as handle:
            handle['solution/u'][...] += 0.1
        solver = mesh_solver + '''
if [[ "$PWD" == */cold_fixed/* ]]; then
  cp "$(dirname "$0")/mismatched_result.h5" outputs/result.h5
fi
'''
    setup.harness.install_solver(solver, 'serial')
    with pytest.raises(BundleError, match='producer failed|invalid producer output|golden validation failed|endpoint comparison failed'):
        golden.refresh('legacy_case', setup.values, setup.workspace)
    report = json.loads((setup.workspace / 'refresh.json').read_text())
    assert report['status'] == 'failed'
    if failure == 'validation':
        assert all(item['status'] == 'passed' for item in report['producers'])
        assert report['checks'][0]['status'] == 'failed'
    else:
        assert any(item['status'] == 'failed' for item in report['producers'])
        assert not (setup.workspace / 'candidate').exists()
    with pytest.raises(BundleError, match='successful, validated refresh'):
        publish(setup)
    assert not setup.output.exists()


def test_publish_rejects_changed_outputs_evidence_and_missing_reason(setup):
    path = setup.root / 'workflows.json'
    workflows = json.loads(path.read_text())
    workflows['workflows']['warm']['transport_overrides'] = {'diff_n_min_phys': 0.1}
    path.write_text(json.dumps(workflows))
    _, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    with pytest.raises(BundleError, match='reason and provenance'):
        golden.publish(setup.workspace, setup.output, 'v1', ' ', 'Review')
    refresh_path = setup.workspace / 'refresh.json'
    original_report = refresh_path.read_bytes()
    report['producers'][0]['validation']['convergence']['passed'] = False
    refresh_path.write_text(json.dumps(report))
    with pytest.raises(BundleError, match='successful recorded validation'):
        publish(setup)
    refresh_path.write_bytes(original_report)
    candidate = setup.workspace / 'candidate'
    manifest_path = candidate / 'manifest.json'
    original = manifest_path.read_bytes()
    manifest_path.write_bytes(original + b'\n')
    with pytest.raises(BundleError, match='evidence changed'):
        publish(setup)
    manifest_path.write_bytes(original)
    run = Path(report['producers'][-1]['run_directory'])
    namelist = run / 'inputs/transport_model.nml'
    original = namelist.read_bytes()
    namelist.write_bytes(original + b'\n')
    with pytest.raises(BundleError, match='evidence changed'):
        publish(setup)
    namelist.write_bytes(original)
    (run / 'outputs/result.h5').write_bytes(b'changed')
    with pytest.raises(BundleError, match='outputs changed'):
        publish(setup)
    assert not setup.output.exists()


@pytest.fixture
def trajectory_setup(setup):
    source = setup.harness.source
    with h5py.File(source / 'trajectory.h5', 'w') as handle:
        handle['time'] = [0., .02, .04]
        handle['grid/r'], handle['grid/z'] = [2., 3.], [-1., 1.]
        handle['equilibrium/psiSep'] = [1., 1., 1.]
        for name in ('psi', 'Bphi', 'Br', 'Bz', 'Jtor'):
            handle[f'equilibrium/{name}'] = [[[n, n+1], [n+2, n+3]] for n in range(3)]
    with h5py.File(source / 'restart.h5', 'r+') as handle:
        handle['simulation_parameters/adimensionalization/time_scale'] = 2e-7
    path = setup.root / 'cases/legacy_case.json'
    case = json.loads(path.read_text())
    case['files']['required']['equilibrium_series'] = 'trajectory.h5'
    case['workflows']['moving'] = {'extends': 'warm',
        'moving_equilibrium': {'series': 'equilibrium_series'}, 'parameter_overrides': {'nts': 1}}
    path.write_text(json.dumps(case))
    setup.source_bundle = setup.harness.root / 'trajectory_bundle'
    create_bundle('legacy_case', source, setup.source_bundle, setup.root / 'cases')
    setup.values = {**setup.values, 'MHDG_REGRESSION_DATA_ROOT': str(setup.source_bundle)}
    # An unrelated unroled artifact must still be removed by compact collection.
    path = setup.source_bundle / 'manifest.json'
    manifest = json.loads(path.read_text())
    retired = setup.source_bundle / 'inputs/retired.h5'
    retired.write_bytes(b'retired field')
    manifest['artifacts']['retired_field'] = {'path': 'inputs/retired.h5',
        **file_identity(retired), 'media_type': 'application/x-hdf5'}
    path.write_text(json.dumps(manifest))
    return setup


def test_trajectory_survives_refresh_publication_and_source_removal(trajectory_setup):
    setup = trajectory_setup
    source_manifest = json.loads((setup.source_bundle / 'manifest.json').read_text())
    case = load_case_definition('legacy_case', setup.root / 'cases')
    artifact_id = source_manifest['roles']['equilibrium_series']
    assert not series_artifact_ids(setup.source_bundle, source_manifest, case)
    golden.refresh('legacy_case', setup.values, setup.workspace)
    publish(setup)
    for bundle in (setup.workspace / 'candidate', setup.output):
        manifest = json.loads((bundle / 'manifest.json').read_text())
        assert manifest['artifacts'][artifact_id] == source_manifest['artifacts'][artifact_id]
        assert file_identity(bundle / manifest['artifacts'][artifact_id]['path']) == {
            key: source_manifest['artifacts'][artifact_id][key] for key in ('sha256', 'size_bytes')}
        assert 'retired_field' not in manifest['artifacts']
        assert not (bundle / 'inputs/retired.h5').exists()
    assert not any(path.is_symlink() for path in setup.output.rglob('*'))
    for directory in (setup.workspace, setup.source_bundle, setup.harness.source, setup.harness.bundle):
        shutil.rmtree(directory)
    validate_bundle_root(setup.output, setup.root / 'cases', workflows=['moving'])
    from regression_tests.prepare import prepare_run
    values = {**setup.values, 'MHDG_REGRESSION_DATA_ROOT': str(setup.output)}
    prepared = prepare_run(values, 'legacy_case', 'moving', 'serial_omp1',
                           setup.root / 'cases', setup.root / 'layouts.json', require_reference=False)
    with h5py.File(prepared.path / 'inputs/equilibrium.h5') as handle:
        assert handle['time'][:].tolist() == [0., .02]
        assert handle['equilibrium/Jtor'].shape == (2, 2, 2)


@pytest.mark.parametrize('damage', ['missing', 'changed'])
def test_publication_rejects_damaged_trajectory(trajectory_setup, damage):
    setup = trajectory_setup
    golden.refresh('legacy_case', setup.values, setup.workspace)
    candidate = setup.workspace / 'candidate'
    manifest = json.loads((candidate / 'manifest.json').read_text())
    artifact_id = manifest['roles']['equilibrium_series']
    path = candidate / manifest['artifacts'][artifact_id]['path']
    if damage == 'missing':
        path.unlink()
    else:
        path.write_bytes(b'changed field')
    with pytest.raises(BundleError, match='missing|size_bytes|sha256|does not exist'):
        publish(setup)
    assert not setup.output.exists()


@pytest.fixture
def stopped_setup(trajectory_setup):
    setup = trajectory_setup
    source = setup.harness.source
    path = setup.root / 'cases/legacy_case.json'
    case = json.loads(path.read_text())
    workflow = case['workflows']['moving']
    workflow['moving_equilibrium']['expected_failure'] = {
        'reason': 'divergence_refinement_limit', 'completed_steps': 1}
    workflow['parameter_overrides'].update(nts=2, max_divergence_refinements=0)
    workflow['parameter_namelists'] = {'max_divergence_refinements': 'adapt_lst'}
    path.write_text(json.dumps(case))
    (setup.root / 'golden.json').write_text(json.dumps({
        'cases': {'legacy_case': {'producers': ['moving'], 'checks': ['moving']}}}))
    (setup.root / 'suites.json').write_text(json.dumps({
        'schema_version': 2, 'defaults': {'case': 'legacy_case', 'layout': 'serial_omp1'},
        'suites': {'moving': {'description': 'Check expected stop', 'workflows': ['moving']}}}))
    for number in range(3):
        path = setup.harness.serial_executable.parent / f'step{number}.h5'
        shutil.copy2(path.with_name('seed.h5'), path)
        with h5py.File(path, 'r+') as handle:
            handle['simulation_parameters/adimensionalization/time_scale'] = 2e-7
            for key, value in {'Current_time_step_number': number, 'Current_time': number * .02 / 2e-7,
                               'Current_time_step': .02 / 2e-7}.items():
                handle[f'simulation_parameters/time/{key}'] = value
    # The old reference has the same completed time, but intentionally different fields.
    with h5py.File(source / 'reference_mpi4_omp4.h5', 'r+') as handle:
        handle['simulation_parameters/adimensionalization/time_scale'] = 2e-7
        for key, value in {'Current_time_step_number': 1, 'Current_time': .02 / 2e-7,
                           'Current_time_step': .02 / 2e-7}.items():
            handle[f'simulation_parameters/time/{key}'] = value
    setup.source_bundle = setup.harness.root / 'stopped_bundle'
    create_bundle('legacy_case', source, setup.source_bundle, setup.root / 'cases')
    setup.values['MHDG_REGRESSION_DATA_ROOT'] = str(setup.source_bundle)
    setup.harness.install_solver('''#!/usr/bin/env bash
set -euo pipefail
cp "$(dirname "$0")/step0.h5" outputs/result_0000.h5
printf 'Output written to file outputs/result_0000.h5\\nTime iteration = 1\\nNR iteration: 1\\nError: 1.0E-8\\n'
cp "$(dirname "$0")/step1.h5" outputs/result_0001.h5
printf 'Output written to file outputs/result_0001.h5\\nTime iteration = 2\\nNR iteration: 1\\nError: 1.0E9\\n'
cp "$(dirname "$0")/step2.h5" outputs/result_NR0001.h5
printf 'Output written to file outputs/result_NR0001.h5\\nNR refinement recovery stopped: per-timestep refinement limit reached: 0\\n'
exit 1
''', 'serial')
    return setup


@pytest.mark.parametrize('terminal', ['divergence_refinement_limit', 'oscillation_refinement_limit'])
def test_expected_stop_refresh_collects_only_last_completion_and_publishes(stopped_setup, terminal):
    setup = stopped_setup
    if terminal == 'oscillation_refinement_limit':
        path = setup.root / 'cases/legacy_case.json'
        case = json.loads(path.read_text())
        workflow = case['workflows']['moving']
        workflow['moving_equilibrium']['expected_failure']['reason'] = 'refinement_limit'
        workflow['parameter_overrides']['max_oscillation_refinements'] = 1
        workflow['parameter_namelists']['max_oscillation_refinements'] = 'adapt_lst'
        path.write_text(json.dumps(case))
        solver = setup.harness.serial_executable.read_text().replace(
            'NR refinement recovery stopped: per-timestep refinement limit reached: 0',
            'Oscillation refinement stopped: per-timestep refinement limit reached: 1')
        setup.harness.install_solver(solver, 'serial')
    _, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    producer = report['producers'][0]
    run = Path(producer['run_directory'])
    assert report['status'] == 'ready' and producer['run_status'] == 'solver_failed'
    assert producer['old_reference']['status'] == 'failed'  # Compared, rather than marked unavailable.
    assert golden._solution(run).name == 'result_NR0001.h5'
    manifest = json.loads((setup.workspace / 'candidate/manifest.json').read_text())
    reference = setup.workspace / 'candidate' / manifest['artifacts'][manifest['roles']['warm_reference']]['path']
    assert reference.read_bytes() == (run / 'outputs/result_0001.h5').read_bytes()
    checks = json.loads((setup.workspace / report['checks'][0]['summary']).read_text())
    assert checks['status'] == 'passed' and checks['results'][0]['run_status'] == 'solver_failed'
    publish(setup)
    receipt = json.loads((setup.output / 'provenance/publication.json').read_text())
    assert receipt['producers'][0]['expected_failure'] == {
        'reason': terminal, 'completed_steps': [1], 'last_attempted_step': 2}
    assert receipt['producers'][0]['run_status'] == 'solver_failed'
    assert not list(setup.output.rglob('*NR*.h5'))
    assert json.loads((run / 'run_metadata.json').read_text())['status'] == 'solver_failed'


@pytest.mark.parametrize('target', ['producer', 'check'])
def test_expected_stop_publication_rejects_changed_checkpoint(stopped_setup, target):
    setup = stopped_setup
    _, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    result = report['producers'][0] if target == 'producer' else json.loads(
        (setup.workspace / report['checks'][0]['summary']).read_text())['results'][0]
    (Path(result['run_directory']) / 'outputs/result_0001.h5').write_bytes(b'changed checkpoint')
    with pytest.raises(BundleError, match='outputs changed'):
        publish(setup)
    assert not setup.output.exists()


@pytest.mark.parametrize('damage', ['crash', 'reason', 'completion'])
def test_expected_stop_refresh_rejects_unexpected_outcome(stopped_setup, damage):
    setup = stopped_setup
    solver = setup.harness.serial_executable.read_text()
    if damage == 'crash':
        solver = solver.replace('exit 1', 'exit 139')
    else:
        solver = solver.replace('NR refinement recovery stopped: per-timestep refinement limit reached: 0', 'unexpected termination')
        if damage == 'completion':
            solver = solver.replace('exit 1', 'exit 0')
    setup.harness.install_solver(solver, 'serial')
    with pytest.raises(BundleError, match='invalid producer output'):
        golden.refresh('legacy_case', setup.values, setup.workspace)
    assert not (setup.workspace / 'candidate').exists()
