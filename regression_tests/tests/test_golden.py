"""Exercise refresh/publication with real files, preparation and comparisons."""

import json
import shutil
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import h5py
import pytest

from regression_tests import cli, golden
from regression_tests.bundles import validate_bundle_root
from regression_tests.support import BundleError
from regression_tests.files import file_identity
from regression_tests.tests.fixtures.harness import create_harness, REGRESSION_ROOT
from regression_tests.tests.fixtures.solutions import write_solution

SOLVER = '''#!/usr/bin/env bash
set -euo pipefail
if [[ -e inputs/restart.h5 ]]; then
  cp inputs/restart.h5 outputs/result.h5
else
  cp "$(dirname "$0")/seed.h5" outputs/result.h5
fi
printf 'Error: 1.0E-8\\nOutput written to file outputs/result.h5\\n'
'''


def write_off_solution(path, offset=0):
    write_solution(path, solution_offset=offset)
    with h5py.File(path, 'r+') as handle:
        handle['simulation_parameters/switches/balance_diagnostics_mode'] = b'off'


@pytest.fixture
def setup(tmp_path, monkeypatch):
    root = tmp_path / 'catalog'
    root.mkdir()
    shutil.copytree(REGRESSION_ROOT / 'schemas', root / 'schemas')
    shutil.copytree(REGRESSION_ROOT / 'cases', root / 'cases')
    for name in ('workflows.json', 'tolerances.json'):
        shutil.copy2(REGRESSION_ROOT / name, root / name)
    workflows = json.loads((root / 'workflows.json').read_text())
    for name in ('warm', 'cold_fixed'):
        workflows['workflows'][name]['layout'] = 'serial_omp1'
    (root / 'workflows.json').write_text(json.dumps(workflows))
    (root / 'layouts.json').write_text(json.dumps({'schema_version': 2, 'layouts': ['serial_omp1', 'serial_omp16']}))
    (root / 'suites.json').write_text(json.dumps({
        'schema_version': 2, 'defaults': {'case': 'legacy_case', 'layout': 'serial_omp1'},
        'suites': {'warm': {'description': 'Recheck collected reference', 'workflows': ['warm']}}}))
    recipe = {'matrix': {'relations': ['all_pairs'], 'tolerance_profile': 'cold_cross_layout',
                         'layout_comparison_policy': 'fixed_hdf5'},
              'cases': {'legacy_case': {'producers': [{'workflow': 'cold_fixed', 'matrix': True}, 'warm'],
                                        'checks': ['warm']}}}
    (root / 'golden.json').write_text(json.dumps(recipe))
    monkeypatch.setattr(golden, 'ROOT', root)
    harness = create_harness(tmp_path, solver=SOLVER, reference_writer=write_off_solution)
    write_off_solution(harness.serial_executable.parent / 'seed.h5', 1)
    build_path = harness.build_manifest
    build = Mock(return_value=SimpleNamespace(metadata_path=build_path))
    monkeypatch.setattr(golden, 'build_solver', build)
    return SimpleNamespace(harness=harness, root=root, build=build, build_path=build_path,
                           values=harness.values, workspace=tmp_path / 'refresh', output=tmp_path / 'golden')


def publish(setup):
    return golden.publish(setup.workspace, setup.output, 'reviewed-1', 'Accepted solver change', 'Review record 42')


@pytest.mark.parametrize('has_old_reference', [True, False])
def test_refresh_handoff_review_and_explicit_self_contained_publication(setup, has_old_reference):
    manifest_path = setup.harness.bundle / 'manifest.json'
    if not has_old_reference:
        manifest = json.loads(manifest_path.read_text())
        for role in ('warm_reference', 'warm_restart'):
            artifact = manifest['artifacts'].pop(manifest['roles'].pop(role))
            (setup.harness.bundle / artifact['path']).unlink()
        manifest_path.write_text(json.dumps(manifest))
    original = manifest_path.read_bytes()
    if has_old_reference:
        machine = setup.workspace.parent / 'settings.json'
        machine.write_text(json.dumps({'defaults': {'bundles': {'legacy_case': str(setup.harness.bundle)}}}))
        assert cli.main(['golden', 'refresh', 'legacy_case', '--settings', str(machine),
                         '--workspace', str(setup.workspace), '--jobs', '2']) == 0
        path = setup.workspace / 'refresh.json'
        report = json.loads(path.read_text())
    else:
        path, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    assert path.is_file() and report['status'] == 'ready' and not setup.output.exists()
    setup.build.assert_called_once()
    assert setup.build.call_args.kwargs['variants'] == {'serial'}
    assert report['parallel_checks'] and all(item['status'] == 'passed' for item in report['parallel_checks'])
    assert report['producers'][-1]['old_reference']['status'] == ('failed' if has_old_reference else 'unavailable')
    warm = Path(report['producers'][-1]['run_directory'])
    restart = warm / 'inputs/restart.h5'
    cold = Path(report['producers'][0]['run_directory'])
    assert restart.is_symlink() and restart.resolve().is_relative_to(cold)
    assert not (warm / 'inputs/reference.h5').exists()
    result = cli.main(['golden', 'publish', str(setup.workspace), '--output', str(setup.output),
                      '--bundle-version', 'reviewed-1', '--reason', 'Accepted solver change', '--provenance', 'Review record 42'])
    assert result == 0
    with pytest.raises(BundleError, match='output already exists'):
        publish(setup)
    assert manifest_path.read_bytes() == original
    assert not any(item.is_symlink() for item in setup.output.rglob('*'))
    manifest = json.loads((setup.output / 'manifest.json').read_text())
    assert manifest['bundle_class'] == 'golden'
    reference = setup.output / manifest['artifacts'][manifest['roles']['warm_reference']]['path']
    assert reference.read_bytes() == restart.read_bytes()
    receipt = json.loads((setup.output / 'provenance/publication.json').read_text())
    assert receipt['reason'] == 'Accepted solver change' and receipt['provenance'] == 'Review record 42'
    shutil.rmtree(setup.workspace)
    shutil.rmtree(setup.harness.bundle)
    assert validate_bundle_root(setup.output, setup.root / 'cases').case_id == 'legacy_case'


@pytest.mark.parametrize('failure', ['exit', 'nonconvergence', 'nonfinite', 'validation'])
def test_failed_producer_or_validation_blocks_publication(setup, failure):
    solver = SOLVER
    if failure == 'exit':
        solver = '#!/usr/bin/env bash\nexit 7\n'
    elif failure == 'nonconvergence':
        solver = SOLVER.replace('1.0E-8', '1.0E5')
    elif failure == 'validation':
        solver = SOLVER.replace('set -euo pipefail', 'set -euo pipefail\n[[ "$PWD" != *-verify ]] || exit 7')
    else:
        with h5py.File(setup.harness.serial_executable.parent / 'seed.h5', 'r+') as handle:
            handle['solution/u'][0] = float('nan')
    setup.harness.install_solver(solver, 'serial')
    build = json.loads(setup.build_path.read_text())
    build['artifacts']['serial'].update(file_identity(setup.harness.serial_executable))
    setup.build_path.write_text(json.dumps(build))
    with pytest.raises(BundleError, match='producer failed|producer did not converge|invalid producer output|golden validation failed'):
        golden.refresh('legacy_case', setup.values, setup.workspace)
    report = json.loads((setup.workspace / 'refresh.json').read_text())
    assert report['status'] == 'failed'
    if failure == 'validation':
        assert all(item['status'] == 'passed' for item in report['producers'])
        assert report['checks'][0]['status'] == 'failed'
    else:
        assert len(report['producers']) == 1
        assert not (setup.workspace / 'candidate').exists()
    with pytest.raises(BundleError, match='successful, validated refresh'):
        publish(setup)
    assert not setup.output.exists()


def test_publish_rejects_changed_outputs_evidence_and_missing_reason(setup):
    _, report = golden.refresh('legacy_case', setup.values, setup.workspace)
    with pytest.raises(BundleError, match='reason and provenance'):
        golden.publish(setup.workspace, setup.output, 'v1', ' ', 'Review')
    candidate = setup.workspace / 'candidate'
    manifest_path = candidate / 'manifest.json'
    original = manifest_path.read_bytes()
    manifest_path.write_bytes(original + b'\n')
    with pytest.raises(BundleError, match='evidence changed'):
        publish(setup)
    manifest_path.write_bytes(original)
    run = Path(report['producers'][-1]['run_directory'])
    (run / 'outputs/result.h5').write_bytes(b'changed')
    with pytest.raises(BundleError, match='outputs changed'):
        publish(setup)
    assert not setup.output.exists()
