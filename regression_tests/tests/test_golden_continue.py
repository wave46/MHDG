"""Continue real tiny campaigns without regenerating successful producers."""

import json
from pathlib import Path
import shutil
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from regression_tests import cli, golden, golden_continue, suites
from regression_tests.files import file_identity
from regression_tests.support import BundleError
from regression_tests.tests.test_golden import setup, SOLVER


@pytest.fixture
def failed_campaign(setup, monkeypatch):
    path = setup.root / 'cases/legacy_case.json'
    case = json.loads(path.read_text())
    case['workflows']['repeat'] = {'extends': 'warm', 'outputs': []}
    path.write_text(json.dumps(case))
    path = setup.root / 'suites.json'
    selection = json.loads(path.read_text())
    selection['suites']['warm']['workflows'] = ['warm', 'repeat']
    path.write_text(json.dumps(selection))
    setup.harness.install_solver(SOLVER + '''
if [[ "$PWD" == */repeat/* ]]; then
  printf 'Error opening source file: fixture failure\n'
fi
''', 'serial')
    with pytest.raises(BundleError, match='golden validation failed'):
        golden.refresh('legacy_case', setup.values, setup.workspace)
    # Keep old executable/build files intact; use a separate corrected build.
    directory = setup.workspace.parent / 'corrected-build'
    shutil.copytree(setup.harness.build_manifest.parent, directory)
    (directory / setup.harness.serial_executable.name).write_text(SOLVER)
    path = directory / setup.harness.build_manifest.name
    build = json.loads(path.read_text())
    record = build['artifacts']['NGammaTiTeNeutral/serial']
    record.update(file_identity(directory / record['path']))
    path.write_text(json.dumps(build))
    setup.continuation_build = Mock(return_value=SimpleNamespace(metadata_path=path))
    monkeypatch.setattr(golden_continue, 'build_solver', setup.continuation_build)
    return setup


def test_continue_reuses_producers_and_passed_cells_and_publishes_provenance(failed_campaign, monkeypatch):
    setup = failed_campaign
    original = json.loads((setup.workspace/'refresh.json').read_text())
    producer = Path(original['producers'][0]['run_directory'])
    old_metadata = file_identity(producer/'run_metadata.json')
    old_build = file_identity(setup.workspace/'build.json')
    run = Mock(wraps=suites.run_cell)
    monkeypatch.setattr(suites, 'run_cell', run)
    assert cli.main(['golden', 'continue', str(setup.workspace), '--settings', str(setup.harness.settings),
                     '--rerun', 'repeat']) == 0
    report = json.loads((setup.workspace/'refresh.json').read_text())
    assert report['status'] == 'ready'
    assert [call.args[1] for call in run.call_args_list] == ['repeat']
    assert file_identity(producer/'run_metadata.json') == old_metadata
    assert file_identity(setup.workspace/'build.json') == old_build
    summary = json.loads((setup.workspace/report['checks'][0]['summary']).read_text())
    assert summary['results'][0]['reused_execution']
    assert summary['results'][1]['run_directory'] != original['producers'][0]['run_directory']
    assert report['continuations'][0]['reused_producers'] == ['warm']
    golden.publish(setup.workspace, setup.output, 'continued', 'Reviewed continuation', 'Fixture review', catalog_root=setup.root)
    receipt = json.loads((setup.output/'provenance/publication.json').read_text())
    assert receipt['build']['record'] == old_build
    assert receipt['continuations'][0]['validation_build']['artifacts']['NGammaTiTeNeutral/serial'] == file_identity(
        setup.workspace.parent/'corrected-build'/setup.harness.serial_executable.name)


@pytest.mark.parametrize('changed', ['producer', 'candidate', 'recipe'])
def test_continue_rejects_changed_evidence_before_building(failed_campaign, changed):
    setup = failed_campaign
    report = json.loads((setup.workspace/'refresh.json').read_text())
    if changed == 'producer':
        path = Path(report['producers'][0]['run_directory'])/'outputs/result.h5'
    elif changed == 'candidate':
        path = setup.workspace/'candidate/inputs/reference_mpi4_omp4.h5'
    else:
        path = setup.root/'golden.json'
    with path.open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(BundleError):
        golden_continue.continue_refresh(setup.workspace, setup.values, ['repeat'], catalog_root=setup.root)
    setup.continuation_build.assert_not_called()
