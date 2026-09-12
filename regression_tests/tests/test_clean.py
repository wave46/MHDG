"""Protect real dependency shapes without solver data or catalog snapshots."""

import json
import pytest

from regression_tests import cli
from regression_tests.clean import cleanup
from support.errors import BundleError, DocumentError


def record(directory, name, **values):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(json.dumps(values))
    return path


def run(root, name, **values):
    directory = root / name
    record(directory, 'run_plan.json', workflow_id='warm', **values)
    record(directory, 'run_metadata.json', status='completed')
    return directory


def test_cli_preview_and_explicit_removal_do_not_follow_links(tmp_path, capsys):
    root = tmp_path / 'storage'
    stale = run(root, 'stale')
    outside = tmp_path / 'source'
    outside.mkdir()
    payload = outside / 'input.h5'
    payload.write_bytes(b'immutable data')
    (stale / 'inputs').symlink_to(outside, target_is_directory=True)
    assert cli.main(['clean', str(root), 'stale']) == 0
    output = capsys.readouterr().out
    assert str(stale) in output and 'bytes' in output and 'preview' in output
    assert stale.exists()
    assert cli.main(['clean', str(root), 'stale', '--delete']) == 0
    assert not stale.exists() and payload.read_bytes() == b'immutable data'
    with pytest.raises(BundleError, match='explicit directory'):
        cleanup(root, delete=True)


def test_retained_suite_and_run_protect_dependencies_until_selected_together(tmp_path):
    root = tmp_path / 'storage'
    build = root / 'build'
    exe = record(build, 'build_metadata.json', status='completed')
    producer = run(root, 'producer', command=[str(exe)])
    kept = run(root, 'consumer')
    (kept / 'restart.h5').symlink_to(producer / 'run_metadata.json')
    suite = root / 'suite'
    record(suite, 'suite_summary.json', status='passed', results=[{'run_directory': str(kept)}])
    report = cleanup(root, [producer, build])
    assert report['blocked']
    with pytest.raises(BundleError, match='needed by'):
        cleanup(root, [producer, kept, build], delete=True)
    assert all(p.exists() for p in (producer, kept, build))
    cleanup(root, [suite, kept, producer, build], delete=True)
    assert not any(root.iterdir())


def test_unknown_active_defaults_and_external_consumers_are_protected(tmp_path):
    root = tmp_path / 'storage'
    active = root / 'active'
    record(active, 'refresh.json', status='running')
    unknown = root / 'source'
    unknown.mkdir()
    (unknown / 'mesh.msh').write_text('source')
    default = root / 'golden'
    record(default, 'manifest.json', bundle_class='golden')
    settings = record(tmp_path, 'settings.json', defaults={'bundles': {'case': str(default)}})
    old = run(root, 'old')
    external = run(tmp_path, 'external', artifacts={'restart': str(old / 'result.h5')})
    for selected, options in ((active, {}), (unknown, {}), (default, {'settings': settings}),
                              (old, {'keep': [external]})):
        with pytest.raises(BundleError, match='protected cleanup selection'):
            cleanup(root, [selected], delete=True, **options)
        assert selected.exists()
    # A standalone published golden's copied provenance is not a live dependency.
    record(default / 'provenance', 'refresh.json', status='ready', workspace=str(old))
    cleanup(root, [old], settings=settings, delete=True)
    assert default.exists() and external.exists()
    # Explicit selection of an old golden is allowed after removing its default.
    settings.write_text('{}')
    cleanup(root, [default], settings=settings, delete=True)
    assert not default.exists()


def test_invalid_selection_and_unreadable_records_cannot_remove_data(tmp_path):
    root = tmp_path / 'storage'
    old = run(root, 'old')
    alias = root / 'alias'
    alias.symlink_to(old, target_is_directory=True)
    for selection in (alias, tmp_path, old / 'run_metadata.json', root):
        with pytest.raises(BundleError):
            cleanup(root, [selection], delete=True)
    (old / 'run_metadata.json').write_text('{invalid')
    with pytest.raises(DocumentError, match='invalid JSON'):
        cleanup(root, [old], delete=True)
    assert old.exists()
