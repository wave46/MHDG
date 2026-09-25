"""Exercise Git change discovery and useful recommendations, not table snapshots."""

import subprocess

import pytest

from regression_tests import cli, reporting
from regression_tests.catalog import load_selection
from regression_tests.suggest import ROOT, changed_paths, recommendations
from regression_tests.support import BundleError


def git(root, *arguments):
    return subprocess.run(['git', *arguments], cwd=root, text=True, check=True,
                          capture_output=True).stdout.strip()


def write(root, path, text='initial\n'):
    target = root / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return target


def commit(root):
    git(root, 'add', '.')
    git(root, '-c', 'user.name=Harness test', '-c', 'user.email=harness@example.invalid',
        'commit', '-qm', 'Test fixture')
    return git(root, 'rev-parse', 'HEAD')


def test_git_base_local_layers_renames_and_cli_are_read_only(tmp_path, monkeypatch, capsys):
    git(tmp_path, 'init', '-q')
    old = 'src/MPI_OMP/old name.f90'
    cancelled = 'src/Adaptivity/adaptivity_projection.f90'
    write(tmp_path, old)
    write(tmp_path, cancelled)
    write(tmp_path, '.gitignore', '*.log\n')
    base = commit(tmp_path)
    committed = 'src/Models/NGammaTiTe/transport_1d/transport_models_1d.f90'
    write(tmp_path, committed)
    commit(tmp_path)
    renamed = 'src/MPI_OMP/new name.f90'
    git(tmp_path, 'mv', old, renamed)
    write(tmp_path, cancelled, 'staged\n')
    git(tmp_path, 'add', cancelled)
    write(tmp_path, cancelled)  # Index/worktree changes cancel relative to HEAD.
    untracked = 'src/unmapped.f90'
    write(tmp_path, untracked)
    write(tmp_path, 'ignored.log')
    assert set(changed_paths(tmp_path, base)) == {committed, old, renamed, cancelled, untracked}
    assert committed not in changed_paths(tmp_path)
    before = git(tmp_path, 'status', '--porcelain')
    monkeypatch.setattr(cli, 'ROOT', tmp_path / 'regression_tests')
    assert cli.main(['suggest', '--base', base]) == 0
    output = capsys.readouterr().out
    assert f'unmapped: {untracked}' in output and 'nothing executed' in output
    assert git(tmp_path, 'status', '--porcelain') == before
    assert git(tmp_path, 'rev-parse', 'HEAD') != base
    with pytest.raises(BundleError, match='cannot inspect Git'):
        changed_paths(tmp_path, '--not-a-revision')


def test_focused_recommendations_select_the_claimed_scientific_checks():
    paths = ['src/Models/NGammaTiTe/transport_1d/transport_models_1d.f90',
             'src/MPI_OMP/Communications.f90']
    report = recommendations(paths)
    selections = {}
    for check in report['checks']:
        command = check['command']
        if command[:4] != ['python', '-m', 'regression_tests', 'check']:
            continue
        name = command[4]
        _, suites = load_selection(name, ROOT / 'suites.json', ROOT / 'layouts.json', ROOT / 'cases')
        selections[name] = suites[0]
    assert 'transport' in selections['transport']['workflow_ids']
    assert 'transport_short' in selections['transport_hybrid']['workflow_ids']
    assert selections['transport_hybrid']['layout_comparisons']
    assert selections['adaptive_parallel']['layout_comparisons']
    assert all('--build' in check['command'] for check in report['checks'])

    build = recommendations(['lib/Makefile'])['checks']
    assert len(build) == 1 and build[0]['command'][-2:] == ['full', '--build']


def test_broader_profile_absorbs_feature_checks_and_keeps_explanations():
    paths = ['src/InOut/read_input.f90', 'src/MPI_OMP/Communications.f90',
             'src/Models/neutral_flux_limiter.f90', 'src/Utils/Diagnostics/balance_diagnostics_output.f90',
             'regression_tests/compare.py', 'src/Models/Laplace/physics.f90', 'README.md']
    report = recommendations(paths)
    scientific = [check for check in report['checks'] if 'check' in check['command']]
    assert len(scientific) == 1 and scientific[0]['command'][-2:] == ['full', '--build']
    assert any(check['command'][2] == 'pytest' for check in report['checks'])
    assert all(path in scientific[0]['reasons'] for path in paths[:4])
    assert report['unmapped'] == ['src/Models/Laplace/physics.f90']
    assert report['documentation'] == ['README.md']


def test_harness_scope_selects_real_evidence_without_rebuilding_and_groups_output(capsys):
    paths = ['regression_tests/prepare.py', 'regression_tests/compare.py', 'regression_tests/tests/test_golden.py']
    report = recommendations(paths)
    commands = [check['command'] for check in report['checks']]
    assert ['python', '-m', 'regression_tests', 'check', 'full'] in commands
    assert any(command[2] == 'pytest' for command in commands)
    assert all('--build' not in command for command in commands)
    reporting.suggestions({'paths': paths, 'base': None, **report})
    output = capsys.readouterr().out
    assert output.count('python -m regression_tests check full') == 1
    assert all(path in output for path in paths)
    report = recommendations(['regression_tests/suggest.py', 'regression_tests/clean.py'])
    assert len(report['checks']) == 1 and report['checks'][0]['command'][2] == 'pytest'
    report = recommendations(['regression_tests/README.md'])
    assert not report['checks'] and not report['unmapped']
