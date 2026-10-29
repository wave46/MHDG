"""Exercise the observer through make with alternative compiler commands."""
import json
import os
from pathlib import Path
import shlex
import sys

import pytest

from regression_tests.files import file_identity
from regression_tests.newton_probe import build_probe
from regression_tests.support import BundleError


@pytest.fixture
def configured_build(tmp_path):
    repository = tmp_path / 'repository'
    library = repository / 'lib'
    (library / 'Make.inc').mkdir(parents=True)
    (repository / 'test').mkdir()
    source = repository / 'src/Utils/Main_utils.f90'
    source.parent.mkdir(parents=True)
    source.write_text((Path(__file__).resolve().parents[2] / 'src/Utils/Main_utils.f90').read_text())
    compiler = tmp_path / 'alternative-fc'
    compiler.write_text('#!' + sys.executable + '\n' + COMPILER)
    compiler.chmod(0o755)
    (library / 'Make.inc/arch.make').write_text(
        f'FC = {compiler}\nFCFLAGS = --configured-flag\nLIB = --configured-library\n')
    (library / 'Makefile').write_text(MAKEFILE)
    for name in ('Main_utils', 'hdg_Residual', 'MPI_OMP', 'driver'):
        (library / (name + '.o')).write_text('ordinary-' + name)
        (library / (name.lower() + '.mod')).write_text('ordinary-module')
    return repository, compiler


def library_identity(repository):
    return {str(path.relative_to(repository)): file_identity(path)
            for path in (repository / 'lib').rglob('*') if path.is_file()}


@pytest.mark.parametrize('wrapped', [False, True])
def test_make_owns_compiler_flags_linking_and_module_isolation(configured_build, tmp_path, wrapped):
    repository, compiler = configured_build
    before = library_identity(repository)
    compiler_command = shlex.quote(str(compiler))
    if wrapped:
        compiler_command = 'env PROBE_WRAPPER=enabled ' + compiler_command + ' --wrapped'
    commands = []
    for variant in ('serial', 'parall'):
        destination = tmp_path / 'build/bin' / ('MHDG-' + variant)
        destination.parent.mkdir(parents=True, exist_ok=True)
        record = build_probe(repository, destination, os.environ,
                             [f'FC={compiler_command}', f'MODE={variant}', destination.name],
                             tmp_path / (variant + '.log'))
        assert record['kind'] == 'mixed_linear_residual_observation'
        assert destination.read_text() == 'observed driver|ordinary-driver|observer'
        work = Path(record['main_utils']['path']).parent
        calls = [json.loads(line) for line in (work / 'compiler.jsonl').read_text().splitlines()]
        assert all('--configured-flag' in call['arguments'] for call in calls)
        assert calls[-1]['arguments'][-2:] == ['--configured-library', 'newton_linear_probe.o']
        assert all(call['wrapper'] == ('enabled' if wrapped else None) for call in calls)
        assert (work / 'main_utils.mod').read_text() == 'private-module'
        commands.append(record['command'])
    assert commands[0] != commands[1]
    assert library_identity(repository) == before


def test_failed_observer_build_keeps_normal_build_and_log(configured_build, tmp_path):
    repository, _ = configured_build
    before = library_identity(repository)
    destination = tmp_path / 'build/bin/MHDG-serial'
    destination.parent.mkdir(parents=True)
    destination.write_text('normal executable')
    log = tmp_path / 'build.log'
    with pytest.raises(BundleError, match='observation build failed; see'):
        build_probe(repository, destination, {**os.environ, 'FAIL_PROBE': '1'},
                    ['MODE=serial', destination.name], log)
    assert 'deliberate compiler failure' in log.read_text()
    assert destination.read_text() == 'normal executable'
    assert library_identity(repository) == before


MAKEFILE = '''include Make.inc/arch.make
Main_utils.F: $(SDIR)/Utils/Main_utils.f90
\tcp $< $@
Main_utils.o: Main_utils.F hdg_Residual.o MPI_OMP.o
\t$(FC) -c $(FCFLAGS) $< -o $@
MHDG-$(MODE): Main_utils.o driver.o
\t$(FC) $(FCFLAGS) Main_utils.o driver.o -o $@ $(LIB)
\tcp $@ ../test/
'''

COMPILER = '''import json, os, pathlib, sys
arguments = sys.argv[1:]
with open('compiler.jsonl', 'a') as log:
    log.write(json.dumps({'arguments': arguments, 'wrapper': os.getenv('PROBE_WRAPPER')})+'\\n')
if os.getenv('FAIL_PROBE'):
    sys.exit('deliberate compiler failure')
output = pathlib.Path(arguments[arguments.index('-o')+1])
if '-c' in arguments:
    source = pathlib.Path(next(x for x in arguments if x.endswith(('.F', '.f90'))))
    if source.name == 'Main_utils.F':
        assert 'USE newton_linear_probe' in source.read_text()
        assert pathlib.Path('newton_linear_probe.mod').exists()
        output.write_text('observed driver')
    else:
        output.write_text('observer')
    pathlib.Path(source.stem.lower()+'.mod').write_text('private-module')
else:
    output.write_text('|'.join(pathlib.Path(x).read_text() for x in arguments if x.endswith('.o')))
'''
