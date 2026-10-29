"""Build an observational test binary without editing production sources."""
from pathlib import Path
import shlex
import shutil
import subprocess

from .files import file_identity
from .support import BundleError


def build_probe(repository, destination, environment, make_arguments, log_path):
    library = repository / 'lib'
    probe = destination.parent.parent / 'newton_probe' / destination.name
    # Reuse the configured Makefile, objects and compiler modules together.
    # Compiling in this copy also keeps compiler-generated modules private.
    work = probe / 'lib'
    shutil.copytree(library, work)
    (probe / 'test').mkdir()
    source = (repository / 'src/Utils/Main_utils.f90').read_text()
    replacements = {
        '\nMODULE Main_utils\n': '\nMODULE Main_utils\n  USE newton_linear_probe\n',
        '    previous_u=sol%u;previous_q=sol%q;previous_trace=sol%u_tilde\n':
            '    CALL capture_newton_linear_residual()\n'
            '    previous_u=sol%u;previous_q=sol%q;previous_trace=sol%u_tilde\n',
        '    CALL finish_newton_update(status, previous_u, previous_q, previous_trace, numer%dumpnr)\n':
            '    IF(status%outcome==solve_success) CALL report_newton_linear_residual()\n'
            '    CALL finish_newton_update(status, previous_u, previous_q, previous_trace, numer%dumpnr)\n',
    }
    for old, new in replacements.items():
        if source.count(old) != 1:
            raise BundleError('Newton probe insertion point changed; review the test instrumentation')
        source = source.replace(old, new)
    observed = work / 'Main_utils.F'
    observed.write_text(source)
    fixture = Path(__file__).parent / 'probes/newton_linear_probe.f90'
    rules = fixture.with_name('newton_probe.mk')
    shutil.copy2(fixture, work / fixture.name)
    shutil.copy2(rules, work / rules.name)
    # Keep the observed generated source rather than recopying its original.
    # Removing its copied object forces compilation even with make's -o flag.
    (work / 'Main_utils.o').unlink()
    command = ['make', '-f', 'Makefile', '-f', rules.name, '-o', observed.name,
               f'SDIR={repository / "src"}', *make_arguments]
    with log_path.open('w') as log:
        log.write(shlex.join(command)+'\n'); log.flush()
        try:
            subprocess.run(command, cwd=work, env={**environment, 'PWD': str(work)},
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            raise BundleError(f'Newton observation build failed; see {log_path}') from exc
    shutil.copy2(work / destination.name, destination)
    return {'kind': 'mixed_linear_residual_observation',
            'command': command,
            'fixture': {'path': str(fixture.resolve()), **file_identity(fixture)},
            'rules': {'path': str(rules.resolve()), **file_identity(rules)},
            'main_utils': {'path': str(observed), **file_identity(observed)}}
