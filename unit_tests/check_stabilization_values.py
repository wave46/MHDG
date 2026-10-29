"""Check production automatic tau values and geometry in 2D and TOR3D.

Uses existing support modules in --build-dir; each physics variant is compiled
in isolation, without rebuilding or running the solver.
"""
from pathlib import Path
import argparse
import subprocess
import tempfile

from check_stabilization import FLAGS, ROOT


def check_dimension(build, ndim):
    flags = [*FLAGS, *(['-DTOR3D'] if ndim == 3 else [])]
    with tempfile.TemporaryDirectory(prefix=f'tau_values_{ndim}d_') as directory:
        tmp = Path(directory)
        includes = ['-I', str(tmp), '-I', str(build), '-J', str(tmp)]

        def run(command):
            subprocess.run(command, cwd=tmp, check=True)

        objects = []
        for name in ['physics', 'physics_stabilization']:
            source = ROOT/'src/Models/NGammaTiTe'/f'{name}.f90'
            obj = tmp/f'{name}.o'
            run(['gfortran', *flags, *includes, '-c', str(source), '-o', str(obj)])
            objects.append(str(obj))
        run(['gfortran', *flags, *includes,
             str(ROOT/'unit_tests/fortran/test_stabilization_values.f90'),
             *objects, str(build/'globals.o'), str(build/'neutral_flux_limiter.o'),
             '-Wl,--gc-sections', '-o', str(tmp/'probe')])
        result = subprocess.run([str(tmp/'probe')], cwd=tmp, check=True,
                                capture_output=True, text=True)
        print(result.stdout, end='')
        if 'Stabilization values: PASS' not in result.stdout:
            raise RuntimeError(f'{ndim}D coefficient check did not finish')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, default=ROOT/'lib')
    args = parser.parse_args()
    for ndim in (2, 3):
        check_dimension(args.build_dir.resolve(), ndim)


if __name__ == '__main__':
    main()
