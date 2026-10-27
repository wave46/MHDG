"""Compile production stabilization and face assemblers for a focused Jacobian check.

Requires existing solver module files/objects in lib (or --build-dir).
The numerical assembly statements are extracted without modification;
inactive diagnostic/pinch branches are omitted to avoid unrelated linkage.
"""
from pathlib import Path
import argparse
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FLAGS = ['-cpp', '-DNGAMMA', '-DTEMPERATURE', '-DNEUTRAL', '-DNEUTRALP',
         '-DAMJUELSPLINES', '-DTHREEBODYREC', '-DEXPANDEDCX', '-DTHERMALCX',
         '-fdefault-real-8', '-fdefault-double-8', '-ffree-line-length-none',
         '-ffunction-sections', '-fdata-sections', '-O0', '-g', '-fcheck=all',
         '-ffpe-trap=invalid,zero,overflow']


def extract(source, name):
    return re.search(r'^\s*SUBROUTINE '+name+r'\(.*?^\s*END\s*SUBROUTINE '+name+r'\b[^\n]*',
                     source, re.I | re.M | re.S)[0]


def omit_inactive(source):
    result, depth = [], 0
    for line in source.splitlines():
        code = line.strip().lower()
        if not depth and any('if ('+name+') then' in code for name in
                             ['switch%transport_1d', 'diagnostics_on', 'face_diagnostics_on']):
            depth = 1
            continue
        if depth:
            if re.match(r'if\s*\(.*\)\s*then', code):
                depth += 1
            if re.match(r'end\s*if\b', code):
                depth -= 1
            continue
        result.append(line)
    assert depth == 0
    return '\n'.join(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, default=ROOT/'lib')
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--constant-neutral', action='store_true')
    args = parser.parse_args()
    if args.optimized:
        FLAGS[:] = [x for x in FLAGS if x not in ['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']] + ['-Ofast','-march=native']
    if args.constant_neutral:
        FLAGS.append('-DCONSTANTNEUTRALDIFF')
    build = args.build_dir.resolve()
    with tempfile.TemporaryDirectory(prefix='stabilization_') as directory:
        tmp = Path(directory)
        def run(command, **kwargs):
            return subprocess.run(command, cwd=tmp, check=True, **kwargs)
        preprocess = lambda path: subprocess.check_output(
            ['gfortran', *FLAGS, '-E', '-P', str(path)], text=True)
        jac = preprocess(ROOT/'src/HDG/hdg_ComputeJacobian.f90')
        bc = preprocess(ROOT/'src/HDG/hdg_BC.f90')
        routines = '\n'.join(omit_inactive(extract(source, name)) for source, name in [
            (jac, 'assemblyIntFacesContribution'), (jac, 'assemblyExtFacesContribution'),
            (bc, 'assembly_bohm_bc')])
        wrapper = '''module extracted_faces
use globals
use physics
use smooth_transition
use bohm_boundary
use hdg_stabilization
implicit none
integer,parameter::neq=5,ndim=2,npfl=1
integer::fl=1
real*8::coefi,coefe
logical::face_diagnostics_on=.false.
type balance_accumulator_type
integer::unused
end type
type balance_plasma_face_type
integer::unused
end type
type balance_neutral_face_type
integer::unused
end type
type balance_plasma_bc_type
integer::unused
end type
contains
'''+routines+'\nend module\n'
        (tmp/'faces.f90').write_text(wrapper)
        inc = ['-I', str(tmp), '-I', str(build), '-J', str(tmp)]
        objects = []
        for path in [ROOT/'src/Models/NGammaTiTe/physics.f90',
                     ROOT/'src/Models/NGammaTiTe/physics_stabilization.f90',
                     ROOT/'src/Models/bohm_boundary.f90', ROOT/'src/HDG/hdg_Stabilization.f90', tmp/'faces.f90']:
            obj = tmp/(path.stem+'.o')
            run(['gfortran', *FLAGS, *inc, '-c', str(path), '-o', str(obj)])
            objects.append(str(obj))
        run(['gfortran', *FLAGS, *inc, str(ROOT/'unit_tests/fortran/test_stabilization.f90'),
             *objects, *[str(build/name) for name in ['globals.o', 'neutral_flux_limiter.o', 'smooth_transition.o', 'interpolation.o']],
             '-Wl,--gc-sections', '-o', str(tmp/'probe')])
        result = subprocess.run([str(tmp/'probe')], cwd=tmp, capture_output=True, text=True)
        print(result.stdout, end='')
        if result.returncode:
            print(result.stderr, end='')
            result.check_returncode()
        if 'Stabilization PASS' not in result.stdout:
            raise RuntimeError('Stabilization fixture did not reach its final check')


if __name__ == '__main__':
    main()
