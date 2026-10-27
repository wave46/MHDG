"""Check production thermal blocks in volume, face and prescribed-trace assembly.

The extracted numerical statements are unmodified. Other physical terms in
these blocks are zero in the fixture, so the test isolates parallel conduction.
Requires the current temperature-model build and NumPy.
"""
from pathlib import Path
import argparse
import ctypes
from itertools import product
import re
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FLAGS = ['-cpp', '-DNGAMMA', '-DTEMPERATURE', '-DNEUTRAL', '-DNEUTRALP',
         '-fdefault-real-8', '-fdefault-double-8', '-ffree-line-length-none',
         '-ffunction-sections', '-fdata-sections', '-O0', '-g', '-fcheck=all']
B = np.array([.2, .4])
NORMAL = np.array([.6, .8])
MREF = 1.3
CAPS = np.array([.6, .9])
PARALLEL = np.array([.4, .7])
FRACTIONS = np.array([.2, .3])


def routine(source, name):
    return re.search(r'^\s*SUBROUTINE '+name+r'\(.*?^\s*END\s*SUBROUTINE '+name+r'\b[^\n]*',
                     source, re.M | re.S | re.I)[0]


def production_wrapper(flags=FLAGS):
    def preprocess(path):
        return subprocess.check_output(['gfortran', *flags, '-E', '-P', str(path)], text=True)
    jacobian = preprocess(ROOT/'src/HDG/hdg_ComputeJacobian.f90')
    boundary = preprocess(ROOT/'src/HDG/hdg_BC.f90')
    sources = [('VOLUME', jacobian, 'assemblyVolumeContribution'),
               ('INTERIOR', jacobian, 'assemblyIntFacesContribution'),
               ('EXTERIOR', jacobian, 'assemblyExtFacesContribution'),
               ('BOHM', boundary, 'assembly_bohm_bc'),
               ('PRESCRIBED', boundary, 'assembly_dirichletstr_bc')]
    wrapper = (ROOT/'unit_tests/fortran/test_heat_flux_assembly.f90').read_text()
    for label, source, name in sources:
        body = routine(source, name)
        call = re.search(r'CALL compute_spitzer_harm_flux_ions\(.*?CALL compute_heat_flux_limiters\(.*?\)', body, re.S | re.M)[0]
        if label == 'BOHM':
            block = body.split('! Parallel diffusion for temperature', 1)[1].split('!endif', 1)[0]
        elif label == 'PRESCRIBED':
            block = body.split('! assemblyExtFacesContribution supplies the constant heat-flux RHS term.', 1)[1]
            block = block[:block.lower().index('end subroutine')]
        else:
            block = body.split('! Parallel diffusion for the temperature', 1)[1]
            block = re.split(r'ELSEIF\s*\(i == inn\)', block, maxsplit=1)[0] + 'ENDIF\n'
        wrapper = wrapper.replace('! PRODUCTION_'+label, call+'\n'+block)
    return wrapper


def heat(state, gradient):
    density, velocity = state[0], state[1]/state[0]
    specific = state[2:4]/density - np.array([.5*velocity**2, 0.])
    alpha = np.minimum(np.maximum(specific, 1.e-20), 1.5*MREF*CAPS)**2.5
    vectors = np.zeros((5, 2))
    vectors[0] = [-state[2]/density**2+state[1]**2/density**3, -state[3]/density**2]
    vectors[1, 0] = -state[1]/density**2
    vectors[2, 0] = vectors[3, 1] = 1./density
    return PARALLEL*(2./(3.*MREF))**3.5*alpha*(B @ gradient @ vectors)


def limiter(state, gradient, enabled):
    if not enabled:
        return np.ones(2)
    density, velocity = state[0], state[1]/state[0]
    specific = state[2:4]/density-np.array([.5*velocity**2, 0.])
    floor = 1.e-5*1.5*MREF
    fs = density*(2./3.*np.maximum(specific, floor))**1.5
    fs[1] *= np.sqrt(3.35e-27/9.109e-31)
    return 1./(1.+abs(heat(state, gradient))/(FRACTIONS*np.maximum(fs, 1.e-10)))


def reference(route, state, gradient, phi):
    result = np.zeros(5)
    factor = 1. if route == 1 else -(B @ NORMAL)
    result[2:4] = factor*phi*heat(state, gradient)
    return result


def check_fixture(probe, state, gradient, route, enabled):
    array = lambda x: np.asfortranarray(x, dtype=np.float64)
    ju, jq, rhs, actual_phi = array(np.zeros((5, 5))), array(np.zeros((5, 10))), np.zeros(5), np.zeros(2)
    probe(array(state), array(gradient), route, enabled, ju, jq, rhs, actual_phi)
    phi = limiter(state, gradient, enabled)
    np.testing.assert_allclose(actual_phi, phi, rtol=2.e-13, atol=1.e-15)
    expected = reference(route, state, gradient, phi)
    if route == 5:
        # Known trace contributes its actual Jacobian column times U_hat.
        np.testing.assert_allclose(rhs, -expected, rtol=2.e-12, atol=2.e-12)
        return 0
    reconstructed = ju @ state + jq @ gradient.ravel(order='F')-rhs
    np.testing.assert_allclose(reconstructed, expected, rtol=3.e-12, atol=2.e-12)
    jacobian = np.column_stack((ju, jq))
    x = np.r_[state, gradient.ravel(order='F')]
    for j in range(len(x)):
        step = 1.e-7*max(abs(x[j]), .1*state[0])
        plus, minus = x.copy(), x.copy()
        plus[j] += step
        minus[j] -= step
        def value(z):
            return reference(route, z[:5], z[5:].reshape((2, 5), order='F'), phi)
        numeric = (value(plus)-value(minus))/(2.*step)
        np.testing.assert_allclose(numeric, jacobian[:, j], rtol=2.e-6, atol=2.e-7)
    return len(x)


def run_checks(build, optimized):
    flags = FLAGS
    if optimized:
        flags = [flag for flag in FLAGS if flag not in ['-O0', '-g', '-fcheck=all']]
        flags += ['-Ofast', '-march=native']
    with tempfile.TemporaryDirectory(prefix='heat_assembly_') as directory:
        tmp = Path(directory)
        (tmp/'extracted.f90').write_text(production_wrapper(flags))
        objects = []
        for path in [ROOT/'src/Definitions/globals.f90',
                     ROOT/'src/Models/NGammaTiTe/physics.f90', tmp/'extracted.f90']:
            obj = tmp/(path.stem+'.o')
            subprocess.run(['gfortran', *flags, '-fPIC', '-c', '-I', str(tmp), '-I', str(build),
                            '-J', str(tmp), str(path), '-o', str(obj)], cwd=tmp, check=True)
            objects.append(str(obj))
        # Export only the probe; discard unrelated physics and solver sections.
        (tmp/'exports.map').write_text('{ global: probe; local: *; };')
        subprocess.run(['gfortran', '-shared', *objects, '-Wl,--gc-sections',
                        '-Wl,--version-script='+str(tmp/'exports.map'),
                        '-o', str(tmp/'probe.so')], cwd=tmp, check=True)
        library = ctypes.CDLL(str(tmp/'probe.so'))
        array = np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS')
        library.probe.argtypes = [array, array, ctypes.c_int, ctypes.c_int, array, array, array, array]
        library.probe.restype = None
        states = [np.array([.8, -.15, 1.6, 2.2, .1]),
                  np.array([1., .2, .5, .8, .1]),
                  np.array([1., 0., -1.e-6, 1., .1])]
        gradient = np.linspace(-.25, .4, 10).reshape((2, 5), order='F')
        columns, cases = 0, 0
        for state, sign, route, enabled in product(states, [-1, 0, 1], range(1, 6), [False, True]):
            try:
                columns += check_fixture(library.probe, state, sign*gradient, route, enabled)
            except AssertionError as error:
                raise AssertionError(f'route={route}, enabled={enabled}, sign={sign}, state={state}') from error
            cases += 1
        label = 'optimized' if optimized else 'bounds-checked'
        print(f'Heat assembly PASS ({label}): {cases} cases, {columns} held-phi derivative columns; all five paths')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, default=ROOT/'lib')
    args = parser.parse_args()
    for optimized in [False, True]:
        run_checks(args.build_dir.resolve(), optimized)


if __name__ == '__main__':
    main()
