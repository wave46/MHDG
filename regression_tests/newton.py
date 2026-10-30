"""Compare complete first updates and observed mixed linear residuals."""
from contextlib import ExitStack
from pathlib import Path
import math
import re

import h5py
import numpy as np

from .compare_common import calculate_error_norms, load_fixed_tolerances
from .compare_fixed import compare_mesh, read_solution_field, required_scalar
from .documents import load_json
from .files import file_identity
from .support import ComparisonError

LINEAR = re.compile(r'Newton linear check: initial=\s*(\S+); trial=\s*(\S+); relative=\s*(\S+)')
RESIDUAL = re.compile(r'HDG residual before update: norm=\s*([^;\s]+)')
ERROR = re.compile(r'^\s*Error:\s*(\S+)', re.M)


def _snapshot(run, pattern):
    paths = sorted((run/'outputs').glob(pattern))
    if len(paths) != 1:
        raise ComparisonError(f'expected one {pattern} snapshot in {run}; found {len(paths)}')
    return paths[0]


def _history(run):
    text = (run/'stdout.log').read_text()
    return {'nonlinear_residuals': [float(x) for x in RESIDUAL.findall(text)],
            'update_errors': [float(x.replace('D', 'E')) for x in ERROR.findall(text)],
            'linear': [dict(zip(('initial', 'trial', 'relative'), map(float, row)))
                       for row in LINEAR.findall(text)]}


def compare_newton_runs(baseline, candidate, tolerances_path, profile, *, catalog=None):
    workflow = {'comparison': {'profile': profile}, 'layout': 'same'}
    _, tolerances = load_fixed_tolerances(tolerances_path, workflow, 'same', None, catalog=catalog)
    failures, report = [], {}
    try:
        plans = [load_json(run/'run_plan.json', 'Newton run plan') for run in (baseline, candidate)]
        controls = [{key: value for key, value in plan['parameter_overrides'].items()
                     if key.lower() != 'newton_formulation'} for plan in plans]
        if controls[0] != controls[1] or plans[0]['layout'] != plans[1]['layout']:
            failures.append('Newton pair has different numerical/physics controls or layout')
        for name in ('restart.h5', 'geometry.geo', 'equilibrium.h5'):
            if file_identity(baseline/'inputs'/name) != file_identity(candidate/'inputs'/name):
                failures.append(f'Newton pair has different input: {name}')
        report['updates'] = _compare_first_updates(baseline, candidate, tolerances, failures)
        report['histories'] = histories = {name: _history(run) for name, run in
                                         [('absolute', baseline), ('incremental', candidate)]}
        for name, history in histories.items():
            numbers = [*history['nonlinear_residuals'], *history['update_errors'],
                       *(x for row in history['linear'] for x in row.values())]
            if not numbers or not all(math.isfinite(x) and x >= 0 for x in numbers):
                failures.append(f'{name}: missing or nonfinite Newton observations')
            if not history['update_errors']:
                failures.append(f'{name}: Newton update monitoring is missing')
            if not history['nonlinear_residuals']:
                failures.append(f'{name}: nonlinear residual monitoring is missing')
            if not history['linear']:
                failures.append(f'{name}: linear check missing; build this suite with --build')
            elif history['linear'][0]['relative'] > tolerances['linear_residual_relative_max']:
                failures.append(f'{name}: first mixed linear equation mismatch exceeds tolerance')
    except (OSError, ValueError, KeyError, ComparisonError) as exc:
        failures.append(str(exc))
    return {**report, 'status': 'failed' if failures else 'passed', 'failures': failures}


def _compare_first_updates(baseline, candidate, tolerances, failures):
    initial = [_snapshot(run, '*_0000.h5') for run in (baseline, candidate)]
    updated = [_snapshot(run, '*_NR0001.h5') for run in (baseline, candidate)]
    report = {'initial': list(map(str, initial)), 'first_updated': list(map(str, updated))}
    with ExitStack() as stack:
        old = [stack.enter_context(h5py.File(path)) for path in initial]
        new = [stack.enter_context(h5py.File(path)) for path in updated]
        for handle in [old[1], *new]:
            compare_mesh(old[0], handle, tolerances, failures)
        equations = int(required_scalar(old[0], 'simulation_parameters', 'Neq'))
        for handle, expected in zip(new, ('absolute', 'incremental')):
            value = np.asarray(handle['simulation_parameters/numerics/newton_formulation'][()]).reshape(-1)[0]
            if isinstance(value, bytes):
                value = value.decode()
            if str(value).strip() != expected:
                failures.append(f'saved Newton formulation is not {expected}')
        for field in ('u', 'q', 'u_tilde'):
            starts = [read_solution_field(handle, field, equations) for handle in old]
            ends = [read_solution_field(handle, field, equations) for handle in new]
            if not np.array_equal(starts[0], starts[1]):
                failures.append(f'Newton pair has different assembled initial {field}')
            changes = [end-start for end, start in zip(ends, starts)]
            components = changes if field != 'q' else [np.moveaxis(x, -2, -1) for x in changes]
            first, second = components
            checks = [calculate_error_norms(*changes)]
            checks.extend(calculate_error_norms(first[..., j], second[..., j])
                          for j in range(equations))
            passed = all(norms.available and norms.relative_l2 <= tolerances['update_relative_l2_max']
                         and norms.normalized_linf <= tolerances['update_normalized_linf_max']
                         for norms in checks)
            report[field] = {'passed': passed, 'relative_l2': checks[0].relative_l2,
                             'normalized_linf': checks[0].normalized_linf,
                             'equation_norms': [{'relative_l2': n.relative_l2,
                                                'normalized_linf': n.normalized_linf} for n in checks[1:]],
                             'absolute_update_l2': float(np.linalg.norm(changes[0].ravel()))}
            if not passed:
                failures.append(f'first Newton {field} updates differ beyond tolerance')
    return report
