"""Resolve a pinned limited bootstrap and derive moving control inputs."""

from pathlib import Path
import math
import shutil

import h5py
import numpy as np

from .bundles import artifact_path, load_manifest
from .compare_fixed import read_discrete_mesh, read_solution_field, required_scalar
from .documents import load_json
from .files import file_identity
from .support import BundleError


def bootstrap_pin(bundle, case_directory):
    """Identify the selected state; publication/relocation may preserve its bytes."""
    bundle = Path(bundle).resolve()
    _, manifest, identity = _selected_bootstrap(bundle, case_directory)
    return {'case_id': 'legacy_case', 'role': 'bootstrap_reference', 'solution': identity,
            'origin': {'bundle_id': manifest['bundle_id'], 'bundle_version': manifest['bundle_version'],
                       'manifest': file_identity(bundle/'manifest.json')}}


def resolve_bootstrap(descriptor, settings, case_directory):
    """Verify the selected dependency before linking its state into a run."""
    try:
        pin = load_json(descriptor, 'moving bootstrap dependency')
    except UnicodeError as exc:
        raise BundleError('moving bootstrap must be a JSON dependency descriptor') from exc
    if pin.get('case_id') != 'legacy_case' or pin.get('role') != 'bootstrap_reference':
        raise BundleError('moving bootstrap must identify legacy_case/bootstrap_reference')
    root = settings.get('MHDG_REGRESSION_BOOTSTRAP_BUNDLE')
    if not root:
        raise BundleError('select the limited bootstrap with --bootstrap-bundle or defaults.bundles.legacy_case')
    root = Path(root).resolve()
    state, _, identity = _selected_bootstrap(root, case_directory)
    if identity != pin.get('solution'):
        raise BundleError('selected limited bootstrap differs from the pinned moving initial state')
    return state, {'pin': pin, 'descriptor': file_identity(descriptor),
                   'source_manifest': {'path': str(root/'manifest.json'), **file_identity(root/'manifest.json')}}


def _selected_bootstrap(bundle, case_directory):
    """Verify this dependency without loading unrelated feature references."""
    manifest = load_manifest(bundle, case_directory, case_id='legacy_case')
    if 'bootstrap_reference' not in manifest['roles']:
        raise BundleError('limited bundle needs bootstrap_reference')
    record = manifest['artifacts'][manifest['roles']['bootstrap_reference']]
    state = artifact_path(bundle, record['path'], 'bootstrap_reference')
    identity = file_identity(state)
    if identity != {key: record[key] for key in ('sha256', 'size_bytes')}:
        raise BundleError('limited bootstrap checksum/size differs from its manifest')
    return state, manifest, identity


def prepare_bootstrap_controls(bootstrap, artifacts, controls, output):
    """Anchor main/short controls to the same original state, including restarts."""
    measured = measure_initial_line_density(bootstrap, artifacts['control_target_density'])
    with h5py.File(bootstrap) as handle:
        puff = required_scalar(handle, 'simulation_parameters/physics', 'puff')
    if not math.isfinite(puff) or puff <= 0:
        raise BundleError('bootstrap puff must be finite and positive')
    record = {'puff_particles_per_s': puff, 'measured_line_density_per_m2': measured, 'curves': {}}
    generated = {}
    for name, start in (('puff', puff), ('target_density', 2*measured)):
        role = controls.get(name)
        if role not in (f'control_{name}', f'short_control_{name}'):
            continue
        time, values = _control(artifacts[f'control_{name}'], name)
        if time[0] != 0 or values[-1] == values[0]:
            raise BundleError(f'{name} needs a nonconstant ramp starting at time zero')
        progress = (values-values[0])/(values[-1]-values[0])
        anchored = start + progress*(values[-1]-start)
        anchored[0], anchored[-1] = start, values[-1]
        if role.startswith('short_'):
            short_time, _ = _control(artifacts[role], name)
            if short_time[0] < time[0] or short_time[-1] > time[-1]:
                raise BundleError(f'short {name} timestamps are outside the main curve')
            anchored = np.interp(short_time, time, anchored)
        target = output/f'{name}.h5'
        shutil.copy2(artifacts[role], target)
        with h5py.File(target, 'r+') as handle:
            handle[name][...] = anchored
        generated[f'inputs/{name}.h5'] = file_identity(target)
        record['curves'][name] = {'old_start': float(values[0]), 'new_start': float(start),
                                  'final_target': float(values[-1]), 'final_time_s': float(time[-1])}
    return generated, record


def _control(path, dataset):
    try:
        with h5py.File(path) as handle:
            time, values = (np.asarray(handle[key][()], dtype=float) for key in ('time', dataset))
    except (OSError, KeyError, ValueError, TypeError) as exc:
        raise BundleError(f'cannot read {dataset} control curve {path}: {exc}') from exc
    if (time.ndim != 1 or len(time) < 2 or values.shape != time.shape
            or not np.isfinite(time).all() or not np.isfinite(values).all() or np.any(np.diff(time) <= 0)):
        raise BundleError(f'invalid {dataset} control curve: {path}')
    return time, values




def measure_initial_line_density(seed, target_curve):
    """Use the controller's 1000-point nearest-node chord quadrature.

    Mirrors compute_line_integrated_density/lineintegration, including the
    native fixed length conversion and the one-sided mesh-entry/exit handling.
    """
    try:
        with h5py.File(target_curve) as handle:
            x0, x1, y0, y1 = (float(np.asarray(handle[key][()]).item()) for key in
                               ('x_lower', 'x_upper', 'y_lower', 'y_upper'))
    except (OSError, KeyError, ValueError, TypeError) as exc:
        raise BundleError(f'cannot read density-control chord {target_curve}: {exc}') from exc
    with h5py.File(seed) as handle:
        mesh = read_discrete_mesh(handle)
        equations = int(required_scalar(handle, 'simulation_parameters', 'Neq'))
        density = read_solution_field(handle, 'u', equations)[..., 0]
        length = required_scalar(handle, 'simulation_parameters/adimensionalization', 'length_scale')
        scale = required_scalar(handle, 'simulation_parameters/adimensionalization', 'density_scale')
    if not all(math.isfinite(value) and value > 0 for value in (length, scale)):
        raise BundleError('bootstrap density and length scales must be positive')
    points = np.column_stack((x0 + (x1-x0)*np.arange(1000)/1000,
                              y0 + (y1-y0)*np.arange(1000)/1000)) / length
    coordinates = mesh['X'].T / length
    connectivity = mesh['T'].T.astype(int) - 1
    samples, inside = np.zeros(1000), np.zeros(1000, dtype=bool)
    for element, nodes in enumerate(connectivity):
        a, b, c = coordinates[nodes[:3]]
        area = _triangle_area(a, b, c)
        subareas = _triangle_area(points, b, c) + _triangle_area(a, points, c) + _triangle_area(a, b, points)
        selected = np.flatnonzero(~inside & (np.abs(area - subareas) < 1e-10))
        if len(selected):
            distance = np.sum((points[selected, None] - coordinates[nodes][None])**2, axis=2)
            samples[selected] = density[element, np.argmin(distance, axis=1)]
            inside[selected] = True
    dl = np.linalg.norm(np.diff(points, axis=0), axis=1) * 1.901e-3
    measured = float(np.sum((0.5*(samples[1:]+samples[:-1])*dl)[inside[1:]]) * scale)
    if not math.isfinite(measured) or measured <= 0:
        raise BundleError('bootstrap chord measurement must be finite and positive')
    return measured


def _triangle_area(a, b, c):
    return 0.5 * np.abs(a[..., 0]*(b[..., 1]-c[..., 1]) +
                        b[..., 0]*(c[..., 1]-a[..., 1]) + c[..., 0]*(a[..., 1]-b[..., 1]))
