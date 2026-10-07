"""Verify production ME restart step counts and checkpoint continuation on a prepared case.

The source case is read-only. Only frozen-field test fixtures in --run-root are written.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import h5py
import numpy as np

parser = argparse.ArgumentParser(
    description='Check ME step limits with frozen fields and controls from a prepared external case.')
parser.add_argument('--case-directory', type=Path, required=True)
parser.add_argument('--solver', type=Path, required=True)
parser.add_argument('--run-root', type=Path, required=True)
parser.add_argument('--mpi-ranks', type=int, default=2)
args = parser.parse_args()
REPO = Path(__file__).resolve().parents[1]
CASE = args.case_directory.resolve()
ROOT = args.run_root.resolve()
SOLVER = args.solver.resolve()
if args.mpi_ranks < 1:
    parser.error('--mpi-ranks must be positive')
with h5py.File(CASE/'magnetic/equilibrium_0001.h5') as f:
    frame_dt = float(np.asarray(f['dt'][()]).ravel()[0])
with h5py.File(CASE/'inputs/restart.h5') as f:
    time_scale = float(f['simulation_parameters/adimensionalization/time_scale'][0])
if ROOT.exists():
    raise SystemExit('Check directory already exists; use a new ROOT.')
ROOT.mkdir()
(ROOT/'magnetic').mkdir()
(ROOT/'controls').mkdir()
for index in [1,2,3,51,52,53,100,101]:
    for prefix in ['equilibrium','current_density']:
        destination = ROOT/f'magnetic/{prefix}_{index:04d}.h5'
        shutil.copy2(CASE/f'magnetic/{prefix}_0001.h5',destination)
        with h5py.File(destination,'r+') as f:
            f['time'][...] = (index-1)*frame_dt
for name in ['target_density','puff','zeff','impurity_concentration']:
    destination = ROOT/f'controls/{name}.h5'
    shutil.copy2(CASE/f'controls/{name}.h5',destination)
    with h5py.File(destination,'r+') as f:
        f[name][...] = f[name][0]
(ROOT/'impurity_model.nml').write_text(
    (CASE/'inputs/impurity_model.nml').read_text().replace(str(CASE/'controls'), str(ROOT/'controls')))
template = (CASE/'runs/full/param.txt').read_text()


def prepare(name, step, nts, target=0, restart=None):
    run = ROOT/name
    run.mkdir()
    (run/'output').mkdir()
    (run/'positionFeketeNodesTri2D.h5').symlink_to(REPO/'test/positionFeketeNodesTri2D.h5')
    if restart is None:
        restart = run/'restart.h5'
        shutil.copy2(CASE/'inputs/restart.h5',restart)
        if step:
            with h5py.File(restart,'r+') as f:
                p=f['simulation_parameters']
                p['switches/ME'][...] = 1
                p['time/Current_time_step_number'][...] = step
                p['time/Current_time'][...] = step*frame_dt/time_scale
                for key,value in [('feedback_integral_error',123.),('feedback_previous_error',456.)]:
                    p['physics'].create_dataset(key,data=[value])
    changes = {
        'nts': str(nts), 'target_variable': str(target),
        # Step-limit verification uses a fixed mesh even when the source case
        # enables oscillation-triggered refinement for its future physics run.
        'adaptivity': '.false.',
        'field_path': repr(str(ROOT/'magnetic/equilibrium')),
        'jtor_path': repr(str(ROOT/'magnetic/current_density')),
        'puff_path': repr(str(ROOT/'controls/puff.h5')),
        'target_density_path': repr(str(ROOT/'controls/target_density.h5')),
        'zeff_path': repr(str(ROOT/'controls/zeff.h5')),
        'impurity_model_path': repr(str(ROOT/'impurity_model.nml')),
        'save_folder': repr(str(run/'output')+'/'), 'freqsave': '1',
    }
    if step == nts:
        # A completed restart must succeed without any magnetic input files.
        changes['field_path'] = repr(str(ROOT/'absent/equilibrium'))
        changes['jtor_path'] = repr(str(ROOT/'absent/current_density'))
    text = template
    for key, value in changes.items():
        text, count = re.subn(rf'(?im)^(\s*{key}\s*=).*$', rf'\g<1> {value}', text)
        assert count == 1, (key, count)
    (run/'param.txt').write_text(text)
    return run, Path(restart).with_suffix('')


env = os.environ.copy()
env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
results = []


def execute(name, step, nts, expected_steps, target=0, restart=None):
    run, restart = prepare(name, step, nts, target, restart)
    with (run/'solver.log').open('w') as log:
        proc = subprocess.run(
            ['mpirun.openmpi', '-n', str(args.mpi_ranks), str(SOLVER), 'me_nts', str(restart)],
            cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=300)
    output=(run/'solver.log').read_text()
    actual=[int(x) for x in re.findall(r'Time iteration\s*=\s*(\d+)',output)]
    frames=[int(x) for x in re.findall(r'Magnetic field loaded from file:.*equilibrium_(\d+)\.h5',output)]
    result = {'name': name, 'returncode': proc.returncode, 'advances': actual,
              'loaded_frames': frames, 'log': str(run/'solver.log')}
    results.append(result)
    if expected_steps is None:
        assert proc.returncode != 0 and 'outside [0, nts' in output and not actual, result
    else:
        assert proc.returncode == 0 and actual == expected_steps, result
        expected_frames = ([step+1] + [s+1 for s in expected_steps if s < nts]) if expected_steps else []
        assert frames == expected_frames, result
        assert 'HDF5-DIAG' not in output, result
        if not expected_steps:
            assert 'ME final step already reached' in output, result
            assert not list((run/'output').iterdir()), result
            assert 'INITIALIZING PUFF' not in output, result
        else:
            finals = [p for p in (run/'output').glob('Sol2D_*.h5') if not re.search(r'_\d{4}\.h5$',p.name)]
            assert len(finals) == 1, finals
            final = finals[0]
            with h5py.File(final) as f:
                assert f['simulation_parameters/time/Current_time_step_number'][0] == nts
                assert np.isfinite(f['solution/u'][()]).all()
            result['final_solution']=str(final)
    print(json.dumps(result), flush=True)
    (ROOT/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    return result

execute('completed_100',100,100,[],target=1)
execute('beyond_100',100,99,None,target=1)
fresh=execute('fresh_2',0,2,[1,2])
middle=execute('restart_50',50,52,[51,52])
execute('restart_99',99,100,[100])
first_checkpoint = next((ROOT/'fresh_2/output').glob('*_0001.h5'))
resumed=execute('resume_1',1,2,[2],restart=first_checkpoint)
comparisons=[]
for other in [middle,resumed]:
    with h5py.File(fresh['final_solution']) as f,h5py.File(other['final_solution']) as g:
        errors={}
        for key in ['u','q','u_tilde']:
            a,b=f[f'solution/{key}'][()],g[f'solution/{key}'][()]
            error=float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-30))
            # Sparse factorization repeatability is weaker for the gradient.
            # Keep this far tighter than the nonlinear convergence tolerance.
            assert error<1e-8,(other['name'],key,error)
            errors[key]=error
        comparisons.append({'against':other['name'],'relative_errors':errors})
(ROOT/'comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n')
print(json.dumps(comparisons),flush=True)

feedback = execute('feedback_fresh_2',0,2,[1,2],target=1)
feedback_checkpoint = next((ROOT/'feedback_fresh_2/output').glob('*_0001.h5'))
feedback_resumed = execute('feedback_resume_1',1,2,[2],target=1,restart=feedback_checkpoint)
with h5py.File(feedback['final_solution']) as f,h5py.File(feedback_resumed['final_solution']) as g:
    errors = {}
    density_scale = float(f['simulation_parameters/physics/n_li'][0])
    for group,key in [('solution','u'),('solution','q'),('solution','u_tilde'),
                      ('simulation_parameters/physics','puff'),
                      ('simulation_parameters/physics','feedback_integral_error'),
                      ('simulation_parameters/physics','feedback_previous_error')]:
        a,b=f[f'{group}/{key}'][()],g[f'{group}/{key}'][()]
        scale = max(np.linalg.norm(a), 1e-30)
        # Controller errors subtract nearly equal line densities. Normalize
        # by their physical density/time scales, not the small residual itself.
        if key == 'feedback_previous_error':
            scale = max(scale, density_scale)
        elif key == 'feedback_integral_error':
            scale = max(scale, density_scale * frame_dt)
        error = float(np.linalg.norm(a-b) / scale)
        assert error<1e-8,(key,error)
        errors[key]=error
    comparisons.append({'against':'feedback_resume_1','normalized_errors':errors})
(ROOT/'comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n')
print('ME step-limit executable checks: 8 PASS',flush=True)
