# Restart and solution output

A saved solution can initialize a new physical transient or continue an
existing one. Saved metadata selects the action, with the same rules for
fixed and moving magnetic equilibrium.

## Starting from a saved solution

Configure `param.txt`, then supply the saved solution as the second argument.
For example, from the executable's working directory:

```bash
./MHDG-NGammaTiTeNeutral-serial-2D /path/to/mesh.msh checkpoint.h5
```

Full filenames are preferred; legacy mesh and solution stems remain accepted.
To use the mesh stored in the solution, set `readMeshFromSol = .true.` in
`SWITCH_LST`; the first command-line argument is still required. Otherwise
supply the matching mesh. Projection onto a different mesh through an
additional command-line argument is unsupported.

| Saved solution | Action with an explicit physical `time_mode` |
| --- | --- |
| Steady or pseudo-time solution | Initialize at configured `t_start`, with numerical index zero and fresh feedback history |
| Physical transient checkpoint | Continue the saved physical time, numerical index and active feedback state |

No restart-mode switch is required. A physical step-1 checkpoint continues
normally; its index does not make it an initialization solution. Startup
reports the selected action. Changing `t_start` does not reset a physical
checkpoint's clock.

Use a compatible plasma model and reference scales: conservative arrays are
loaded in solver units, without conversion between reference normalizations.
BDF1 continues from the loaded accepted solution; earlier BDF states are not
reconstructed. See [time stepping](time-stepping.md) for mode configuration.

### Final time and prescribed inputs

Set `t_final` later than the checkpoint time to advance. A checkpoint beyond
`t_final` is rejected. One already at `t_final`, within roundoff tolerance,
exits without advancing or writing another solution, before equilibrium/
heating loading and normal startup control updates.

For active continuation, equilibrium, heating and prescribed controls are
sampled at the restored time, then at subsequent proposed endpoints. Input
histories must cover those times. Feedback-controlled actuators and history
are restored; configured constants and schedules follow the
[control-selection rules](../inputs/controls.md#configuration).
Feedback restart details, including older derivative-history reconstruction,
are in the [controls guide](../inputs/controls.md#restart-and-saved-history).

### Continuing the timestep schedule

Current timestep mode, limits and growth settings remain authoritative:

| Configured mode | Schedule restored from a physical checkpoint |
| --- | --- |
| `'fixed'` | Use configured `dt0`, clipped at `t_final` |
| `'prescribed'` | Restore pending subdivisions when failure reduction remains enabled; otherwise target the next supplied magnetic frame |
| `'adaptive'` | Restore the next proposed interval and, when enabled and available, growth history |

A shortened final interval does not replace the adaptive proposal. A run
with a 20 ms proposal and a 10 ms final remainder still proposes 20 ms when
continued to a later final time. A positive configured `dt_max` caps the
restored proposal; minimum-step rules also apply.

Prescribed checkpoints retain unfinished subdivisions. At 0.2125 s, pending
times stored as `[0.3, 0.25, 0.225]` are visited as 0.225, 0.25, 0.3 s.
A shortened `t_final` clips this schedule and clears it on completion.

Without saved schedule metadata, adaptive mode uses configured `dt0` and
prescribed mode targets the next supplied frame. Missing growth history
starts a fresh comparison. Steady/pseudo-time initialization does not restore
transient scheduling state. Newton prediction starts with fresh history
after restart; see its [continuation limits](time-stepping.md#predictor-details).

## Output configuration and filenames

Set `save_folder` in `INPUT_LST`, for example `'results/'` for an existing
directory. A positive `freqsave` in `UTILS_LST` writes every specified number
of accepted steps; startup and final output are independent of this frequency.
`saveNR` in `SWITCH_LST` enables intermediate Newton output.

For explicit physical modes, the filename's step label is the accepted
numerical index, with at least four digits:

| Output | Step label | Meaning |
| --- | --- | --- |
| Startup, including continuation | `_0000` | Initial output of this launch; metadata records the actual restored time and index |
| Accepted step 6 | `_0006` | Numerical step 6, regardless of save frequency or magnetic-frame index |
| Run finishes at accepted step 7 | `_0007` | Final solution is always written, even if `freqsave` skips step 7 |
| Accepted step 10000 | `_10000` | Indices expand beyond four digits |
| Newton iteration 2 | `_NR0002` | Intermediate iterate; overwritten at that iteration in subsequent steps |

For example, with `freqsave = 5` and a run ending at step 7, regular output
is written at step 5 and the final solution at step 7 (`_0007`). If the run
instead ends at step 10, the regular and final saves use the same `_0010`
filename.

The solver adds the mesh/model-dependent prefix and `.h5`; toroidal parallel
output can also have partition tags. Use accepted-step files for physical
continuation. Newton files contain intermediate iterates. HDF5 metadata
identifies physical time, independently of the filename.

Without explicit `time_mode`, output retains legacy sequence numbers (step
divided by save frequency), startup `_0000` and an unsuffixed final filename.
Legacy indices also expand beyond four digits.

## Time and schedule metadata

These datasets are under `/simulation_parameters/time`:

| Dataset | Meaning and units |
| --- | --- |
| `Physical_time_evolution` | Integer boolean: both saved `steady` and `psdtime` are false |
| `Current_time_seconds` | Saved clock in seconds; discharge time only for a physical transient |
| `Current_time` | Saved clock in solver units |
| `Current_time_step` | Saved numerical interval in solver units |
| `Initial_time_step`, `Final_time` | Configured initial interval and final time in solver units |
| `Current_time_step_number` | Numerical index, independent of equilibrium-frame count |
| `Time_mode` | Timestep mode; written for explicit modes |
| `Equilibrium_frame_left`, `Equilibrium_frame_right` | Zero-based bracketing frames; explicit modes with rectangular-grid equilibrium input |
| `Equilibrium_fraction` | Interpolation fraction between those frames |

Frames `6,7` and fraction `0.5` describe a midpoint sample. At a supplied
frame, the indices are equal and the fraction is zero; fixed equilibrium
uses frame `0`. Nodal input omits these provider fields. On restart, the
physical clock selects equilibrium position, independently of saved indices.

The same group stores continuation state:

| Dataset or subgroup | Written when | Contents |
| --- | --- | --- |
| `Proposed_time_step_seconds` | Adaptive mode | Next proposed interval before final-time clipping, in seconds |
| `Pending_time_targets_seconds` | Unfinished prescribed subdivisions | Remaining endpoints in seconds, next target last |
| `Growth_history` | Adaptive growth enabled, comparison available | `easy_steps`, `newton_iterations`, `first_newton_error`, `step_seconds` |
| `Newton_predictor`, `Newton_predictor_rate_limit` | Explicit modes | Predictor settings; extra solution history is not saved |

`Growth_history/step_seconds` is in seconds; the first error uses the
[growth metric](time-stepping.md#growth-metric). Current input selects the
predictor flag and rate limit on continuation. Feedback state is saved
separately under `/simulation_parameters/physics`.

## Older-file compatibility and implementation

When `Physical_time_evolution` is absent, classification uses saved flags:

| Older metadata | Interpretation |
| --- | --- |
| Either saved `steady` or `psdtime` is true | Initialize a new physical run |
| Both flags present and false | Continue a physical checkpoint |
| Flags incomplete | A true legacy moving-equilibrium marker identifies a physical checkpoint; otherwise initialize a new run |

When `Current_time_seconds` is absent, convert saved `Current_time` using
`/simulation_parameters/adimensionalization/time_scale`, falling back to the
current time reference if absent. This does not renormalize saved plasma.

When present, `Current_time_step` is also restored using the saved time
scale, preserving the checkpoint's interval in startup output. The next
proposal follows the [continuation schedule](#continuing-the-timestep-schedule).

Continuation restores the accepted clock before selecting another endpoint;
time-history storage grows with the numerical index. Runs without `time_mode`
retain their previous restart interpretation, including the legacy step-1
bootstrap convention and [frame/step rules](../inputs/equilibrium.md#legacy-time-stepping).
