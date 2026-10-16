# Physical-time checkpoint interpretation

For explicit `time_mode`, supplying a saved solution has two meanings, selected
automatically from the saved metadata:

- A saved steady or pseudo-time solution initializes the plasma for a new run
  at configured `t_start`, with numerical index zero and fresh feedback history.
- A physical transient checkpoint continues its saved time, numerical index and
  active feedback state, with either fixed or moving magnetic equilibrium.

There is no required restart-mode switch. A saved physical step 1 is a checkpoint;
its index does not identify it as steady initialization. Startup prints the
selected action. Existing legacy configurations retain their previous restart
interpretation, including the old step-1 bootstrap convention.

New output adds two datasets to `/simulation_parameters/time`:

| Dataset | Meaning |
| --- | --- |
| `Physical_time_evolution` | Integer boolean: the saved run has `steady=false` and `psdtime=false` |
| `Current_time_seconds` | Saved clock converted to seconds using the saved run's time reference |

Explicit modes also write `Time_mode`. When a rectangular-grid equilibrium
provider is active, they write `Equilibrium_frame_left`, `Equilibrium_frame_right`
(zero-based input frame indices) and `Equilibrium_fraction`. For example,
`left=6`, `right=7`, `fraction=0.5` means halfway between supplied frames 6 and 7.
At a supplied frame, the indices are equal and the fraction is zero. A fixed
equilibrium has both indices zero. These fields describe the trajectory position;
the physical clock remains authoritative on restart, so numerical indices never
select equilibrium frames. Nodal inputs omit this provider metadata.
The actuator value is saved in `physics/puff`, with feedback histories saved
separately.

The same time group stores the numerical schedule needed for continuation:

| Dataset | Written when | Meaning |
| --- | --- | --- |
| `Proposed_time_step_seconds` | Adaptive mode | Next proposed interval, before clipping to the final time |
| `Pending_time_targets_seconds` | Prescribed mode has subdivided an interval | Remaining target times, with the next target last |

For example, a checkpoint at 0.2125 s with pending times `[0.3, 0.25, 0.225]`
continues toward 0.225 s, then 0.25 s and 0.3 s. The numerical step index does
not select any of these times. Pending targets are restored when prescribed
mode and failure reduction remain enabled. Otherwise the configured mode
selects the next endpoint normally. A newly shortened `t_final` clips pending
targets too; reaching it clears the remaining targets.

Adaptive continuation restores its saved proposal instead of resetting to
`dt0`. A shortened final interval does not overwrite that proposal: a run with
a 20 ms proposal and a 10 ms final remainder still proposes 20 ms when restarted
with a later final time. Configured timestep limits remain enforced. Older
files without schedule metadata start from configured `dt0` in adaptive mode,
or the next supplied frame in prescribed mode. Saved steady/pseudo-time states
continue to initialize a new run and do not restore transient scheduling state.

The existing `Current_time`, `Current_time_step_number`, and other datasets keep
their names and units. The step number is the numerical index, independently of
the equilibrium frame count. For nonphysical runs, the classification flag
prevents treating their scaled clock as elapsed physical discharge time.

When the classification dataset is absent, saved `steady`/`psdtime` flags identify
old steady, pseudo-time and physical runs. If those flags are incomplete, a saved
true `ME` flag identifies a physical checkpoint; otherwise the solution initializes
a new run. The normalized legacy `Current_time` is converted using saved
`adimensionalization/time_scale`, falling back to the current time reference when
that scale is absent. New files use `Current_time_seconds` directly.

Continuation restores the accepted clock and sets the proposed endpoint initially
to that clock. The next-step selector then chooses a later endpoint or reports
completion at `t_final`; checkpoint times beyond `t_final` are rejected. Time
history storage grows to accommodate the numerical index. BDF1 starts from the
loaded accepted solution; earlier BDF history is not required or reconstructed.
The existing solution reader loads conservative variables in solver units, so
their plasma reference scales must remain compatible.

If the checkpoint has already reached `t_final`, startup exits without loading
equilibrium/heating inputs, updating controllers or writing another solution.

Derivative history is stored with each feedback channel so continuation preserves
the next derivative contribution when timestep changes. Older-file reconstruction
and the timing convention are described in the
[control guide](../Models/control_series.md#derivative-timing).

Feedback histories and the actuators owned by feedback are saved/restored for
explicit modes independently of `ME`. Controller formulas and configured gains
are retained. Configured constants and source schedules keep their existing
ownership rules from the [control guide](../Models/control_series.md).

For explicit physical-time modes, solution filenames end with the accepted
numerical index, with at least four digits: `_0000`, `_0006`, `_10000`. The save
frequency selects which steps are written; it does not divide the filename index.
Startup output remains a service file named `_0000`, including on continuation;
its metadata still records the restored physical time and numerical index.
Final output uses the same indexed name as a regular output of that step,
including when the final step falls outside the save frequency. Newton diagnostics keep their existing
Newton-only names, such as `_NR0002`, and overwrite those files on subsequent
steps. HDF5 timestamps remain authoritative.

Legacy configurations retain their output sequence numbers (step divided by save
frequency), startup `_0000` and unsuffixed final filename. All indices can grow
beyond four digits.

Fixed, prescribed and adaptive physical-time modes are enabled for BDF1.
Reduction is opt-in with `dt_reduce_on_failure=true`; it requires a positive
`dt_min`, in the same `time_units` as `dt0`. Repeated halving stops when the next
retry would be smaller than this minimum; no separate rejection limit is needed.
Existing Newton mesh-recovery triggers and budgets are tried within each trial;
recoverable failure then restores the last accepted mesh, plasma and controllers
before halving. Finite Newton-budget exhaustion also permits temporal retry.
Fatal solver errors terminate immediately.

Prescribed mode retains pending endpoints during subdivision: two failures of
0→1 and 0→0.5 can produce accepted times 0.25, 0.5 and 1. After reaching that
frame, the next supplied frame is targeted in full. Adaptive mode starts at
`dt0`, retains a successful reduced interval, and may cross equilibrium frames.
Both sample equilibrium and controls at each trial endpoint in memory. Fixed
mode retains its configured timestep and does not permit adaptation.

A final remainder may be shorter than `dt_min` to reach `t_final` exactly, but a
failed remainder cannot be halved below `dt_min`. The code also stops if rounding
prevents representing a time strictly between the accepted and failed endpoints.
These limits terminate without saving the rejected trial as an
accepted step. Legacy inputs retain their existing stepping and failure policy.
The Newton predictor remains a later slice.

Optional growth is enabled with `time_mode='adaptive'`,
`dt_grow_on_easy=true` and an explicit positive `dt_max`. Reduction and growth
are independent options. For example:

```fortran
&TIME_LST
    time_mode = 'adaptive'
    time_units = 'seconds'
    dt0 = 0.02
    t_final = 1.36
    dt_reduce_on_failure = .true.
    dt_min = 1.e-6
    dt_grow_on_easy = .true.
    dt_max = 0.10
/
```

| Setting | Default | Meaning |
| --- | --- | --- |
| `dt_growth_factor` | 1.2 | Multiply the proposal after an easy-step streak, capped by `dt_max` |
| `dt_growth_easy_steps` | 2 | Consecutive easy steps required |

A clean accepted step qualifies when both its Newton iteration count and its
first NR correction are unchanged or smaller than in the previous accepted
step. Thus `14 -> 13 -> 13` and `15 -> 15 -> 15` supply two easy steps with the
default settings if the first correction does not increase. Comparisons require
the same numerical interval; the configured Newton budget does not determine
growth eligibility. The first correction is the existing
absolute RMS conservative-variable update divided by Newton damping; it is not
an equation residual or a timestep-error estimate.

The first clean step establishes a comparison. A changed interval establishes a
new comparison too, so the smaller first correction caused by halving cannot
immediately trigger regrowth. Rejection, Newton remeshing, scheduled time-level
remeshing and final-step clipping reset the easy-step history. A clean step that
does not qualify resets the streak but remains the next comparison. No absolute
Newton-count threshold triggers growth.

Adaptive checkpoints save `Growth_history` under the time group when growth is
enabled and a comparison is available. It contains `easy_steps`,
`newton_iterations`, `first_newton_error` and `step_seconds`. These restore the
growth decision across restart; older files without the history start with a
fresh comparison.
Current configured limits and growth settings remain authoritative. A saved
proposal above a newly configured `dt_max` is capped on loading.

Growth measures nonlinear-solve efficiency. It does not estimate temporal error
or constrain the rate of change of equilibrium and control signals. Choose
`dt_max` to resolve those dynamics; an accuracy estimator remains separate work.
