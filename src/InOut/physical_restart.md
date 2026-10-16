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

Fixed and prescribed physical-time modes are enabled. Adaptive mode and the
reduction/growth options remain explicitly rejected until temporal recovery is
implemented. Predictor history also follows later.
