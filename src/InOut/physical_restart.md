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

Derivative history is stored with each feedback channel so continuation preserves
the next derivative contribution when timestep changes. Older-file reconstruction
and the timing convention are described in the
[control guide](../Models/control_series.md#derivative-timing).

Feedback histories and the actuators owned by feedback are saved/restored for
explicit modes independently of `ME`. Controller formulas and configured gains
are retained. Configured constants and source schedules keep their existing
ownership rules from the [control guide](../Models/control_series.md).

This checkpoint slice leaves explicit solver execution guarded. Output naming,
initial/restart output handling and end-to-end loop activation follow in the next
slice; temporal rejection, timestep growth and predictor history follow later.
