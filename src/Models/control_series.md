# Physical-time sources and feedback

With explicit `time_mode`, source schedules and the existing feedback controllers
are updated once before each timestep solve, independently of whether the magnetic
equilibrium moves. Curves are sampled at the proposed physical endpoint in seconds.
Feedback uses the numerical interval `time%dt * refval_time`, also in seconds.
At startup, prescribed actuator values and Zeff are sampled at the accepted
physical time. Feedback-owned actuator values and their histories are retained;
no PID update occurs until an advance is attempted.

`target_variable` retains its existing meanings:

| Value | Puff | Recycling | First impurity concentration, when radiation is enabled |
| --- | --- | --- | --- |
| 0 | Constant or prescribed | Configured value | Constant or prescribed |
| 1 | Density feedback | Configured value | Constant or prescribed |
| 2 | Constant or prescribed | Density feedback | Constant or prescribed |
| 3 | Density feedback | Configured value | XPR feedback |

For explicit modes, empty `puff_path`, `zeff_path`, and
`impurity_concentration_path` retain the corresponding configured constants.
A supplied path prescribes that quantity unless feedback owns it in the table.
Zeff may be constant or prescribed with every controller selection. A concentration
schedule or XPR feedback changes only entry 1; other mixture entries remain fixed.
The recycling controller changes `Re`, retaining the existing treatment of pump
recycling. The proportional/integral formulas, gains, saturation and anti-windup
rules are retained. The derivative uses the separate sampling intervals below.

Density feedback requires `target_density_path`; XPR feedback additionally requires
`target_density_xpr_path`. These supply targets, not actuator values. The existing
division of each stored density target by two and the existing line-integration
implementation are preserved.

Control files accept full `.h5` filenames or legacy stems. They contain rank-one
`time` and value datasets with matching lengths of at least two samples:

| Input path | Value dataset |
| --- | --- |
| `puff_path` | `puff` |
| `zeff_path` | `zeff` |
| `impurity_concentration_path` | `impurity_concentration` |
| Either target-density path | `target_density` |

Target-density files also contain the scalar chord coordinates `x_lower`,
`x_upper`, `y_lower`, and `y_upper` in metres. Times must be finite and strictly
increasing; values and chord coordinates must be finite. A curve must cover each
requested endpoint. Interpolation is linear, including exact first/final samples;
there is no extrapolation. Each curve keeps its own timestamps.

The solver reads each curve and its chord once, closes the file, and retains the
dimensional arrays in memory. Mesh adaptation retains these arrays and evaluates
feedback measurements on the current mesh. Inputs must remain unchanged during a
run. A changed path reloads that curve; normal shutdown releases the caches.

Mutable source/controller values remain in `phys`: puff, recycling, Zeff,
impurity concentrations, the two integral/previous-error pairs, and the last
line-integrated density. The curves contain no mutable feedback history. Temporal
retry will restore those values with the accepted solution; the current slice
does not yet implement rejected-step recovery.

## Derivative timing

At accepted plasma time `t`, the target is sampled at the proposed endpoint
`t+h`, while the measurement `y` comes from the accepted plasma at `t`.
The preceding controller update stored target `r_previous`, measurement
`y_previous`, and the physical time `t_previous` of that measurement. Use:

```text
D = Kd * ((r(t+h) - r_previous)/h - (y(t) - y_previous)/(t - t_previous))
```

The target increment spans the upcoming interval; the measurement increment
spans the previous accepted interval. Dividing both increments by the upcoming
interval creates artificial derivative spikes when that interval changes.
For equal intervals, this is algebraically the existing error-difference formula.
Both explicit and legacy controllers use this correction, including XPR feedback.
The first update has zero derivative because no preceding measurement exists.

Each channel has a small `Feedback_derivative_history_type` in `phys`, separate
from the integral/error scalars. Density feedback (puff or recycling) uses one
record and XPR feedback uses another. The records keep the previous target,
measurement, measurement time in seconds and an initialization flag. They must
be included when restoring rejected trial state, along with actuator values,
integral accumulators and previous errors. Changing timestep does not rescale
the accumulated integral; its tentative increment remains `error*h`.

Checkpoints store these records in `physics/feedback_derivative` and
`physics/feedback_derivative_xpr`. Older checkpoints reconstruct the preceding
target from the current target curve at checkpoint time, the measurement from
that target minus the saved previous error, and its time from the saved last
timestep. This assumes unchanged target history over the completed interval.
When an old file lacks usable timestep metadata, the first derivative is zero;
the actuator and integral history are retained. Startup reports either action.

Configurations without `time_mode` retain their existing ME-triggered source
updates, file requirements, and frame-based controller clock. Their interpolation
and input reading use the same cache and endpoint-safe locator.
