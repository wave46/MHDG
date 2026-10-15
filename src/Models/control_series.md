# Physical-time sources and feedback

With explicit `time_mode`, source schedules and the existing feedback controllers
are updated once before each timestep solve, independently of whether the magnetic
equilibrium moves. Curves are sampled at the proposed physical endpoint in seconds.
Feedback uses the numerical interval `time%dt * refval_time`, also in seconds.
Explicit solver execution remains guarded pending the restart/output slice.

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
recycling. No controller formulas, gains, saturation or anti-windup rules change.

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

Configurations without `time_mode` retain their existing ME-triggered source
updates, file requirements, and frame-based controller clock. Their interpolation
and input reading use the same cache and endpoint-safe locator.
