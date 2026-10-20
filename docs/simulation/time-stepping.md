# Physical time stepping

Physical transients advance to a final time with fixed, prescribed or adaptive
intervals. Numerical steps are independent of magnetic-frame indices. Fixed
equilibrium can also support a plasma transient.

Adaptive stepping responds to nonlinear-solve difficulty. It does not estimate
temporal error or limit changes in equilibrium and sources; choose intervals
and an upper bound that resolve the physics of interest.

## Choosing and configuring a mode

| `time_mode` | Numerical intervals | Failure reduction | Growth |
| --- | --- | --- | --- |
| `'fixed'` | Configured `dt0`, shortened at the final time | Not allowed | Not allowed |
| `'prescribed'` | Next magnetic timestamp, clipped at `t_final`; needs at least two equilibrium frames | Optional subdivision | Not allowed |
| `'adaptive'` | Start at `dt0`, retain reduced intervals and optionally grow | Optional | Optional |

For example, set these entries in the existing `param.txt` blocks:

```fortran
&SWITCH_LST
  steady = .false.
  psdtime = .false.
/
&TIME_LST
  time_mode = 'fixed'
  time_units = 'seconds'
  tis = 1
  dt0 = 0.02
  t_start = 0.0
  t_final = 0.08
/
```

All three physical modes require `steady = .false.`, `psdtime = .false.`
and BDF1 (`tis = 1`). They stop at `t_final`, even if a time-convergence
criterion is reached earlier; `nts` does not set their final step number.

| Entry in `TIME_LST` | Requirement or default |
| --- | --- |
| `time_mode` | One of the three physical modes; omission retains legacy stepping |
| `time_units` | `'seconds'` or `'normalized'`; default `'normalized'` |
| `t_start` | Initial time for a new run; default `0`; continuation restores checkpoint time |
| `t_final` | Required and greater than `t_start` |
| `dt0` | Positive for fixed/adaptive modes; prescribed intervals come from magnetic timestamps |

Times must be finite. Moving equilibrium, time-dependent heating and control
histories must cover their sampling times. One-frame equilibrium/heating
input is fixed at any finite time. See the
[equilibrium](../inputs/equilibrium.md#time-sampling),
[heating](../inputs/heating.md#time-sampling-and-restart) and
[controls](../inputs/controls.md#controller-timing) sampling rules.

### Optional reduction, growth and prediction

Reduction, growth and prediction default to inactive. Reduction and growth
are independent; their controls and the predictor settings belong in `TIME_LST`.

| Entry | Default | Meaning |
| --- | --- | --- |
| `dt_reduce_on_failure` | `.false.` | Retry recoverable failures with half the interval; prescribed/adaptive only |
| `dt_min` | `0` | Positive minimum required for reduction |
| `dt_grow_on_easy` | `.false.` | Grow after consecutive easy steps; adaptive only |
| `dt_max` | `0` | Positive upper bound required for growth; also caps adaptive steps without growth if supplied |
| `dt_growth_factor` | `1.2` | Growth multiplier; must exceed one when growth is enabled |
| `dt_growth_easy_steps` | `2` | Consecutive qualifying comparisons required; positive integer |
| `newton_predictor` | `.false.` | Extrapolate the initial Newton guess after sufficiently smooth evolution |
| `newton_predictor_rate_limit` | `0.35` | Allowed relative change between accepted rates; finite and nonnegative when prediction is enabled |

For adaptive reduction and growth, change the example's `time_mode` to
`'adaptive'` and add:

```fortran
  dt_reduce_on_failure = .true.
  dt_min = 1.e-6
  dt_grow_on_easy = .true.
  dt_max = 0.10
```

`dt_min` and `dt_max` must be finite and nonnegative. Adaptive `dt0` must
be at least `dt_min`; a positive `dt_max` must be at least both `dt0` and
`dt_min`.

## Time units and endpoints

`time_units` applies to namelist `dt0`, `dt_min`, `dt_max`, `t_start` and
`t_final`. With `'seconds'`, inputs are normalized once during initialization.
With `'normalized'`, supply solver time units:

$$
t_{\mathrm{physical}} = t_0\widehat{t},
\qquad
\Delta t_{\mathrm{physical}} = t_0\Delta\widehat{t}.
$$

The time reference $t_0$ is saved as
`/simulation_parameters/adimensionalization/time_scale`. Equilibrium,
heating and control-file timestamps are always in seconds. General terminal
time reporting and HDF5 fields such as `Final_time` remain in solver units;
see the [output time fields](restart-and-output.md#time-and-schedule-metadata).

Prescribed mode targets the next magnetic timestamp after the accepted time.
From 0.012 s with frames at 0, 0.02, 0.07 and 0.11 s and `t_final = 0.10`,
the endpoints are 0.02, 0.07 and 0.10 s. Fixed/adaptive steps may cross frame
timestamps, interpolating input at their proposed endpoints.

All physical modes shorten the last interval to reach `t_final`. A final
remainder may be shorter than `dt_min`; if it fails, it cannot be halved
below that minimum. See [restart and output](restart-and-output.md) for
continuation and saved timestep proposals.

## Failure recovery

With `dt_reduce_on_failure = .true.`, recoverable nonlinear failure restores
the last accepted mesh, plasma and controls before retrying half the interval.
Equilibrium, heating and controls are reevaluated at the new endpoint from
retained input data. Rejected trials do not become accepted output.

Existing Newton mesh-recovery triggers and budgets are tried within each
trial before temporal reduction. Newton remeshing projects the accepted
BDF state onto the new mesh, keeping it distinct from the current iterate.
Finite Newton-budget exhaustion also permits
reduction; fatal solver/backend errors terminate immediately. Halving stops
when another retry would fall below `dt_min`, or floating-point precision
cannot represent an intermediate time.

Prescribed mode retains failed endpoints during subdivision. Two failures
of 0→1 s and 0→0.5 s can produce accepted times 0.25, 0.5 and 1 s. After
reaching that frame, the next supplied frame is targeted in full. Adaptive
mode instead retains the successful reduced interval as its next proposal.
Fixed mode does not permit timestep adaptation.

## Timestep growth

An accepted step qualifies as easy when both its Newton iteration count and
first Newton correction are unchanged or smaller than in the preceding
comparable step. Comparisons require equal numerical intervals and clean
solves. The first clean step establishes a comparison; two qualifying
comparisons trigger growth with the defaults.

For example, `14 → 13 → 13` iterations qualify if the first correction does
not increase. Growth multiplies the proposal by
`dt_growth_factor`, capped at `dt_max`. The Newton budget is not a growth
threshold.

A changed interval starts a fresh comparison, preventing immediate regrowth
from a smaller correction caused by halving. Rejection, Newton remeshing,
scheduled time-level remeshing, predictor fallback and final clipping reset
history. A clean nonqualifying step resets the streak but becomes the next
comparison. This measures solve efficiency, without estimating temporal error.

## Optional Newton prediction

Set `newton_predictor = .true.` to extrapolate the initial guess. It works
with all three physical modes and fixed/moving equilibrium after two clean
accepted intervals establish history. It leaves the BDF1 equations unchanged.

For each field $y$ in the conservative variables, gradients and independent
face trace, the guess uses one common weight:

$$
y_{\mathrm{guess}} = y_n + w(y_n-y_{n-1}),
\qquad
w = \min\!\left(1,\frac{h}{t_n-t_{n-1}}\right).
$$

$h$ is the proposed interval. The cap limits extrapolation to one previous
accepted change when the timestep grows. Prediction is attempted only when
the two latest accepted rates change by no more than
`newton_predictor_rate_limit` (default `0.35`, or 35%):

$$
v_{\mathrm{previous}}=\frac{y_{n-1}-y_{n-2}}{t_{n-1}-t_{n-2}},
\qquad
v_{\mathrm{current}}=\frac{y_n-y_{n-1}}{t_n-t_{n-1}},
$$

$$
C_y=\frac{\lVert v_{\mathrm{current}}-v_{\mathrm{previous}}\rVert}
{\max(\lVert v_{\mathrm{current}}\rVert,\lVert v_{\mathrm{previous}}\rVert)}.
$$

The largest $C_y$ must satisfy the limit; if both rates vanish, $C_y=0$.
Otherwise ordinary Newton runs while history continues updating. This is a
convergence heuristic: recent smoothness does not guarantee an easy next solve.

A predicted solve needing nonlinear or mesh recovery first retries once
from the accepted plasma at the same endpoint and interval. Controls are not
advanced again, mesh-recovery budgets are not reset, and failed predicted
iterations count toward total Newton work. If the ordinary solve also fails,
reduction applies when enabled; fatal backend errors remain fatal. A fallback
cannot trigger growth.

Prediction can change finite-tolerance results and the nonlinear branch
reached. Compare runs at the same timesteps before combining it with growth.

## Implementation notes

### Growth metric

The first Newton correction is the absolute RMS conservative-variable update
divided by Newton damping, rather than an equation residual or temporal-error
estimate. Its history is [saved for continuation](restart-and-output.md#time-and-schedule-metadata).

### Predictor details

The arrays are `u` (conservative variables), `q` (gradients) and `u_hat`
(face trace, stored as `u_tilde`). Rate comparisons use Euclidean coefficient
norms separately per conservative field and array, combining gradient
components for that field. MPI sums count owned elements/faces once.

If a guessed array is nonfinite or fails existing owned-node density and
thermal-energy checks, the common weight is halved up to three times.
Optional neutral-positivity checks also apply. If no guess is admissible,
Newton uses the accepted state. BDF history `u0` is unchanged; feedback uses
accepted plasma before prediction, and boundary conditions are assembled on
every Newton iteration.

History is separate from BDF history: one preceding accepted `u/q/u_hat`
state and the preceding interval's rates. Remeshing, rejection or fallback
releases it; two new clean accepted intervals reestablish it. Reinitialization
and restart begin without predictor history. Checkpoints record settings,
not extra arrays; current input settings apply on continuation. Newton counts
and finite-tolerance results may therefore differ after restart.

## Legacy compatibility

Omitting `time_mode` retains the existing failure policy. Fixed-equilibrium
runs use step-count/convergence stopping; moving trajectories complete their
steps through `nts`, following the
[moving-frame/step rules](../inputs/equilibrium.md#legacy-time-stepping).
Reduction, growth and prediction require an explicit physical mode.

`tfi` remains a deprecated alias for `t_final` in explicit modes; conflicting
values are rejected. Without `time_mode`, it produces a warning and does not
enforce a final time. Supplying `t_final` or nonzero `t_start` without an
explicit mode is rejected.
