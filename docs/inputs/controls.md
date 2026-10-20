# Source schedules and density feedback

Prescribe gas puffing, effective charge and the first impurity concentration
as functions of time, or use density feedback to adjust puff, recycling or
impurity concentration. With an explicit `time_mode`, these controls work
with fixed or moving magnetic equilibrium.

## Configuration

Choose the control arrangement with `target_variable` in `SWITCH_LST`:

| `target_variable` | Puff | Recycling | First impurity concentration, when radiation is enabled |
| --- | --- | --- | --- |
| `0` | Constant or prescribed | Configured value | Constant or prescribed |
| `1` | Density feedback | Configured value | Constant or prescribed |
| `2` | Constant or prescribed | Density feedback | Constant or prescribed |
| `3` | Density feedback | Configured value | X-point-region density feedback (XPR) |

**Moving equilibrium does not require control schedules.** With explicit
`time_mode`, leave a schedule path empty to use its configured constant
where allowed below. Feedback overrides schedules for quantities it controls.

| Quantity | Constant setting | Schedule path | Can use a constant? |
| --- | --- | --- | --- |
| Puff | `puff` | `puff_path` | Yes for `0,2`; feedback for `1,3` |
| Effective charge | `Zeff` | `zeff_path` | Yes for every selection |
| Recycling | `Re` | None | Yes for `0,1,3`; feedback for `2` |
| Pump recycling | `Re_pump` | None | Always; separate from recycling feedback |
| First impurity concentration | `impurity_concentrations(1)` | `impurity_concentration_path` | Yes for `0,1,2`; XPR feedback for `3` |
| Other impurity concentrations | `impurity_concentrations(2:)` | None | Always |

In `param.txt`, constants belong in `PHYS_LST` and paths in `INPUT_LST`.
The exception is impurity configuration: concentrations and their schedule
path belong in `IMPURITY_RADIATION_LST`, selected by `impurity_model_path`.
Concentrations are fractions of electron density, used when impurity radiation
is enabled. See the [sample impurity namelist](../../test/impurity_model.nml).

Feedback targets have no constant-value namelist alternative:

| Target path in `INPUT_LST` | Required when | Otherwise |
| --- | --- | --- |
| `target_density_path` | `target_variable = 1,2,3` | Unused; can be omitted |
| `target_density_xpr_path` | `target_variable = 3` and impurity radiation is enabled | Unused; can be omitted |

Targets use the [two-way chord convention](#density-target-convention).
Paths are relative to the working directory unless absolute. Runs without
`time_mode` retain the
[legacy file requirements](#implementation-and-legacy-compatibility).

### Feedback gains

The [feedback equation](#feedback-update) adds proportional, integral and
derivative contributions to the current actuator value. Set its gains in
`PHYS_LST`:

| Contribution | Puff or recycling feedback | XPR impurity feedback |
| --- | --- | --- |
| $K_p$ (proportional) | `feedback_propotional_gain` | `feedback_propotional_gain_xpr` |
| $K_i$ (integral) | `feedback_integral_gain` | `feedback_integral_gain_xpr` |
| [$K_d$ (derivative)](#derivative-timing) | `feedback_derivative_gain` | `feedback_derivative_gain_xpr` |

The spelling `propotional` is part of the input format. Gains use dimensional
density errors and seconds; choose their values for the actuator's units and
scale.

## HDF5 format

Normally, supply **separate files**, for example `puff_path = 'puff.h5'` and
`zeff_path = 'zeff.h5'`. Each needs `/time` and the dataset selected by its path:

| Path | Value dataset | Quantity and units |
| --- | --- | --- |
| `puff_path` | `/puff` | Puff rate in particles/s, in the same convention as `PHYS_LST`'s `puff` |
| `zeff_path` | `/zeff` | Effective charge $Z_{\mathrm{eff}}$, dimensionless |
| `impurity_concentration_path` | `/impurity_concentration` | First impurity's fraction of electron density, dimensionless |
| `target_density_path` or `target_density_xpr_path` | `/target_density` | Two-way chord-integrated density target in m⁻² |

`/time` contains physical timestamps in seconds. Both arrays have shape
`(nt,)` in Python/HDF5 notation, with matching lengths and at least two
samples, even for a constant curve. All datasets are at the file root.

**Shared files are supported:** point several paths to one file containing
their value datasets and a common `/time` vector. Both density-target paths
read `/target_density` and the chord fields below: different targets or
chords require separate files.

Density-target files additionally require the measurement chord:

| Dataset | Shape | Quantity and units |
| --- | --- | --- |
| `/x_lower` | Scalar | First endpoint's radial coordinate in metres |
| `/x_upper` | Scalar | Second endpoint's radial coordinate in metres |
| `/y_lower` | Scalar | First endpoint's vertical coordinate in metres |
| `/y_upper` | Scalar | Second endpoint's vertical coordinate in metres |

All values and chord coordinates must be finite; timestamps must be finite
and strictly increasing. Lengths are inferred; legacy `*_dimension` inputs
do not set array sizes. Full `.h5` filenames and legacy stems are accepted.
File timestamps remain in seconds regardless of the namelist's `time_units`.
A requested file with missing or unreadable required datasets stops the run;
omitting an optional schedule path is the way to retain its configured constant.

Separate files may use independent timestamps, unrelated to magnetic or
heating input. Sampling is linear, including exact endpoints. Cover startup,
continuation and subsequent requested times; there is
[no extrapolation beyond endpoint tolerance](equilibrium.md#time-sampling).

## Density-target convention

Supply a **two-way chord-integrated density signal**, as in a double-pass
experimental measurement. The controller divides this supplied signal by two
before comparing it with the simulated one-way ion-density integral:

$$
r(t) = \frac{N_{\mathrm{target,file}}(t)}{2},
\qquad
y(t) \simeq \int_{\mathrm{chord}} n_i(R,Z,t)\,\mathrm{d}\ell,
\qquad
e = r-y.
$$

Here $n_i$ is the evolved ion density and $r,y$ are in m⁻². The chord is the
straight segment between the supplied endpoints. Double a one-way target
before storing it in `/target_density`. This convention applies to both
target files; revisiting it remains open for a future controller refactor.

## Controller timing

Before each attempted timestep, prescribed sources and targets are sampled
at the proposed endpoint $t+h$, where $t$ is the accepted plasma time and
$h$ is the proposed interval in seconds. Feedback uses the accepted plasma
at $t$, giving error $e=r(t+h)-y(t)$.

Startup and continuation sample prescribed values at the accepted time.
Feedback-controlled actuators and their restored histories are retained
until a timestep is attempted.

Rejected steps restore actuators, integral accumulators, previous errors and
derivative histories with the accepted solution. Retries evaluate the new
endpoint without accumulating rejected feedback updates.

### Feedback update

Each controller adds its contributions to the current actuator value $a$:

$$
a_{\mathrm{trial}} = a + K_p e + K_i I + D,
\qquad
a_{\mathrm{new}} = \max(a_{\mathrm{trial}},0).
$$

$I$ is the stored integral accumulator and $D$ is the derivative contribution
defined below. Actuators have a lower bound of zero; the controller imposes
no upper bound. Integral back-calculation handles lower saturation. When the
new actuator is positive, the integral also accumulates $e h$. Changing the
timestep does not rescale the accumulated history.

## Derivative timing

The target increment spans the upcoming interval $h$, while the measured
density increment spans the preceding accepted interval. If the previous
controller evaluation stored target $r_{\mathrm{prev}}$, measurement
$y_{\mathrm{prev}}$ and measurement time $t_{\mathrm{prev}}$, the derivative
contribution is

$$
D = K_d\left[
\frac{r(t+h)-r_{\mathrm{prev}}}{h}
-\frac{y(t)-y_{\mathrm{prev}}}{t-t_{\mathrm{prev}}}
\right].
$$

Separate intervals avoid artificial derivative spikes when the timestep
changes; equal intervals recover the error-difference formula. The first
derivative is zero. Puff and recycling use the density-channel history;
XPR has a separate history. Legacy controllers use the same derivative timing.

## Restart and saved history

Checkpoints preserve the actuator values controlled by feedback and each
channel's integral, previous error and derivative history. General
continuation rules are in the [restart guide](../simulation/restart-and-output.md).

Derivative records are stored under
`/simulation_parameters/physics/feedback_derivative` and, for the XPR channel,
`/simulation_parameters/physics/feedback_derivative_xpr`:

| Dataset within each record | Meaning |
| --- | --- |
| `initialized` | Whether a preceding measurement is available |
| `target` | Previously sampled target, after division by two |
| `measurement` | Previously measured one-way ion-density integral |
| `measurement_time_seconds` | Physical time of that measurement |

Older checkpoints reconstruct the preceding target from the current curve
at checkpoint time, the measurement as that target minus the saved previous
error, and its time as checkpoint time minus the saved last timestep. Keep
the target history unchanged over that interval. Without usable timestep
metadata, the first derivative is zero; actuator and integral history are
retained. Startup reports either action.

## Implementation and legacy compatibility

Curves and chords are cached in dimensional units until shutdown. Keep files
unchanged during a run: only a changed path reloads a curve. Mesh adaptation
retains curves and evaluates density on the current mesh.

The chord measurement currently samples 1000 points using the nearest node
in the containing element and a trapezoidal sum. Outside-mesh points do not
contribute, and MPI sums omit ghost elements. Its length conversion uses the
fixed 1.901 × 10⁻³ m reference, matching the current physical normalization;
this is an implementation constraint if the normalization is changed.

The stored integral receives a back-calculation increment
$\gamma(a_{\mathrm{new}}-a_{\mathrm{trial}})$, with $\gamma=0.01$ for puff and
$\gamma=0.1$ for recycling and impurity concentration, before the conditional
$e h$ increment above.

Without `time_mode`, source updates and the controller clock remain tied to
moving magnetic frames. Empty schedule paths do not provide the constant
fallback above. File reading and interpolation use the same cached curves
and endpoint-safe sampling.
